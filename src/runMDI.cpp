// runMDI.cpp
// =============================================================================
# include "runMDI.h"

using namespace Rcpp ;
using namespace arma ;

// =============================================================================
// runMDI function implementation
Rcpp::List runMDI(
    arma::uword R,
    arma::uword thin,
    arma::field< arma::mat > Y,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::umat labels,
    arma::umat fixed,
    arma::field< arma::vec > proposal_windows,
    bool save_parameters,
    bool save_imputed,
    arma::vec prior,
    arma::vec density_prior,
    arma::uvec save_allocation_probabilities,
    bool save_pointwise,
    bool phi_slice,
    arma::vec betas,
    arma::uword swap_scheme,
    arma::uword swap_every,
    arma::uword split_merge
) {
  
  if(thin < 1) {
    Rcpp::stop("thin must be a positive integer.");
  }
  
  // Parallel tempering ladder. T = 1 (no ladder) is the ordinary sampler.
  const uword T = std::max<uword>(betas.n_elem, 1);
  if(betas.n_elem == 0) {
    betas = arma::vec({1.0});
  }
  if(!betas.is_finite() || any(betas < 0.0) || any(betas > 1.0)) {
    Rcpp::stop("Inverse temperatures must lie in [0, 1].");
  }
  if(T > 1) {
    if(betas(T - 1) != 1.0) {
      Rcpp::stop("The last (coldest) inverse temperature must be 1, the posterior.");
    }
    for(uword t = 1; t < T; t++) {
      if(!(betas(t) > betas(t - 1))) {
        Rcpp::stop("Inverse temperatures must be strictly increasing.");
      }
    }
  }
  if(swap_scheme > 1) {
    Rcpp::stop("swap_scheme must be 0 (deterministic even-odd) or 1 (stochastic even-odd).");
  }
  if(swap_every < 1) {
    Rcpp::stop("swap_every must be a positive integer.");
  }
  
  const uword L = Y.n_elem, n_saved = R / thin + 1;
  uword save_ind = 0;
  
  if(save_allocation_probabilities.n_elem == 1) {
    save_allocation_probabilities = arma::uvec(L, arma::fill::value(save_allocation_probabilities(0)));
  }
  if(save_allocation_probabilities.n_elem != L) {
    Rcpp::stop("save_allocation_probabilities must have one entry per view.");
  }
  
  // One replica per temperature. Replica t starts at temperature t; after 
  // exchanges at_temp(t) is the replica currently at temperature t, so the 
  // coldest, at_temp(T - 1), is the one that is recorded.
  std::vector< std::unique_ptr<mdi> > replicas;
  replicas.reserve(T);
  for(uword t = 0; t < T; t++) {
    replicas.push_back(std::unique_ptr<mdi>(
      new mdi(Y, mixture_types, outlier_types, K, labels, fixed, prior, density_prior)
    ));
    replicas[t]->phi_slice = phi_slice;
    for(uword l = 0; l < L; l++) {
      // Only Gaussian process views use proposal windows
      if(mixture_types[l] == 3) {
        replicas[t]->mixtures[l]->density_ptr->receiveHyperParametersProposalWindows(proposal_windows[l]);
      }
    }
    if(betas(t) != 1.0) {
      replicas[t]->setBeta(betas(t));
    }
    replicas[t]->setSplitMerge(split_merge);
  }
  arma::uvec at_temp = arma::regspace<arma::uvec>(0, T - 1);
  auto cold = [&]() -> mdi& { return *replicas[at_temp(T - 1)]; };
  mdi& my_mdi_initial = cold();
  const uword N = my_mdi_initial.N;
  const uword LC2 = my_mdi_initial.LC2, K_max = my_mdi_initial.K_max;
  
  // Exchange diagnostics (T > 1 only)
  arma::uvec swap_attempts(T > 1 ? T - 1 : 0, arma::fill::zeros), 
    swap_accepts(T > 1 ? T - 1 : 0, arma::fill::zeros);
  arma::vec swap_accept_prob_sum(T > 1 ? T - 1 : 0, arma::fill::zeros);
  arma::uvec last_extreme(T, arma::fill::value(2));   // 2: not yet at the hottest
  arma::uvec round_trips(T, arma::fill::zeros);
  uword swap_round = 0;
  const uword n_pt_saved = (T > 1) ? R / thin + 1 : 0;
  arma::mat pt_data_ll(n_pt_saved, T, arma::fill::zeros);
  arma::umat pt_replica(n_pt_saved, T, arma::fill::zeros);
  if(T > 1) {
    last_extreme(at_temp(0)) = 0;
  }
  
  vec complete_likelihood_record(n_saved, arma::fill::zeros), 
    observed_likelihood_record(n_saved, arma::fill::zeros),
    joint_likelihood_record(n_saved, arma::fill::zeros),
    evidence(n_saved, arma::fill::zeros);
  mat pointwise_record;
  if(save_pointwise) {
    pointwise_record.zeros(n_saved, N);
  }
  
  mat phis_record(n_saved, LC2, arma::fill::zeros), 
    mass_record(n_saved, L, arma::fill::zeros),
    outlier_weight_record(n_saved, L, arma::fill::zeros);
  
  ucube class_record(n_saved, N, L), outlier_record(n_saved, N, L);
  class_record.zeros();
  outlier_record.zeros();
  
  cube weight_record(n_saved, K_max, L, arma::fill::zeros);
  ucube N_k_record(K_max, L, n_saved, arma::fill::zeros);
  
  field< cube > alloc(L);
  field< mat > hyper_record(L), parameter_record(L), imputed_record(L), pooled_record(L);
  field< umat > missing_cells(L);
  field< vec > acceptance_count(L);
  
  for(uword l = 0; l < L; l++) {
    // Only needed for semi-supervised views; N x K x draws is large
    alloc(l) = save_allocation_probabilities(l) ? zeros<cube>(N, K(l), n_saved) : zeros<cube>(0, 0, 0);
    hyper_record(l) = zeros< mat >(n_saved, 3 * K(l));
    pooled_record(l) = zeros< mat >(n_saved, my_mdi_initial.mixtures[l]->density_ptr->pooledHyperparameters().n_elem);
    acceptance_count(l) = zeros< vec >(3 * K(l));
    
    if(save_parameters) {
      parameter_record(l) = zeros< mat >(n_saved, my_mdi_initial.mixtures[l]->density_ptr->parameters().n_elem);
    }
    
    // The positions of missing entries (0-based row, column); fixed for the run
    const umat& has_missing = my_mdi_initial.mixtures[l]->density_ptr->has_missing;
    uvec missing_linear = find(has_missing == 1);
    missing_cells(l) = zeros<umat>(missing_linear.n_elem, 2);
    for(uword i = 0; i < missing_linear.n_elem; i++) {
      missing_cells(l)(i, 0) = missing_linear(i) % N;
      missing_cells(l)(i, 1) = missing_linear(i) / N;
    }
    if(save_imputed) {
      imputed_record(l) = zeros< mat >(n_saved, missing_linear.n_elem);
    }
  }
  
  auto record = [&](uword s) {
    // The coldest replica holds the posterior draw
    mdi& my_mdi = cold();
    for(uword l = 0; l < L; l++) {
      class_record.slice(l).row(s) = my_mdi.labels.col(l).t();
      weight_record.slice(l).row(s) = my_mdi.w.col(l).t();
      if(save_allocation_probabilities(l)) {
        alloc(l).slice(s) = my_mdi.mixtures[l]->alloc;
      }
      outlier_record.slice(l).row(s) = my_mdi.mixtures[l]->outliers.t();
      outlier_weight_record(s, l) = my_mdi.mixtures[l]->outlierComponent_ptr->outlier_weight;
      
      const auto& density_ptr = my_mdi.mixtures[l]->density_ptr;
      if(mixture_types(l) == 3) {
        hyper_record(l).row(s) = density_ptr->hypers.t();
      }
      if(pooled_record(l).n_cols > 0) {
        pooled_record(l).row(s) = density_ptr->pooledHyperparameters().t();
      }
      if(save_parameters) {
        parameter_record(l).row(s) = density_ptr->parameters().t();
      }
      if(save_imputed) {
        for(uword i = 0; i < missing_cells(l).n_rows; i++) {
          imputed_record(l)(s, i) = density_ptr->X(missing_cells(l)(i, 0), missing_cells(l)(i, 1));
        }
      }
    }
    const arma::vec pointwise = my_mdi.pointwiseLogLikelihood();
    joint_likelihood_record(s) = accu(pointwise);
    if(save_pointwise) {
      pointwise_record.row(s) = pointwise.t();
    }
    complete_likelihood_record(s) = my_mdi.complete_likelihood;
    observed_likelihood_record(s) = my_mdi.observed_likelihood;
    evidence(s) = my_mdi.Z_start;
    mass_record.row(s) = my_mdi.mass.t();
    phis_record.row(s) = my_mdi.phis.t();
    N_k_record.slice(s) = my_mdi.N_k;
  };
  
  // The data log-likelihood of every replica, the quantity the exchange
  // acceptance depends on, and the record of who sits where
  arma::vec replica_ll(T, arma::fill::zeros);
  auto record_pt = [&](uword s) {
    for(uword t = 0; t < T; t++) {
      pt_data_ll(s, t) = replicas[at_temp(t)]->dataLogLikelihood();
      pt_replica(s, t) = at_temp(t);
    }
  };
  
  // One round of exchange attempts between neighbouring temperatures. Pairs
  // within a round are disjoint, so the likelihoods are not stale between 
  // attempts. The acceptance probability of exchanging the states at 
  // temperatures i and j is 
  //   min(1, exp((beta_i - beta_j) (log L(state_j) - log L(state_i)))),
  // since the priors cancel in the ratio of the products of tempered targets.
  auto exchange = [&]() {
    for(uword t = 0; t < T; t++) {
      replica_ll(t) = replicas[t]->dataLogLikelihood();
    }
    const uword parity = (swap_scheme == 0) ? (swap_round % 2) : (randu() < 0.5 ? 1 : 0);
    for(uword t = parity; t + 1 < T; t += 2) {
      const uword a = at_temp(t), b = at_temp(t + 1);
      const double log_alpha = (betas(t) - betas(t + 1)) * (replica_ll(b) - replica_ll(a));
      swap_attempts(t)++;
      swap_accept_prob_sum(t) += std::min(1.0, std::exp(log_alpha));
      if(std::log(randu()) < log_alpha) {
        swap_accepts(t)++;
        at_temp(t) = b;
        at_temp(t + 1) = a;
        replicas[b]->setBeta(betas(t));
        replicas[a]->setBeta(betas(t + 1));
      }
    }
    swap_round++;
    
    // Round trips: hottest -> coldest -> hottest
    for(uword t = 0; t < T; t++) {
      const uword rep = at_temp(t);
      if(t == 0) {
        if(last_extreme(rep) == 1) {
          round_trips(rep)++;
        }
        last_extreme(rep) = 0;
      } else if(t == T - 1 && last_extreme(rep) == 0) {
        last_extreme(rep) = 1;
      }
    }
  };
  
  // The initial state
  for(uword t = 0; t < T; t++) {
    replicas[t]->complete_likelihood = accu(replicas[t]->complete_likelihood_vec);
    replicas[t]->observed_likelihood = accu(replicas[t]->observed_likelihood_vec);
  }
  record(save_ind);
  if(T > 1) {
    record_pt(save_ind);
  }
  
  for(uword r = 0; r < R; r++) {
    Rcpp::checkUserInterrupt();
    for(uword t = 0; t < T; t++) {
      replicas[t]->sweep(r);
    }
    if(T > 1 && (r + 1) % swap_every == 0) {
      exchange();
    }
    if((r + 1) % thin == 0) {
      save_ind++;
      record(save_ind);
      if(T > 1) {
        record_pt(save_ind);
      }
    }
  }
  
  mdi& my_mdi = cold();
  for(uword l = 0; l < L; l++) {
    if(mixture_types(l) == 3) {
      // Accepted proposals over proposals made, per hyperparameter and component
      // (a component is only updated while it is occupied)
      const auto& d = my_mdi.mixtures[l]->density_ptr;
      acceptance_count(l) = conv_to< vec >::from(d->acceptance_count)
        / arma::max(arma::ones<vec>(d->acceptance_attempts.n_elem), conv_to< vec >::from(d->acceptance_attempts));
    }
  }
  
  // Mass acceptance rate averaged over replicas (a replica visits many
  // temperatures, so a per-temperature rate is not recorded)
  arma::vec mass_acceptance_total(L, arma::fill::zeros);
  for(uword t = 0; t < T; t++) {
    mass_acceptance_total += replicas[t]->mass_acceptance_count;
  }
  
  Rcpp::List pt = Rcpp::List::create();
  if(T > 1) {
    arma::vec rejection_rate(T - 1);
    for(uword t = 0; t + 1 < T; t++) {
      rejection_rate(t) = (swap_attempts(t) > 0) ? 
        1.0 - swap_accept_prob_sum(t) / (double) swap_attempts(t) : arma::datum::nan;
    }
    pt = Rcpp::List::create(
      Named("betas") = betas,
      Named("swap_scheme") = swap_scheme,
      Named("swap_every") = swap_every,
      Named("swap_attempts") = swap_attempts,
      Named("swap_accepts") = swap_accepts,
      Named("swap_accept_prob_sum") = swap_accept_prob_sum,
      Named("rejection_rate") = rejection_rate,
      Named("round_trips") = round_trips,
      Named("data_log_likelihood") = pt_data_ll,
      Named("replica") = pt_replica + 1
    );
  }
  
  Rcpp::List result = List::create(
      Named("allocations") = class_record,
      Named("phis") = phis_record,
      Named("weights") = weight_record,
      Named("mass") = mass_record,
      Named("outliers") = outlier_record,
      Named("outlier_weights") = outlier_weight_record,
      Named("allocation_probabilities") = alloc,
      Named("N_k") = N_k_record,
      Named("complete_likelihood") = complete_likelihood_record,
      Named("observed_likelihood") = observed_likelihood_record,
      Named("joint_likelihood") = joint_likelihood_record,
      Named("pointwise_likelihood") = pointwise_record,
      Named("evidence") = evidence,
      Named("hypers") = hyper_record,
      Named("acceptance_count") = acceptance_count,
      Named("mass_acceptance_rate") = mass_acceptance_total / ((double) std::max<uword>(R, 1) * (double) T),
      Named("parameters") = parameter_record,
      Named("pooled_hyperparameters") = pooled_record,
      Named("imputed") = imputed_record,
      Named("missing_cells") = missing_cells
  );
  result["tempering"] = pt;
  uword sm_attempts = 0, sm_accepts = 0;
  for(uword t = 0; t < T; t++) {
    sm_attempts += replicas[t]->split_merge_attempts;
    sm_accepts += replicas[t]->split_merge_accepts;
  }
  result["split_merge"] = Rcpp::List::create(
    Named("attempts") = (double) sm_attempts,
    Named("accepts") = (double) sm_accepts
  );
  return result;
}
