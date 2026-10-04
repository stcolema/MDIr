// runSMC.cpp
// =============================================================================
# include "runSMC.h"

using namespace Rcpp ;
using namespace arma ;

namespace {

// Conditional ESS of the incremental weights exp(delta * ell) under the current
// normalised weights W (Zhou, Johansen and Aston, 2016, eq. 3.16), as a
// fraction of the number of particles
double conditionalESSFraction(const arma::vec& W, const arma::vec& ell, double delta) {
  const arma::vec d = delta * ell;
  const double m = d.max();
  const arma::vec e = arma::exp(d - m);
  const double a = arma::accu(W % e), b = arma::accu(W % arma::square(e));
  return (a * a) / b;
}

arma::uvec resampleIndices(const arma::vec& W, uword scheme) {
  const uword N = W.n_elem;
  arma::uvec anc(N);
  const arma::vec cum = arma::cumsum(W);
  if(scheme == 0) {
    // systematic: one uniform, N equally spaced positions
    const double u = randu() / (double) N;
    uword i = 0;
    for(uword j = 0; j < N; j++) {
      const double pos = u + (double) j / (double) N;
      while(i + 1 < N && cum(i) < pos) {
        i++;
      }
      anc(j) = i;
    }
  } else {
    // multinomial: N independent draws
    for(uword j = 0; j < N; j++) {
      const double pos = randu() * cum(N - 1);
      uword i = 0;
      while(i + 1 < N && cum(i) < pos) {
        i++;
      }
      anc(j) = i;
    }
  }
  return anc;
}

}

Rcpp::List runMDISMC(
    arma::uword n_particles,
    arma::field<arma::mat> Y,
    arma::uvec K,
    arma::uvec mixture_types,
    arma::uvec outlier_types,
    arma::umat fixed,
    arma::vec prior,
    arma::vec density_prior,
    bool phi_slice,
    arma::vec betas,
    bool adaptive,
    double cess_target,
    double resample_threshold,
    arma::uword resample_scheme,
    arma::uword sweeps_per_step,
    arma::uword max_steps,
    arma::uword final_sweeps,
    arma::uword final_thin,
    double beta_start,
    arma::uword start_sweeps
) {
  if(n_particles < 2) {
    Rcpp::stop("At least two particles are needed.");
  }
  if(any(vectorise(fixed) != 0)) {
    Rcpp::stop("The sequential Monte Carlo sampler supports unsupervised views only "
               "(the prior of the labels cannot be sampled exactly given observed labels).");
  }
  if(sweeps_per_step < 1 || final_thin < 1 || max_steps < 1) {
    Rcpp::stop("sweeps_per_step, final_thin and max_steps must be positive.");
  }
  if(!(beta_start >= 0.0 && beta_start < 1.0)) {
    Rcpp::stop("beta_start must lie in [0, 1).");
  }
  if(!adaptive && betas(0) <= beta_start) {
    Rcpp::stop("A fixed schedule must start above beta_start.");
  }
  if(resample_scheme > 1) {
    Rcpp::stop("resample_scheme must be 0 (systematic) or 1 (multinomial).");
  }
  if(adaptive) {
    if(!(cess_target > 0.0 && cess_target < 1.0)) {
      Rcpp::stop("cess_target must lie in (0, 1).");
    }
  } else {
    if(betas.n_elem < 1 || !betas.is_finite() || any(betas <= 0.0) || any(betas > 1.0)
         || betas(betas.n_elem - 1) != 1.0) {
      Rcpp::stop("A fixed schedule needs inverse temperatures in (0, 1] ending at 1.");
    }
    for(uword t = 1; t < betas.n_elem; t++) {
      if(!(betas(t) > betas(t - 1))) {
        Rcpp::stop("The inverse temperatures must be strictly increasing.");
      }
    }
  }

  const uword L = Y.n_elem, N = Y(0).n_rows;
  const uword P = n_particles;
  arma::umat labels0(N, L, arma::fill::zeros);

  std::vector< std::unique_ptr<mdi> > particles;
  particles.reserve(P);
  for(uword i = 0; i < P; i++) {
    particles.push_back(std::unique_ptr<mdi>(
      new mdi(Y, mixture_types, outlier_types, K, labels0, fixed, prior, density_prior)
    ));
    particles[i]->phi_slice = phi_slice;
    particles[i]->initialiseFromPrior();
    // beta = 0 (this also refuses the models the tempering is not defined for)
    particles[i]->setBeta(0.0);
  }
  const uword K_max = particles[0]->K_max, LC2 = particles[0]->LC2;

  arma::vec W(P, arma::fill::value(1.0 / (double) P)), ell(P);
  arma::uvec root = arma::regspace<arma::uvec>(0, P - 1);

  double beta = 0.0, log_evidence = 0.0;
  uword step = 0, sweep_counter = 0;

  // A start above the prior: equilibrate every particle at beta_start. The
  // weights stay equal, and the evidence estimate is that of the path from
  // beta_start to 1.
  if(beta_start > 0.0) {
    for(uword i = 0; i < P; i++) {
      particles[i]->setBeta(beta_start);
    }
    for(uword s = 0; s < start_sweeps; s++) {
      Rcpp::checkUserInterrupt();
      for(uword i = 0; i < P; i++) {
        particles[i]->sweep(sweep_counter);
      }
      sweep_counter++;
    }
    beta = beta_start;
  }
  bool forced = false;

  std::vector<double> beta_trace, ess_trace, cess_trace, evidence_increment;
  std::vector<uword> resampled_trace, unique_root_trace;

  while(beta < 1.0) {
    Rcpp::checkUserInterrupt();
    for(uword i = 0; i < P; i++) {
      ell(i) = particles[i]->dataLogLikelihood();
    }

    // The next inverse temperature
    double beta_next = 1.0;
    if(adaptive) {
      if(step + 1 >= max_steps) {
        beta_next = 1.0;
        forced = true;
      } else if(conditionalESSFraction(W, ell, 1.0 - beta) >= cess_target) {
        beta_next = 1.0;
      } else {
        double lo = 0.0, hi = 1.0 - beta;
        for(uword it = 0; it < 80; it++) {
          const double mid = 0.5 * (lo + hi);
          if(conditionalESSFraction(W, ell, mid) >= cess_target) {
            lo = mid;
          } else {
            hi = mid;
          }
        }
        // a step that is too small to change the weights would never finish
        const double delta = std::max(lo, 1e-10);
        beta_next = std::min(1.0, beta + delta);
      }
    } else {
      beta_next = betas(step);
    }
    const double delta = beta_next - beta;

    // Incremental weights, evidence increment and the new normalised weights
    const arma::vec lw = delta * ell + arma::log(W);
    const double lse = logSumExp(lw);
    evidence_increment.push_back(lse);
    log_evidence += lse;
    const double cess = conditionalESSFraction(W, ell, delta);
    W = arma::exp(lw - lse);
    const double ess = 1.0 / arma::accu(arma::square(W));

    uword resampled = 0;
    if(ess < resample_threshold * (double) P) {
      const arma::uvec anc = resampleIndices(W, resample_scheme);
      std::vector<mdiState> snapshot;
      snapshot.reserve(P);
      for(uword i = 0; i < P; i++) {
        snapshot.push_back(particles[i]->saveState());
      }
      const arma::uvec old_root = root;
      for(uword j = 0; j < P; j++) {
        particles[j]->loadState(snapshot[anc(j)]);
        root(j) = old_root(anc(j));
      }
      W.fill(1.0 / (double) P);
      resampled = 1;
    }

    for(uword i = 0; i < P; i++) {
      particles[i]->setBeta(beta_next);
    }
    for(uword s = 0; s < sweeps_per_step; s++) {
      for(uword i = 0; i < P; i++) {
        particles[i]->sweep(sweep_counter);
      }
      sweep_counter++;
    }

    arma::uvec sorted_root = arma::sort(root);
    uword n_unique = 1;
    for(uword i = 1; i < P; i++) {
      if(sorted_root(i) != sorted_root(i - 1)) {
        n_unique++;
      }
    }
    beta_trace.push_back(beta_next);
    ess_trace.push_back(ess);
    cess_trace.push_back(cess);
    resampled_trace.push_back(resampled);
    unique_root_trace.push_back(n_unique);

    beta = beta_next;
    step++;
  }

  // Recorded draws: the state at the end of annealing and then thinned extra
  // sweeps at beta = 1. Each particle keeps its weight for all of them.
  const uword n_draws = 1 + final_sweeps / final_thin;
  arma::field<arma::ucube> allocations(P);
  arma::field<arma::cube> weights_record(P);
  arma::field<arma::mat> phis_record(P), mass_record(P);
  arma::mat data_ll(n_draws, P, arma::fill::zeros);
  for(uword i = 0; i < P; i++) {
    allocations(i) = arma::ucube(n_draws, N, L, arma::fill::zeros);
    weights_record(i) = arma::cube(n_draws, K_max, L, arma::fill::zeros);
    phis_record(i) = arma::mat(n_draws, LC2, arma::fill::zeros);
    mass_record(i) = arma::mat(n_draws, L, arma::fill::zeros);
  }
  auto record = [&](uword d) {
    for(uword i = 0; i < P; i++) {
      const mdi& m = *particles[i];
      for(uword l = 0; l < L; l++) {
        allocations(i).slice(l).row(d) = m.labels.col(l).t();
        weights_record(i).slice(l).row(d) = m.w.col(l).t();
      }
      phis_record(i).row(d) = m.phis.t();
      mass_record(i).row(d) = m.mass.t();
      data_ll(d, i) = particles[i]->dataLogLikelihood();
    }
  };
  record(0);
  uword saved = 0;
  for(uword s = 0; s < final_sweeps; s++) {
    Rcpp::checkUserInterrupt();
    for(uword i = 0; i < P; i++) {
      particles[i]->sweep(sweep_counter);
    }
    sweep_counter++;
    if((s + 1) % final_thin == 0) {
      saved++;
      record(saved);
    }
  }

  return Rcpp::List::create(
    Named("allocations") = allocations,
    Named("weights_record") = weights_record,
    Named("phis") = phis_record,
    Named("mass") = mass_record,
    Named("data_log_likelihood") = data_ll,
    Named("particle_weights") = W,
    Named("log_evidence") = log_evidence,
    Named("beta") = beta_trace,
    Named("ess") = ess_trace,
    Named("cess") = cess_trace,
    Named("evidence_increment") = evidence_increment,
    Named("resampled") = resampled_trace,
    Named("unique_roots") = unique_root_trace,
    Named("root") = root + 1,
    Named("forced_final_step") = forced
  );
}
