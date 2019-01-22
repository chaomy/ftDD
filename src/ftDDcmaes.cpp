/*
 * @Author: chaomy
 * @Date:   2018-04-14 13:19:01
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-07-11 15:49:32
 */

#include "ftDD.hpp"

using arma::accu;
using arma::endr;
using arma::linspace;
using arma::mat;
using arma::randu;
using arma::vec;

double ftDD::cmaes(arma::mat& iterate, const vector<string>& kk) {
  int maxIt = 10000;
  double tolerance = 1e-20;

  // Population size.
  // int lambda = (4 + std::round(3 * std::log(iterate.n_elem))) * 10;
  int lambda = (4 + std::round(3 * std::log(iterate.n_elem))) * 30;

  // Parent weights.
  const size_t mu = std::round(lambda / 2);

  arma::vec w = std::log(mu + 0.5) -
                arma::log(arma::linspace<arma::vec>(0, mu - 1, mu) + 1.0);
  w /= arma::sum(w);

  // Number of effective solutions.
  const double muEffective = 1 / arma::accu(arma::pow(w, 2));

  // Step size control parameters.
  arma::vec sigma(3);

  // double upperBound = 10., lowerBound = -10.;
  // sigma(0) = 0.1 * (upperBound - lowerBound);
  sigma(0) = 2.0;

  const double cs = (muEffective + 2) / (iterate.n_elem + muEffective + 5);
  const double ds =
      1 + cs +
      2 * std::max(std::sqrt((muEffective - 1) / (iterate.n_elem + 1)) - 1,
                   0.0);
  const double enn =
      std::sqrt(iterate.n_elem) * (1.0 - 1.0 / (4.0 * iterate.n_elem) +
                                   1.0 / (21 * std::pow(iterate.n_elem, 2)));

  // Covariance update parameters Cumulation for distribution.
  const double cc = (4 + muEffective / iterate.n_elem) /
                    (4 + iterate.n_elem + 2 * muEffective / iterate.n_elem);
  const double h = (1.4 + 2.0 / (iterate.n_elem + 1.0)) * enn;

  const double c1 = 2 / (std::pow(iterate.n_elem + 1.3, 2) + muEffective);
  const double alphaMu = 2;
  const double cmu = std::min(
      1 - c1,
      alphaMu * (muEffective - 2 + 1 / muEffective) /
          (std::pow(iterate.n_elem + 2, 2) + alphaMu * muEffective / 2));

  arma::cube mps(iterate.n_rows, iterate.n_cols, 3);  // meam

  // mps.slice(0) =
  //     lowerBound +
  //     arma::randu(iterate.n_rows, iterate.n_cols) * (upperBound -
  //     lowerBound);

  mps.slice(0) = iterate;

  arma::mat step = arma::zeros(iterate.n_rows, iterate.n_cols);

  // Calculate the first objective function.
  double currentobj =
      (this->*calobj[sparams["ptype"]])(decodestdv(mps.slice(0)), kk);
  double overallobj = currentobj;
  double lastobj = 1e30;

  // Population parameters.
  arma::cube pStep(iterate.n_rows, iterate.n_cols, lambda);
  arma::cube pps(iterate.n_rows, iterate.n_cols, lambda);
  arma::vec pobj(lambda);
  arma::cube ps = arma::zeros(iterate.n_rows, iterate.n_cols, 2);
  arma::cube pc = ps;
  arma::cube C(iterate.n_elem, iterate.n_elem, 2);
  C.slice(0).eye();

  // Covariance matrix parameters.
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec eigvalZero = arma::zeros(iterate.n_elem);

  // The current visitation order (sorted by population objectives).
  arma::uvec idx = arma::linspace<arma::uvec>(0, lambda - 1, lambda);

  for (size_t i = 1; i < maxIt; ++i) {
    const size_t idx0 = (i - 1) % 2;
    const size_t idx1 = i % 2;

    arma::mat covLower;
    if (!arma::chol(covLower, C.slice(idx0), "lower")) break;

    for (size_t j = 0; j < lambda; ++j) {
      if (iterate.n_rows > iterate.n_cols) {
        pStep.slice(idx(j)) =
            covLower * arma::randn(iterate.n_rows, iterate.n_cols);
      } else {
        pStep.slice(idx(j)) =
            arma::randn(iterate.n_rows, iterate.n_cols) * covLower;
      }
      pps.slice(idx(j)) = mps.slice(idx0) + sigma(idx0) * pStep.slice(idx(j));
      pobj(idx(j)) =
          (this->*calobj[sparams["ptype"]])(decodestdv(pps.slice(idx(j))), kk);
    }

    // Sort population.
    idx = sort_index(pobj);

    step = w(0) * pStep.slice(idx(0));
    for (size_t j = 1; j < mu; ++j) step += w(j) * pStep.slice(idx(j));

    mps.slice(idx1) = mps.slice(idx0) + sigma(idx0) * step;

    currentobj =
        (this->*calobj[sparams["ptype"]])(decodestdv(mps.slice(idx1)), kk);

    // for meams
    vector<double> tm(move(decodestdv(mps.slice(0))));
    cout << setprecision(15) << "CMA-ES: i = " << i << ", objective "
         << overallobj << " "
         << " cs " << sigma(idx1) << " " << (lastobj - overallobj) / lastobj
         << " " << tm[0] << endl;

    if (currentobj < overallobj) {
      overallobj = currentobj;
      iterate = mps.slice(idx1);
    }

    // Update Step Size.
    if (iterate.n_rows > iterate.n_cols) {
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * covLower.t() * step;
    } else {
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * step * covLower.t();
    }

    const double psNorm = arma::norm(ps.slice(idx1));
    sigma(idx1) =
        sigma(idx0) * std::pow(std::exp(cs / ds * psNorm / enn - 1), 0.3);

    // Update covariance matrix.
    if ((psNorm / sqrt(1 - std::pow(1 - cs, 2 * i))) < h) {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0) +
                       std::sqrt(cc * (2 - cc) * muEffective) * step;

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t());
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1));
      }
    } else {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0);

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t() +
                              (cc * (2 - cc)) * C.slice(idx0));
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1) +
                              (cc * (2 - cc)) * C.slice(idx0));
      }
    }

    if (iterate.n_rows > iterate.n_cols) {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)) *
                                            pStep.slice(idx(j)).t();
      }
    } else {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)).t() *
                                            pStep.slice(idx(j));
      }
    }

    arma::eig_sym(eigval, eigvec, C.slice(idx1));
    const arma::uvec negativeEigval = find(eigval < 0, 1);
    if (!negativeEigval.is_empty()) {
      if (negativeEigval(0) == 0) {
        C.slice(idx1).zeros();
      } else {
        C.slice(idx1) = eigvec.cols(0, negativeEigval(0) - 1) *
                        arma::diagmat(eigval.subvec(0, negativeEigval(0) - 1)) *
                        eigvec.cols(0, negativeEigval(0) - 1).t();
      }
    }

    if (std::isnan(overallobj) || std::isinf(overallobj)) {
      cout << "CMA-ES: converged to " << overallobj << "; "
           << "terminating with failure.  Try a smaller step size?"
           << std::endl;
      return overallobj;
    }

    if (std::abs(lastobj - overallobj) < tolerance && i > 200) {
      cout << "CMA-ES: minimized within tolerance " << tolerance << "; "
           << "terminating optimization." << std::endl;
      return overallobj;
    }

    lastobj = overallobj;
  }

  return overallobj;
}