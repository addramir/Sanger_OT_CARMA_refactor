#include <pybind11/pybind11.h>
#include <armadillo>

namespace py = pybind11;

double ind_normal_sigma_fixed_marginal_fun_indi(const double &zSigmaz_S, const double &tau,
                         const double & p_S, const double & det_S){
  double result = p_S / 2.00 * log(tau) - 0.50 * log(det_S) + (zSigmaz_S / 2.00);
  return result;
}

double ind_Normal_fixed_sigma_marginal_external(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double &tau,  const double & p_S, const double & y_sigma){
  arma::uvec index_vec = index_vec_input - 1;
  arma::mat Sigma_S = Sigma.submat(index_vec, index_vec);
  arma::mat A = tau * arma::eye(p_S, p_S);
  double det_S = det(Sigma_S + A);
  arma::mat Sigma_S_inv = pinv(Sigma_S + A, 0.00001, "std");
  arma::mat sub_z = z.rows(index_vec);
  arma::mat zSigmaz_S = sub_z.t() * Sigma_S_inv * sub_z;
  double b = ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S(0, 0), tau, p_S, det_S);
  return b;
}

double outlier_ind_Normal_marginal_external(const arma::uvec & index_vec_input, const arma::mat & Sigma,
                           const arma::vec & z, const double &tau,  const double & p_S, const double & y_sigma){
  arma::uvec index_vec = index_vec_input - 1;
  arma::mat Sigma_S = Sigma.submat(index_vec, index_vec);
  arma::mat A = tau * arma::eye(p_S, p_S);
  arma::mat Sigma_S_I_inv = pinv(Sigma_S + A, 0.00001, "std");
  arma::mat Sigma_S_inv = pinv(Sigma_S, 0.00001, "std");
  double det_S = det(Sigma_S_inv);
  double det_I_S = det(Sigma_S_I_inv);
  det_I_S = fabs(det_I_S);
  det_S = fabs(det_S);
  arma::mat sub_z = z.rows(index_vec);
  arma::mat zSigmaz_S = sub_z.t() * Sigma_S_inv * sub_z;
  arma::mat zSigmaz_I_S = sub_z.t() * Sigma_S_I_inv * sub_z;
  double b = 0.5 * (log(det_S) + log(det_I_S)) - 0.5 * (zSigmaz_S(0, 0) - zSigmaz_I_S(0, 0));
  return b;
}

PYBIND11_MODULE(ind_normal_carma, m) {
    m.def("ind_normal_sigma_fixed_marginal_fun_indi", &ind_normal_sigma_fixed_marginal_fun_indi);
    m.def("ind_Normal_fixed_sigma_marginal_external", &ind_Normal_fixed_sigma_marginal_external);
    m.def("outlier_ind_Normal_marginal_external", &outlier_ind_Normal_marginal_external);
}
