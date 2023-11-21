import numpy as np
from scipy.linalg import inv, det, pinv
from scipy.stats import multivariate_normal

def ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S, tau, p_S, det_S):
    result = p_S / 2.0 * np.log(tau) - 0.5 * np.log(det_S) + zSigmaz_S / 2.0
    return result

def ind_Normal_fixed_sigma_marginal_external(index_vec_input, Sigma, z, tau, p_S, y_sigma):
    index_vec = index_vec_input - 1
    Sigma_S = Sigma[np.ix_(index_vec, index_vec)]
    A = tau * np.eye(p_S)
    
    det_S = det(Sigma_S + A)
    Sigma_S_inv = inv(Sigma_S + A)
    
    sub_z = z[index_vec]
    zSigmaz_S = np.dot(sub_z.T, np.dot(Sigma_S_inv, sub_z))
    
    b = ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S, tau, p_S, det_S)
    
    results = b
    
    return results

def outlier_ind_Normal_marginal_external(index_vec_input, Sigma, z, tau, p_S, y_sigma):
    index_vec = index_vec_input - 1
    Sigma_S = Sigma[np.ix_(index_vec, index_vec)]
    A = tau * np.eye(p_S)
    Sigma_S_I_inv = pinv(Sigma_S + A, rcond=0.00001)
    Sigma_S_inv = pinv(Sigma_S, rcond=0.00001)
    det_S = np.abs(det(Sigma_S_inv))
    det_I_S = np.abs(det(Sigma_S_I_inv))
    sub_z = z[index_vec]
    zSigmaz_S = np.dot(sub_z.T, np.dot(Sigma_S_inv, sub_z))
    zSigmaz_I_S = np.dot(sub_z.T, np.dot(Sigma_S_I_inv, sub_z))
    b = 0.5 * (np.log(det_S) + np.log(det_I_S)) - 0.5 * (zSigmaz_S[0, 0] - zSigmaz_I_S[0, 0])
    return b


def outlier_ind_Normal_marginal_external(index_vec_input, Sigma, z, tau, p_S, y_sigma):
    index_vec = index_vec_input - 1
    
    Sigma_S = Sigma[np.ix_(index_vec, index_vec)]
    A = tau * np.eye(p_S)
    
    Sigma_S_I_inv = pinv(Sigma_S + A, rcond=0.00001)
    Sigma_S_inv = pinv(Sigma_S, rcond=0.00001)
    
    det_S = np.abs(det(Sigma_S_inv))
    det_I_S = np.abs(det(Sigma_S_I_inv))
    
    sub_z = z[index_vec]
    zSigmaz_S = np.dot(sub_z, np.dot(Sigma_S_inv, sub_z))
    zSigmaz_I_S = np.dot(sub_z, np.dot(Sigma_S_I_inv, sub_z))
    
    b = 0.5 * (np.log(det_S) + np.log(det_I_S)) - 0.5 * (zSigmaz_S - zSigmaz_I_S)
    results = b
    
    return results