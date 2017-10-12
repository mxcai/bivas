//
//  bivas_aux.cpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/19.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include "bivas_aux.hpp"
using namespace std;
using namespace arma;


void VB_update(const mat& X, const vec& y_bar, const vec& xtx, const vec& group, const vec& glevel, double logodds,
               const double& sb2, const double& se2, const double& alpha, vec& mu, vec& s2, vec& alpha_jk,
               vec& pi_k, vec& pi_p, vec& yt, vec& offDiag, uword n, uword p, uword K){

    for(int k = 0; k < K; k++) {
        uvec idx = find(group == glevel(k));
        uword n_idx = idx.n_elem;

        vec alpha_k = alpha_jk(idx);
        vec mu_k = mu(idx);
        vec s2_k = s2(idx);
        mat X_k = X.cols(idx);
        //        mat X_k = conv_to<mat>::from(XX);
        vec xtx_k = xtx(idx);

        //        vec ytw = matvec(lpfX, idx, r, n, n_idx);
        vec ytw = X_k * (alpha_k % mu_k);

        vec ytb_k = yt - pi_k[k] * ytw;


        for(int j = 0; j < n_idx; j++) {

            //            uword idx_j = idx[j];
            //            vec X_j = getcol(lpfX, idx_j, n);
            vec X_j = X_k.col(j);

            vec ytw_k = ytw - alpha_k[j] * mu_k[j] * X_j;

            mu_k[j] = sum(X_j % (y_bar - ytw_k - ytb_k)) / (xtx_k[j] + se2 / sb2);


            double v = log(alpha / (1 - alpha)) + .5 * pi_k[k] * log(s2_k[j] / sb2) + .5 * pi_k[k] * mu_k[j] * mu_k[j] / s2_k[j];
            alpha_k[j] = 1 / (1 + exp(-v));

            ytw = ytw_k + alpha_k[j] * mu_k[j] * X_j;

        }

        alpha_jk(idx) = alpha_k;

        mu(idx) = mu_k;

        double u = logodds + .5 * sum(alpha_k % (log(s2_k / sb2) + square(mu_k) / s2_k));
        pi_k[k] = 1 / (1 + exp(-u));

        extend_pi_p(pi_k, pi_p, idx, n_idx, k);

        yt = ytb_k + pi_k[k] * ytw;


        mat Xc = (X_k).each_row() % (alpha_k.t() % mu_k.t());
        mat XcX = Xc.t() * Xc;
        offDiag[k] = (pi_k[k] - pow(pi_k[k], 2)) * (accu(XcX) - sum(XcX.diag()));

    }
}



void param_update(const vec& y_bar, const vec& xtx, double& sb2, double& se2, double& alpha, const vec& mu,
                  const vec& s2, const vec& alpha_jk, const vec& pi_p, const vec& yt, const vec& offDiag,
                  const vec& pi_a_s_mu, uword n, uword p){

    se2 = sum(square(y_bar - yt)) / n + sum((pi_a_s_mu - square(pi_p % alpha_jk % mu)) % xtx) / n + sum(offDiag) / n;

    sb2 = sum(pi_a_s_mu) / sum(pi_p % alpha_jk);

    alpha = sum(alpha_jk) / p;

}



double lb_linear(const vec& y_bar, vec& yt, const vec& xtx, vec& pi_p, vec& alpha_jk, vec& mu, vec& pi_a_s_mu,
                 double se2, vec& offDiag){
  uword n = y_bar.n_elem;

  double lb1 = - 0.5 * n * log(se2);
  double lb2 = - 0.5 * sum(square(y_bar-yt)) / se2;
  double lb3 = - 0.5 * sum((pi_a_s_mu - square(pi_p % alpha_jk % mu)) % xtx) / se2;
  double lb4 = - 0.5 * sum(offDiag) / se2;

  return lb1 + lb2 + lb3 + lb4;

}


double lb_gamma(double alpha, vec& alpha_jk){
    double lb = log(alpha) * sum(alpha_jk) + log(1 - alpha) * sum(1 - alpha_jk);
    return lb;
}

double lb_eta(double logodds, vec& pi_k){
    return sum((pi_k - 1) * logodds + (-logpexp(-logodds)));
}

double lb_klbeta(double sb2, double se2, vec& mu, vec& alpha_jk, vec& pi_k, vec& pi_p, vec& s2, vec& pi_a_s_mu){
    double lb1 = - sum(alpha_jk % log(alpha_jk + (alpha_jk == 0)) + (1 - alpha_jk) % log(1 - alpha_jk + (alpha_jk == 1)));
    double lb2 = - sum(pi_k % log(pi_k + (pi_k == 0)) + (1 - pi_k) % log(1 - pi_k + (pi_k == 1)));
    double lb3 = 0.5 * sum(pi_p % alpha_jk % (1 + log(s2 / sb2)));
    double lb4 = - 0.5 * sum (pi_a_s_mu + sb2) / sb2;

    return lb1 + lb2 + lb3 + lb4;
}

vec logpexp(vec& x) {
    vec y = x;
    uvec idx = (x < 16);
    y(idx) = log(1 + exp(x(idx)));
    return y;
}

double logpexp(double x) {
    double y = log(1 + exp(x));
    return y;
}

void getInv(const mat& X, const vec& y, const mat& Z, mat& SZX, vec& SZy, uword p) {
    mat invZZ = (Z.t() * Z).i();
    SZy = invZZ * Z.t() * y;
    SZX = invZZ * Z.t() * X;
}

void extend_pi_p(mat& pi_k, mat& pi_p, const vec& group, const vec& glevel, uword K){
    for(int i = 0; i < K; i++){
        uvec idx = find(group == glevel(i));
        pi_p.each_row(idx) = pi_k.row(i);
    }
}

void extend_pi_p(vec& pi_k, vec& pi_p, uvec idx, uword n_idx, int k) {
    for(int i = 0; i < n_idx; i++) {
        uword j = idx[i];
        pi_p[j] = pi_k[k];
    }
}


void ptr2List(Rcpp::List& Lst, vec *Lst_ptr[]) {
  int ns = Lst.size();
  for (int i = 0; i < ns; i++){
    Lst[i] = *Lst_ptr[i];
  }
}

