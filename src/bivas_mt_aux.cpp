//
//  bivas_mt_aux.cpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/19.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include "bivas_mt_aux.hpp"
// #include "bivas_aux.hpp"
using namespace std;
using namespace arma;


void unlist(Rcpp::List& Xs, mat *pX[]) {
  int l = Xs.size();
  for (int i = 0; i < l; i++){
    SEXP xi = Xs[i];
    Rcpp::NumericMatrix Xr(xi);

    int n = Xr.nrow();
    int k = Xr.ncol();

    pX[i] = new mat(Xr.begin(),n,k,true);
  }
}

void unlist(Rcpp::List& ys, vec *py[]) {
  int l = ys.size();
  for (int i = 0; i < l; i++){
    SEXP yi = ys[i];
    Rcpp::NumericMatrix yr(yi);

    int n = yr.size();

    py[i] = new vec(yr.begin(),n,true);
  }
}

void centerize_mt(mat *pX[], vec *py[], mat& Xm, vec& ym){
  uword l = ym.size();
  for (int i = 0; i < l; i++){
    Xm.col(i) = mean(*pX[i],0).t();
    ym[i] = as_scalar(mean(*py[i]));

    (*pX[i]).each_row() -= Xm.col(i).t();
    (*py[i]) -= ym[i];
  }
}

void getInv_mt(mat *pX[], vec *py[], mat *pZ[], mat *pSZX[], vec *pSZy[], uword l, uword q, uword K){
  for (int i = 0; i < l; i++){
    mat Z = (*pZ[i]);
    mat invZZ = (Z.t() * Z).i();

    vec SZy = invZZ * Z.t() * (*py[i]);
    mat SZX = invZZ * Z.t() * (*pX[i]);
    pSZy[i] = new vec(SZy.begin(), (q+1));
    pSZX[i] = new mat(SZX.begin(), (q+1), K);

  }
}

void xtx_mt(mat& xtx, mat *pX[]){
  uword l = xtx.n_cols;
  for (int i = 0; i < l; i++){
    xtx.col(i) = sum(square(*pX[i]), 0).t();  //a K*l matrix
  }
}

void VB_update(mat *pX[], vec *py_bar[], mat &xtx, double logodds, vec &sb2, vec &se2, double &alpha,
               mat &mu, mat &s2, mat &alpha_jk, vec &pi_k, vec *pyt_j[], vec nn, double K){

  uword l = xtx.n_cols;
  s2 = se2.t() / (xtx.each_row() + (se2 / sb2).t()).each_row();

  // for (int j = 0; j < l; j++){
  //   int n = nn(j);
  //   vec y_bar = vec((*py_bar[j]).memptr(), n, false);
  //
  //   mat X_j = mat((*pX[j]).memptr(), n, K, false);
  //   // vec yt_jk = (*pyt_j[j]) - X_j * (pi_k % alpha_jk.col(j) % mu.col(j));
  //
  //   // s2.col(j) = se2[j] / (xtx.col(j) + se2[j] / sb2[j]);
  //
  //   for (int k = 0; k < K; k++){
  //     // mat X_j = mat((*pX[j]).memptr(), n, K, false);
  //     vec yt_jk = (*pyt_j[j]) - X_j.col(k) * (pi_k[k] * alpha_jk(k,j) * mu(k,j));
  //
  //     // s2(k,j) = se2[j] / (xtx(k,j) + se2[j] / sb2[j]);
  //
  //     mu(k,j) = s2(k,j) * sum(X_j.col(k) % (y_bar - yt_jk)) / se2[j];
  //
  //     double v = log(alpha / (1-alpha)) + .5 * pi_k[k] * log(s2(k,j) / sb2[j]) + .5 * pi_k[k] * mu(k,j) * mu(k,j) / s2(k,j);
  //     alpha_jk(k,j) = 1 / (1 + exp(-v));
  //
  //     // cout << accu(alpha_jk) << endl;
  //
  //     double u = logodds + .5 * sum(alpha_jk.row(k) % (log(s2.row(k) / sb2.t()) + square(mu.row(k)) / s2.row(k)));
  //     pi_k[k] = 1 / (1 + exp(-u));
  //
  //     *pyt_j[j] = yt_jk + X_j.col(k) * (pi_k[k] * alpha_jk(k,j) * mu(k,j));
  //   }
  //   // *pyt_j[j] = yt_jk + X_j * (pi_k % alpha_jk.col(j) % mu.col(j));
  // }

  vec *pyt_jk[l];

  for (int k = 0; k < K; k++){

    for (int j = 0; j < l; j++){
      mat X_j = mat((*pX[j]).memptr(), nn(j), K, false);

      vec yt_jk = (*pyt_j[j]) - X_j.col(k) * (pi_k[k] * alpha_jk(k,j) * mu(k,j));
      pyt_jk[j] = new vec(yt_jk.begin(),nn(j));

      mu(k,j) = s2(k,j) * sum(X_j.col(k) % (*py_bar[j] - *pyt_jk[j])) / se2[j];

      double v = log(alpha / (1-alpha)) + .5 * pi_k[k] * log(s2(k,j) / sb2[j]) + .5 * pi_k[k] * mu(k,j) * mu(k,j) / s2(k,j);
      alpha_jk(k,j) = 1 / (1 + exp(-v));
    }

    double u = logodds + .5 * sum(alpha_jk.row(k) % (log(s2.row(k) / sb2.t()) + square(mu.row(k)) / s2.row(k)));
    pi_k[k] = 1 / (1 + exp(-u));

    for (int j = 0; j < l; j++){
      mat X_j = mat((*pX[j]).memptr(), nn(j), K, false);
      *pyt_j[j] = *pyt_jk[j] + X_j.col(k) * (pi_k[k] * alpha_jk(k,j) * mu(k,j));

      delete pyt_jk[j];
    }
  }
}

void param_update(vec *py_bar[], mat &xtx, vec &sb2, vec &se2, double &alpha, mat &mu, mat &s2,
                  mat &alpha_jk, vec &pi_k, vec *pyt_j[], vec nn, double K){

  uword l = xtx.n_cols;
  for (int j = 0; j < l; j++){
    int n = nn(j);
    double sumyytilde = sum(square(*py_bar[j] - *pyt_j[j]));
    vec term1 = pi_k % alpha_jk.col(j) % (s2.col(j) + square(mu.col(j)));
    double term2 = sum((term1 - square(pi_k % alpha_jk.col(j) % mu.col(j))) % xtx.col(j));

    se2(j) = (sumyytilde + term2) / n;

    sb2(j) = sum(term1) / sum(pi_k % alpha_jk.col(j));

    // cout << sumyytilde << endl;
    // cout << sum(pi_k) << endl;
    // cout << sum(alpha_jk.col(j)) << endl;
    // cout << sum(mu.col(j)) << endl;
    // cout << sum(s2.col(j)) << endl;
    // cout << sum(pi_k % alpha_jk.col(j)) << endl;
    // cout << se2[j] << endl;
    // cout << sb2[j] << endl;
  }
  alpha = accu(alpha_jk) / (l*K);
  // cout << alpha << endl;
}

void cov_update(vec *py[], mat *pZ[], mat &cov, vec *py_bar[], vec *pSZy[], mat *pSZX[],
                vec &pi_k, mat &alpha_jk, mat &mu){
  uword l = alpha_jk.n_cols;
  for (int j = 0; j < l; j++){
    cov.col(j) = *pSZy[j] - *pSZX[j] * (pi_k % alpha_jk.col(j) % mu.col(j));
    *py_bar[j] = *py[j] - (*pZ[j] * cov.col(j));
  }
}

double lb_linear(vec nn, vec &se2){

  double lb = - 0.5 * sum(nn % log(se2)) - sum(nn);
  return lb;
}

double lb_linear(vec *py_bar[], vec *pyt_j[], mat& xtx, mat& alpha_jk, mat& mu, vec& pi_k,
                 vec &se2, mat &s2, vec nn){
  double l = mu.n_cols;
  double lb = 0;
  for(int j = 0; j < l; j++){
    double lb1 = - 0.5 * nn[j] * log(se2[j]);
    double lb2 = - 0.5 * sum(square(*py_bar[j]-*pyt_j[j])) / se2[j];
    double lb3 = - 0.5 * sum((alpha_jk.col(j) % pi_k % (s2.col(j) + square(mu.col(j))) - square(alpha_jk.col(j) % mu.col(j) % pi_k)) % xtx.col(j))/se2[j];
    lb = lb + lb1 + lb2 + lb3;
  }

  return lb;
}

double lb_gamma(double alpha, mat& alpha_jk){
  double lb = log(alpha) * accu(alpha_jk) + log(1 - alpha) * accu(1 - alpha_jk);
  return lb;
}

double lb_klbeta(vec &sb2, mat &s2, vec &pi_k, mat &alpha_jk, mat &mu){
  double K = pi_k.n_elem;
  double l = mu.n_cols;
  mat s2mu2 = s2+square(mu);
  double lb1 = - accu(alpha_jk % log(alpha_jk + (alpha_jk == 0)) + (1 - alpha_jk) % log(1 - alpha_jk + (alpha_jk == 1)));
  double lb2 = - sum(pi_k % log(pi_k + (pi_k == 0)) + (1 - pi_k) % log(1 - pi_k + (pi_k == 1)));
  double lb3 = 0.5 * accu(alpha_jk.each_col() % pi_k % ((1+log(s2.each_row() / sb2.t()))-s2mu2.each_row()/sb2.t()));

  return lb1+lb2+lb3;
}
