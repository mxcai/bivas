//
//  bivas_mt_aux.hpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/19.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#ifndef bivas_mt_aux_hpp
#define bivas_mt_aux_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
// #include <Rcpp.h>

using namespace arma;
using namespace std;



void unlist(Rcpp::List& Xs, mat *pX[]);
void unlist(Rcpp::List& ys, vec *py[]);
void centerize_mt(mat *pX[], vec *py[], mat& Xm, vec& ym);
void getInv_mt(mat *pX[], vec *py[], mat *pZ[], mat *pSZX[], vec *pSZy[], uword l, uword q, uword K);
void xtx_mt(mat& xtx, mat *pX[]);

void VB_update(mat *pX[], vec *py_bar[], mat &xtx, double logodds, vec &sb2, vec &se2, double &alpha,
               mat &mu, mat &s2, mat &alpha_jk, vec &pi_k, vec *pyt[], vec nn, double K);

void param_update(vec *py_bar[], mat &xtx, vec &sb2, vec &se2, double &alpha, mat &mu, mat &s2,
                  mat &alpha_jk, vec &pi_k, vec *pyt_j[], vec nn, double K);

void cov_update(vec *py[], mat *pZ[], mat &cov, vec *py_bar[], vec *pSZy[], mat *pSZX[], vec &pi_k, mat &alpha_jk, mat &mu);


double lb_linear(vec nn, vec &se2);
double lb_linear(vec *py_bar[], vec *pyt_j[], mat& xtx, mat& alpha_jk, mat& mu, vec& pi_k,
                 vec &se2, mat &s2, vec nn);
double lb_gamma(double alpha, mat& alpha_jk);
double lb_klbeta(vec &sb2, mat &s2, vec &pi_k, mat &alpha_jk, mat &mu);
#endif /* bivas_mt_aux_hpp */
