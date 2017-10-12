//
//  bivas_aux.hpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/19.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#ifndef bivas_aux_hpp
#define bivas_aux_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
// #include <Rcpp.h>

using namespace arma;
using namespace std;

void VB_update(const mat& X, const vec& y, const vec& xtx, const vec& group, const vec& glevel,
               double logodds, const double& sb2, const double& se2, const double& alpha, vec& mu,
               vec& s2, vec& alpha_jk, vec& pi_k, vec& pi_p, vec& yt, vec& offDiag, uword n, uword p,
               uword K);

void param_update(const vec& y, const vec& xtx, double& sb2, double& se2, double& alpha, const vec& mu,
                  const vec& s2, const vec& alpha_jk, const vec& pi_p, const vec& yt, const vec& offDiag,
                  const vec& pi_a_s_mu, uword n, uword p);

double lb_linear(const vec& y_bar, vec& yt, const vec& xtx, vec& pi_p, vec& alpha_jk, vec& mu,
                 vec& pi_a_s_mu, double se2, vec& offDiag);

double lb_gamma(double alpha, vec& alpha_jk);

double lb_eta(double logodds, vec& pi_k);

double lb_klbeta(double sb2, double se2, vec& mu, vec& alpha_jk, vec& pi_k, vec& pi_p, vec& s2, vec& pi_a_s_mu);

vec logpexp(vec& x);

double logpexp(double x);

void getInv(const mat& X, const vec& y, const mat& Z, mat& SZX, vec& SZy, uword p);

void extend_pi_p(mat& pi_k, mat& pi_p, const vec& group, const vec& glevel, uword K);

void extend_pi_p(vec& pi_k, vec& pi_p, uvec idx, uword n_idx, int k);

void ptr2List(Rcpp::List& Lst, vec *Lst_ptr[]);

#endif /* bivas_aux_hpp */
