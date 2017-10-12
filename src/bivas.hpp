//
//  bivas.hpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/16.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#ifndef bivas_hpp
#define bivas_hpp

#include <stdio.h>
#include <math.h>
#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include "bivas_aux.hpp"
#include "bivas_mt_aux.hpp"
#include <thread>
#include <mutex>


using namespace std;
using namespace arma;

class Bivas{
public:
  int current_idx=0;

  // std::mutex mtx;
  mat X;
  vec y;
  mat Z;
  mat SZX;
  mat SZy;
  vec xtx;
  mat Xm;
  double ym;
  vec group;
  vec glevel;

  vec logodds;

  vec sb2;
  vec se2;
  vec alpha;
  mat covMat;
  mat muMat;
  mat alpha_jkMat;
  mat pi_kMat;
  mat pi_pMat;
  vec logw;

  uword maxIter;
  double tol;
  int verbose;

  vec **Lq_list;

  Bivas(const mat& X, const vec& y, const mat& Z, const mat& SZX, const vec& SZy, const vec& xtx,
        mat& Xm, double ym, const vec& group, const vec& glevel, vec logodds, vec& sb2, vec& se2,
        vec& alpha, mat& covMat, mat& muMat, mat& alpha_jkMat, mat& pi_kMat, mat& pi_pMat,
        vec& logw, uword maxIter, double tol, int verbose, vec *Lq_list[]){

    this -> X = X;
    this -> y = y;
    this -> Z = Z;
    this -> SZX = SZX;
    this -> SZy = SZy;
    this -> xtx = xtx;
    this -> Xm = Xm;
    this -> ym = ym;
    this -> group = group;
    this -> glevel = glevel;
    this -> logodds = logodds;
    this -> sb2 = sb2;
    this -> se2 = se2;
    this -> alpha = alpha;
    this -> covMat = covMat;
    this -> muMat = muMat;
    this -> alpha_jkMat = alpha_jkMat;
    this -> pi_kMat = pi_kMat;
    this -> pi_pMat = pi_pMat;
    this -> logw = logw;
    this -> maxIter = maxIter;
    this -> tol = tol;
    this -> verbose = verbose;
    this -> Lq_list = Lq_list;

  }

  void outerloop_by_thread(int i);
  void update_by_thread(int thread_id);
  int  next();
};



class Bivas_mt{
public:
  int current_idx=0;

  std::mutex mtx;
  mat **pX;
  vec **py;
  mat **pZ;
  mat **pSZX;
  vec **pSZy;
  mat xtx;
  mat Xm;
  vec ym;

  vec logodds;

  mat sb2Mat;
  mat se2Mat;
  vec alpha;
  cube covMat;
  cube muMat;
  cube alpha_jkMat;
  mat pi_kMat;
  vec logw;

  vec nn;
  uword l;
  uword K;
  uword q;

  uword maxIter;
  double tol;
  int verbose;

  vec **Lq_list;

  Bivas_mt(mat *pX[], vec *py[], mat *pZ[], mat *pSZX[], vec *pSZy[], const mat& xtx,
        mat& Xm, vec& ym, vec logodds, mat& sb2Mat, mat& se2Mat, vec& alpha, cube &covMat, cube& muMat,
        cube& alpha_jkMat, mat& pi_kMat, uword l, vec nn, uword K, uword q, vec& logw, uword maxIter,
        double tol, int verbose, vec *Lq_list[]){

    this -> pX = pX;
    this -> py = py;
    this -> pZ = pZ;
    this -> pSZX = pSZX;
    this -> pSZy = pSZy;
    this -> xtx = xtx;
    this -> Xm = Xm;
    this -> ym = ym;
    this -> logodds = logodds;
    this -> sb2Mat = sb2Mat;
    this -> se2Mat = se2Mat;
    this -> alpha = alpha;
    this -> covMat = covMat;
    this -> muMat = muMat;
    this -> alpha_jkMat = alpha_jkMat;
    this -> pi_kMat = pi_kMat;
    this -> logw = logw;
    this -> maxIter = maxIter;
    this -> tol = tol;
    this -> verbose = verbose;
    this -> Lq_list = Lq_list;

    this -> nn = nn;
    this -> l = l;
    this -> K = K;
    this -> q = q;

  }

  void outerloop_by_thread(int i);
  void update_by_thread(int thread_id);
  int  next();
};


void outerloop(const mat& X, const vec& y, const mat& Z, const mat& SZX, const vec& SZy, const vec& xtx,
               mat& Xm, double ym, const vec& group, const vec& glevel, vec logodds, vec& sb2, vec& se2,
               vec& alpha, mat& covMat, mat& muMat, mat& alpha_jkMat, mat& pi_kMat, mat& pi_pMat,
               vec& logw, uword maxIter, double tol, int verbose, int i, Rcpp::List& Lq_List);


#endif /* bivas_hpp */
