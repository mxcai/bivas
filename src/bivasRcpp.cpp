//
//  bivasRcpp.cpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/3/1.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include <stdio.h>
#include <iomanip>
#include <RcppArmadillo.h>
// #include <omp.h>
#include <R.h>
// #include <Rinternals.h>
#include "bivas_aux.hpp"//;
#include "bivas_mt_aux.hpp"//;
#include "bivas.hpp"//;


// [[ Rcpp :: depends ( RcppArmadillo )]]
using namespace arma;
using namespace Rcpp;
using namespace std;


// // [[Rcpp::export]]
// RcppExport SEXP outerloop_omp(arma::mat& X, arma::vec& y, arma::mat& Z, const arma::vec& group,
//                           const arma::vec& glevel, const arma::vec& logodds, arma::vec sb2,
//                           arma::vec se2, arma::vec alpha, arma::mat mu, arma::mat alpha_jk,
//                           arma::mat pi_k, arma::uword maxIter, double tol, int verbose, int coreNum) {
//
//     uword p = X.n_cols;
//     uword K = glevel.n_elem;
//     uword ns = pi_k.n_cols;
//
//     // centerize X and y
//     mat Xm = mean(X, 0);
//     double ym = mean(y);
//     X.each_row() -= Xm;
//     y -= ym;
//
//     // get invese
//     mat SZX;
//     vec SZy;
//     getInv(X, y, Z, SZX, SZy, p);
//
//     // do some pre-calculations
//     vec xtx = sum(square(X), 0).t();
//     mat pi_p(p, ns);
//     extend_pi_p(pi_k, pi_p, group, glevel, K);
//
//     mat cov(Z.n_cols, ns); // initialize vector for fixed effect estimates
//     vec logw(ns); //initialize vector for logweights
//     Rcpp::List Lq_List(ns); //initialize a list for Lq records of each prior setting
//
//
//     // outerloop for all prior settings
//
//     //set parallel computation
// #pragma omp parallel for num_threads(coreNum)
//     // #pragma omp parallel for if (coreNum > 1) num_threads(coreNum)
//     for(int i = 0; i < ns; i++) {
//
//       if(verbose == 1) {
//         printf("Start %d-th outer loop, logodds = %f \n", i+1, logodds(i));
//       }
//         outerloop(X, y, Z, SZX, SZy, xtx, Xm, ym, group, glevel, logodds, sb2, se2, alpha, cov, mu, alpha_jk, pi_k, pi_p,
//                   logw, maxIter, tol, verbose, i, Lq_List); // one iteration of outer loop
//
//     }
//
//     Rcpp::List ret;
//     ret["sb2"] = sb2;
//     ret["se2"] = se2;
//     ret["alpha"] = alpha;
//     ret["alpha_jk"] = alpha_jk;
//     ret["mu"] = mu;
//     ret["pi_k"] = pi_k;
//     ret["pi_p"] = pi_p;
//     ret["cov"] = cov;
//     ret["logw"] = logw;
//     ret["Lq_List"] = Lq_List;
//
//     return ret;
//
// }

// [[Rcpp::export]]
RcppExport SEXP outerloop_thread(arma::mat& X, arma::vec& y, arma::mat& Z, const arma::vec& group,
                          const arma::vec& glevel, const arma::vec& logodds, arma::vec sb2,
                          arma::vec se2, arma::vec alpha, arma::mat mu, arma::mat alpha_jk,
                          arma::mat pi_k, arma::uword maxIter, double tol, int verbose, const int coreNum) {




  uword p = X.n_cols;
  uword K = glevel.n_elem;
  uword ns = pi_k.n_cols;

  // centerize X and y
  mat Xm = mean(X, 0);
  double ym = mean(y);
  X.each_row() -= Xm;
  y -= ym;

  // get invese
  mat SZX;
  vec SZy;
  getInv(X, y, Z, SZX, SZy, p);

  // do some pre-calculations
  vec xtx = sum(square(X), 0).t();
  mat pi_p(p, ns);
  extend_pi_p(pi_k, pi_p, group, glevel, K);

  mat cov(Z.n_cols, ns); // initialize vector for fixed effect estimates
  vec logw(ns); //initialize vector for logweights
  Rcpp::List Lq_List(ns); //initialize a list for Lq records of each prior setting


  Bivas BivasObj(X, y, Z, SZX, SZy, xtx, Xm, ym, group, glevel, logodds, sb2, se2, alpha, cov, mu,
                    alpha_jk, pi_k, pi_p, logw, maxIter, tol, verbose, Lq_List);

  // outerloop for all prior settings

  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&Bivas::update_by_thread, &BivasObj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  Rcpp::List ret;
  ret["sb2"] = BivasObj.sb2;
  ret["se2"] = BivasObj.se2;
  ret["alpha"] = BivasObj.alpha;
  ret["alpha_jk"] = BivasObj.alpha_jkMat;
  ret["mu"] = BivasObj.muMat;
  ret["pi_k"] = BivasObj.pi_kMat;
  ret["pi_p"] = BivasObj.pi_pMat;
  ret["cov"] = BivasObj.covMat;
  ret["logw"] = BivasObj.logw;
  ret["Lq_List"] = BivasObj.Lq_List;

  return ret;

}




// [[Rcpp::export]]
RcppExport SEXP outerloop_mt_thread(Rcpp::List& Xs, Rcpp::List& ys, Rcpp::List& Zs, const arma::vec& logodds, arma::mat sb2,
                                 arma::mat se2, arma::vec alpha, arma::cube mu, arma::cube alpha_jk, arma::mat pi_k,
                                 arma::vec nn, arma::uword K, arma::uword q, arma::uword maxIter, double tol, int verbose, const int coreNum) {

  uword l = Xs.size();
  uword ns = pi_k.n_cols;

  mat *pX[l];
  vec *py[l];
  mat *pZ[l];

  unlist(Xs, pX);
  unlist(ys, py);
  unlist(Zs, pZ);

  // centerize X and y
  mat Xm(K, l);
  vec ym(l);
  centerize_mt(pX, py, Xm, ym);

  // get invese
  mat *pSZX[l];
  vec *pSZy[l];
  getInv_mt(pX, py, pZ, pSZX, pSZy, l, q, K);

  // do some pre-calculations
  mat xtx(K, l);
  xtx_mt(xtx, pX);

  cube cov(q+1, l, ns); // initialize ns matrices for fixed effect estimates, each is q*l
  vec logw(ns); //initialize vector for logweights
  Rcpp::List Lq_List(ns); //initialize a list for Lq records of each prior setting


  Bivas_mt Bivas_mtObj(pX, py, pZ, pSZX, pSZy, xtx, Xm, ym, logodds, sb2, se2, alpha, cov, mu,
                 alpha_jk, pi_k, l, nn, K, q, logw, maxIter, tol, verbose, Lq_List);

  // outerloop for all prior settings

  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&Bivas_mt::update_by_thread, &Bivas_mtObj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  // set to nullptr
  for(int j = 0; j < l; j++){
    delete pX[j];
    delete py[j];
    delete pZ[j];
    delete pSZX[j];
    delete pSZy[j];
  }

  Rcpp::List ret;
  ret["sb2"] = Bivas_mtObj.sb2Mat;
  ret["se2"] = Bivas_mtObj.se2Mat;
  ret["alpha"] = Bivas_mtObj.alpha;
  ret["alpha_jk"] = Bivas_mtObj.alpha_jkMat;
  ret["mu"] = Bivas_mtObj.muMat;
  ret["pi_k"] = Bivas_mtObj.pi_kMat;
  ret["logw"] = Bivas_mtObj.logw;
  ret["Lq_List"] = Bivas_mtObj.Lq_List;
  ret["cov"] = Bivas_mtObj.covMat;

  return ret;

}





// [[Rcpp::export]]
RcppExport SEXP tt(arma::vec s1, arma::vec s2, arma::mat xtx) {

  cout << s1 << endl;
  cout << s2 << endl;
  cout << xtx << endl;
  cout << xtx.each_row() + (s1 / s2).t() << endl;
  mat ss(xtx.n_rows,xtx.n_cols);
  ss = s1.t() / (xtx.each_row() + (s1 / s2).t()).each_row();
  cout << ss << endl;
  // return 0;
  //
  // cout << XX.slice(1).memptr() << endl;
  // mat pXX = mat(XX.slice(1).memptr(),2,3,false);
  // cout << pXX << endl;
  //
  // int l = Xs.size();
  // mat *mpX[l];
  // mat *mpX1[l];
  // mat **pX;
  // mat **pX1;
  // pX = new mat* [l];
  // // double** pointer_array;
  // // pointer_array = new double* [l];
  //
  // for (int i = 0; i < l; i++){
  //   SEXP xi = Xs[i];
  //   Rcpp::NumericMatrix Xr(xi);
  //
  //   int n = Xr.nrow();
  //   int k = Xr.ncol();
  //   // pointer_array[i] = Xr.begin();
  //
  //   mat t1 = mat(2,3);
  //   t1.fill(i);
  //   mpX1[i] = new mat(t1.begin(),2,3,true);
  //
  //  mpX[i] = new mat(Xr.begin(),n,k,false);
  //  pX[i] = new mat(Xr.begin(),n,k,false);
  //  // cout << pX[i] << endl;
  //  // cout << *pX[i] <<endl;
  //  // cout << *mpX[i] << endl;
  // }
  //
  // cout << *mpX1[0] << endl;
  // cout << *mpX1[1] << endl;
  //
  // pX1 = mpX;
  // cout << pX << endl;    //address of pointer pX
  // cout << pX[0] << endl;   //address of pointer pX[0]
  // cout << *pX[0] << endl;   //content in address pX[0]
  // cout << *pX[1] << endl;
  // cout << mpX << endl;
  // cout << pX1 << endl;
  // cout << pX1[0] << endl;
  // cout << (*pX1[0]).memptr() << endl;
  // cout << (*pX1[0]).begin() << endl;
  // cout << *mpX[0] << endl;
  // cout << *mpX[1] << endl;
  // cout << *pX1[0] << endl;
  // cout << *pX1[1] << endl;
  // // cout << (*mpX[0]).n_cols << endl;
  // // cout << *mpX[0] << endl;
  //
  // // for(int i = 0; i < l ; i++){
  // //   cout <<"index:" << i << endl;
  // //   double* p = pointer_array[i];
  // //   for(int j = 0; j < 3; j++)
  // //     cout << p[j] << " ";
  // //   cout << endl;
  // // }
  //
  Rcpp::List ret;

  return ret;

}

