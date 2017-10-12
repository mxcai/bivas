//
//  bivas.cpp
//  bivas
//
//  Created by CAI Mingxuan on 2017/2/16.
//  Copyright © 2017年 CAI Mingxuan. All rights reserved.
//

#include "bivas.hpp"

using namespace std;
using namespace arma;

// std::mutex mtx;
//
// int next(int current_idx, int total_task_num){
//   std::lock_guard<std::mutex> lock(mtx);
//   if(current_idx >= total_task_num){
//     return -1;
//   }
//
//   current_idx++;
//   std::cout << current_idx << std::endl;
//   return current_idx;
// }
//
//
// void outerloop_by_thread(int thread_id, const mat& X, const vec& y, const mat& Z, const mat& SZX, const vec& SZy,
//                          const vec& xtx, const vec& group, const vec& glevel, vec logodds, vec& sb2, vec& se2,
//                          vec& alpha, mat& covMat, mat& muMat, mat& alpha_jkMat, mat& pi_kMat, mat& pi_pMat,
//                          vec& logw, uword maxIter, double tol, int verbose, int i, Rcpp::List& Lq_List){
//   int total_task_num = logodds.size();
//   while(true){
//     int idx = next(idx, total_task_num);
//     if(idx == -1){
//       break;
//     }
//
//     outerloop(idx);
//   }
// }
//
//
// void cycle_threads(std::thread* threads, int n_thread){
//   for(int i = 0; i < n_thread; i++){
//     threads[i] = std::thread(outerloop_by_thread,i);
//   }
//   for(int i = 0; i < n_thread; i++){
//     threads[i].join();
// }


void outerloop(const mat& X, const vec& y, const mat& Z, const mat& SZX, const vec& SZy, const vec& xtx,
               mat& Xm, double ym, const vec& group, const vec& glevel, vec logodds, vec& sb2, vec& se2,
               vec& alpha, mat& covMat, mat& muMat, mat& alpha_jkMat, mat& pi_kMat, mat& pi_pMat,
               vec& logw, uword maxIter, double tol, int verbose, int i, Rcpp::List& Lq_List){



    uword n = y.n_elem;
    uword p = X.n_cols;
    uword K = glevel.n_elem;
    uword q = Z.n_cols;

    vec cov = vec(covMat.colptr(i), q, false);
    vec mu = vec(muMat.colptr(i), p, false);
    vec alpha_jk = vec(alpha_jkMat.colptr(i), p, false);
    vec pi_k = vec(pi_kMat.colptr(i), K, false);
    vec pi_p = vec(pi_pMat.colptr(i), p, false);

    // track Lq
    vec Lq(maxIter+1,fill::zeros);
    Lq(0) = -datum::inf;

    vec s2(p);
    vec pi_a_s_mu(p);
    vec offDiag(K);
    vec y_bar(n);

    cov = SZy - SZX * (pi_p % alpha_jk % mu);

    vec yt = X * (pi_p % alpha_jk % mu);


    // inner loop; update parameters and posterior
    for(int iter = 0; iter < 1; iter++) {
        s2 = se2(i) / (xtx + se2(i) / sb2(i));
        y_bar = y - Z*cov;

        VB_update(X, y_bar, xtx, group, glevel, logodds(i), sb2(i), se2(i), alpha(i), mu, s2,alpha_jk, pi_k,pi_p,
                  yt, offDiag, n, p, K);

        pi_a_s_mu = pi_p % alpha_jk % (s2 + square(mu));
        param_update(y_bar, xtx, sb2(i), se2(i), alpha(i), mu, s2, alpha_jk, pi_p, yt, offDiag, pi_a_s_mu, n, p);

        cov = SZy - SZX * (pi_p % alpha_jk % mu);


        Lq(iter+1) = lb_linear(y_bar, yt, xtx, pi_p, alpha_jk, mu, pi_a_s_mu, se2(i), offDiag) +
          lb_gamma(alpha(i), alpha_jk) + lb_eta(logodds(i), pi_k) +
          lb_klbeta(sb2(i), se2(i), mu, alpha_jk, pi_k, pi_p, s2, pi_a_s_mu);
        logw(i) = Lq(iter+1);

        if(verbose == 1){ //& (iter+1)%10==0) {
            printf("%d-th iteration: diff(Lq)=%f, alpha = %f, sigma2e = %f, sigma2beta = %f \n", iter+1, Lq(iter+1)-Lq(iter), alpha(i), se2(i), sb2(i));
            if(Lq(iter+1) < Lq(iter)){
                printf("Lowerbound decreasing at iteration %d th iteration \n",iter+1);
            }
        }

        if(abs(Lq(iter+1) - Lq(iter)) < abs(tol * Lq(iter))){
            if(verbose == 1){
                printf("Converge at %d th iteration \n",iter+1);
            }
            Lq = Lq.subvec(0,iter+1);
            break;
        }
    }
    // return the intercept
    // cov(0) += ym - dot(Xm,pi_p % alpha_jk % mu);
    //
    // Lq_List[i] = Lq;


}








void Bivas::outerloop_by_thread(int i){

  uword n = y.n_elem;
  uword p = X.n_cols;
  uword K = glevel.n_elem;
  uword q = Z.n_cols;

  vec cov = vec(covMat.colptr(i), q, false);
  vec mu = vec(muMat.colptr(i), p, false);
  vec alpha_jk = vec(alpha_jkMat.colptr(i), p, false);
  vec pi_k = vec(pi_kMat.colptr(i), K, false);
  vec pi_p = vec(pi_pMat.colptr(i), p, false);

  // track Lq
  vec Lq(maxIter+1,fill::zeros);
  Lq(0) = -datum::inf;

  vec s2(p);
  vec pi_a_s_mu(p);
  vec offDiag(K);
  vec y_bar(n);

  cov = SZy - SZX * (pi_p % alpha_jk % mu);

  vec yt = X * (pi_p % alpha_jk % mu);


  // inner loop; update parameters and posterior
  for(int iter = 0; iter < maxIter; iter++) {
    s2 = se2(i) / (xtx + se2(i) / sb2(i));
    y_bar = y - Z*cov;

    VB_update(X, y_bar, xtx, group, glevel, logodds(i), sb2(i), se2(i), alpha(i), mu, s2,alpha_jk, pi_k,pi_p,
              yt, offDiag, n, p, K);

    pi_a_s_mu = pi_p % alpha_jk % (s2 + square(mu));
    param_update(y_bar, xtx, sb2(i), se2(i), alpha(i), mu, s2, alpha_jk, pi_p, yt, offDiag, pi_a_s_mu, n, p);

    cov = SZy - SZX * (pi_p % alpha_jk % mu);


    Lq(iter+1) = lb_linear(y_bar, yt, xtx, pi_p, alpha_jk, mu, pi_a_s_mu, se2(i), offDiag) +
      lb_gamma(alpha(i), alpha_jk) + lb_eta(logodds(i), pi_k) +
      lb_klbeta(sb2(i), se2(i), mu, alpha_jk, pi_k, pi_p, s2, pi_a_s_mu);
    logw(i) = Lq(iter+1);

    if(verbose == 1){ //& (iter+1)%10==0) {
      printf("%d-th iteration: diff(Lq)=%f, alpha = %f, sigma2e = %f, sigma2beta = %f \n", iter+1, Lq(iter+1)-Lq(iter), alpha(i), se2(i), sb2(i));
      if(Lq(iter+1) < Lq(iter)){
        printf("Lowerbound decreasing at iteration %d th iteration \n",iter+1);
      }
    }

    if(abs(Lq(iter+1) - Lq(iter)) < abs(tol * Lq(iter))){
      if(verbose == 1){
        printf("Converge at %d th iteration \n",iter+1);
      }
      Lq = Lq.subvec(0,iter+1);
      break;
    }
  }
  // return the intercept
  cov(0) += ym - dot(Xm,pi_p % alpha_jk % mu);

  Lq_list[i] = new vec(Lq.begin(),Lq.n_elem,true);
}

// int current_idx=0;
std::mutex mtx;
int Bivas::next(){
  std::lock_guard<std::mutex> lock(mtx);
  if(current_idx >= logodds.n_elem){
    return -1;
  }
  current_idx++;

  return current_idx-1;
}

void Bivas::update_by_thread(int thread_id){
  while(true){
    int idx = next();

    if(idx == -1){
      break;
    }
    outerloop_by_thread(idx);
  }
  // std::lock_guard<std::mutex> lock(mtx);
  // for (int idx = 0; idx < logodds.n_elem; idx++){
  //   if(verbose == 1) {
  //     printf("Start %d-th outer loop, logodds = %f \n", idx+1, logodds(idx));
  //   }
  //   outerloop_by_thread(idx);
  // }
}








void Bivas_mt::outerloop_by_thread(int i){

  printf("Start %d -th outer loop: \n",i+1);
  uword ns = logodds.size();

  mat cov = mat(covMat.slice(i).memptr(), q+1, l, false);
  mat mu = mat(muMat.slice(i).memptr(), K, l, false);
  mat alpha_jk = mat(alpha_jkMat.slice(i).memptr(), K, l, false);
  vec pi_k = vec(pi_kMat.colptr(i), K, false);
  vec se2 = vec(se2Mat.colptr(i), l, false);
  vec sb2 = vec(sb2Mat.colptr(i), l, false);

  // track Lq
  vec Lq(maxIter+1,fill::zeros);
  Lq(0) = -datum::inf;

  mat s2(K, l);
  vec *py_bar[l];
  vec *pyt_j[l];

  for (int j = 0; j < l; j++){

    vec yt = *pX[j] * (pi_k % alpha_jk.col(j) % mu.col(j));
    pyt_j[j] = new vec(yt.begin(),nn(j),true);           //need to be null ptr later

    cov.col(j) = *pSZy[j] - *pSZX[j] * (pi_k % alpha_jk.col(j) % mu.col(j));
    vec y_bar = *py[j] - (*pZ[j] * cov.col(j));
    py_bar[j] = new vec(y_bar.begin(),nn(j),true);          //need to be null ptr later
  }

  // inner loop; update parameters and posterior
  for(int iter = 0; iter < maxIter; iter++) {

    VB_update(pX, py_bar, xtx, logodds(i), sb2, se2, alpha(i), mu, s2, alpha_jk, pi_k, pyt_j, nn, K);

    // Lq(iter+1) = lb_linear(py_bar, pyt_j, xtx, alpha_jk, mu, pi_k, se2, s2, nn) +
    //   lb_gamma(alpha(i), alpha_jk) + l*lb_eta(logodds(i), pi_k) +
    //   lb_klbeta(sb2, s2, pi_k, alpha_jk, mu);

    param_update(py_bar, xtx, sb2, se2, alpha(i), mu, s2, alpha_jk, pi_k, pyt_j, nn, K);

    cov_update(py, pZ, cov, py_bar, pSZy, pSZX, pi_k, alpha_jk, mu);

    // Lq(iter+1) = lb_linear(nn, se2) +
    Lq(iter+1) = lb_linear(py_bar, pyt_j, xtx, alpha_jk, mu, pi_k, se2, s2, nn) +
      lb_gamma(alpha(i), alpha_jk) + lb_eta(logodds(i), pi_k) +
      lb_klbeta(sb2, s2, pi_k, alpha_jk, mu);

    logw(i) = Lq(iter+1);

    if(verbose == 1){ //& (iter+1)%10==0) {
      // cout << Lq(iter+1) << endl;
      printf("%d-th iteration: diff(Lq)=%f, alpha = %f, sigma2e = %f, sigma2beta = %f \n", iter+1, Lq(iter+1)-Lq(iter), alpha(i), se2(0), sb2(0));
      // printf("%d-th iteration: alpha = %f, sigma2e = %f, sigma2beta = %f \n", iter+1, alpha(i), se2(0), sb2(0));
      if(Lq(iter+1) < Lq(iter)){
        printf("Lowerbound decreasing at iteration %d th iteration \n",iter+1);
      }
    }

    if(abs(Lq(iter+1) - Lq(iter)) < abs(tol * Lq(iter))){
      if(verbose == 1){
        printf("Converge at %d -th iteration. \n",iter+1);
      }

      Lq = Lq.subvec(0,iter+1);
      break;
    }
  }
  // set to nullptr
  for (int j = 0; j < l; j++){
    delete pyt_j[j];
    delete py_bar[j];
  }

  // return the intercept
  for( int j = 0; j < l; j++){
    cov(0,j) += ym(j) - dot(Xm.col(j),pi_k % alpha_jk.col(j) % mu.col(j));
  }


  Lq_list[i] = new vec(Lq.begin(),Lq.n_elem,true);
  printf("%d -th outer loop finished \n",i+1);
}

int Bivas_mt::next(){
  std::lock_guard<std::mutex> lockGuard(mtx);
  if(current_idx >= logodds.n_elem){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void Bivas_mt::update_by_thread(int thread_id){
  while(true){
    int idx = next();
    // cout << idx << endl;
    if(idx == -1){
      break;
    }
    outerloop_by_thread(idx);
  }
  // std::lock_guard<std::mutex> lock(mtx);
  // for (int idx = 0; idx < logodds.n_elem; idx++){
  //   if(verbose == 1) {
  //     printf("Start %d-th outer loop, logodds = %f \n", idx+1, logodds(idx));
  //   }
  //   outerloop_by_thread(idx);
  // }
}
