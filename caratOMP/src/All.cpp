#include "RcppArmadillo.h"
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace std;

/*
 * GROUP-1: general functions (10 - 536)
 */


/*
 * same: 
 *   Judge whether two vector are same or not.
 */

// [[Rcpp::export]]
double Rprod(arma::vec v){
  double pd = arma::prod(v); 
  return pd; 
}

bool same(arma::vec a, arma::vec b){
  bool re = TRUE;
  if(a.n_elem == b.n_elem){
    for(unsigned int i = 0; i < a.n_elem; i++){
      if(a[i]!=b[i]){
        re = FALSE;break;
      }
    }
  }
  else{re = FALSE;}
  return re;
}

/*
 * PStrR: 
 *   Remove redundant duplicate columns of "M".
 */

arma::mat PStrR(arma::mat M){
  unsigned n = M.n_cols;
  arma::vec ind(n, arma::fill::zeros);
  unsigned int i, j;
  for(i = 0; i < n; i++){
    for(j = i + 1; j < n; j++){
      if(same(M.col(j), M.col(i))){
        ind(i) = 1; break;
      }
    }
  }
  return M.cols(find(ind == 0));
}

arma::mat permutations(int n){
  arma::mat Z(n, arma::prod(arma::linspace<arma::vec>(1, n, n))); 
  Z(0, 0) = 1;
  for(int i = 1; i < n; i++){
    int nip = arma::prod(arma::linspace<arma::vec>(1, i, i)); 
    Z.cols(0, nip - 1).row(i).fill(i + 1);
    arma::mat xtemp = Z.cols(0, nip - 1).rows(0, i);
    arma::uvec ind(2 * i + 1); 
    ind.subvec(0, i) = arma::linspace<arma::uvec>(0, i, i + 1);
    ind.subvec(i + 1, i + i) = arma::linspace<arma::uvec>(0, i - 1, i); 
    for(int j = 0; j < i; j++){
      Z.cols((j + 1) * nip, (j + 2) * nip - 1).rows(0, i) = xtemp.rows(ind.subvec(j + 1, j + 1 + i));
    }
  }
  return Z;
}

/*
 * Bpert: 
 *   Generate a matrix satisfies: 
 *     (1) with "bsize" rows;
 *     (2) each element is a one of numbers from 1 to "tr_num"
 *     (3) each column contains all numbers from 1 to "tr_num" 
 *         and numbers of 1 to "tr_num" are equal;
 *     (4) no duplications among columns.
 */

// [[Rcpp::export]]
arma::mat Bpert(int bsize = 4, int tr_num = 2){
  arma::mat BPT; 
  arma::mat BP = permutations(bsize);
  for(int i = 1; i < (bsize / tr_num); i++){
    for(int j = 1; j < (tr_num + 1); j++){
      BP(arma::find(BP == (i * tr_num + j))).fill(j);
    }
  }
  BPT = PStrR(BP);
  //arma::umat B = arma::conv_to<arma::umat>::from(BPT);
  return BPT;
}

/*
 * Csample: 
 *   Sample "num" numbers from 1 to n with or without replace with probability "proba".
 */

Rcpp::NumericVector Csample(int n, int num,
                            bool replace,
                            arma::vec proba) {
  Rcpp::NumericVector a(n);
  for(int i =0;i<n;i++){
    a(i) = i + 1;
  }
  Rcpp::NumericVector ret = Rcpp::RcppArmadillo::sample(a, num, replace, proba);
  return ret;
}

arma::rowvec sample_vec(int n, int num, arma::vec proba){
  proba = proba / sum(proba);
  arma::rowvec y = arma::randu<arma::rowvec>(num);
  arma::vec pro = cumsum(proba);
  y.elem(find(y <= pro(0))).fill(1);
  for(int i = 2; i < n; i++){
    y.elem(find(y >= pro(i - 2) && y <= pro(i - 1))).fill(i);
  }
  y.elem(find(y >= pro(n - 2) && y < 1)).fill(n);
  return y;
}

int sample_one(double prob){
  arma::vec T_One = arma::randu<arma::vec>(1);
  return(sum(T_One > prob) + 1);
}
/*
 * PStr:
 *   Generate matrix filled with 1 and 2, and no duplications among columns.
 */

arma::mat PStr(int n){
  arma::mat PStr;
  int i;
  arma::vec one(n), result(n);
  arma::mat final(n, pow(2, n));
  arma::uvec id;
  one(n - 1) = 1; result(n - 1) = 1; final(n - 1, 1) = 1;
  
  final.col(0) = one;
  for(i = 1; i < (pow(2, n) - 1); i++){
    result = result + one;
    while(sum(result == 2) > 0){
      id = find(result == 2);
      result(id - 1) = result(id - 1) + 1;
      result(id).fill(0);
    }
    final.col(i) = result;
  }
  final(find(final != 1.0)).fill(2);
  return final;
}

/*
 * PStrGen: 
 *   (1) Inputs: length of level_num should equal to cov_num;
 *   (2) Generate a matrix satisfis: 
 *       (a) with cov_num rows; 
 *       (b) elements in row i(i = 1,..., cov_num) are 
 *           a number among 1 to number which equals to i-th element of level_num.
 */

// [[Rcpp::export]]
arma::mat PStrGen(unsigned int cov_num, arma::vec level_num){
  arma::mat PS, P;
  arma::uvec ind = find(level_num >= 2.000001);
  arma::vec lel = level_num(ind);
  int leng = lel.n_elem;
  arma::mat PStr_temp = PStr(cov_num);
  int dimen = PStr_temp.n_cols;
  PS = arma::zeros<arma::mat>(cov_num, arma::prod(lel - 1) * dimen);
  if(leng == 0){
    return PStr_temp;
  }
  else{
    PS.cols(0, dimen - 1) = PStr_temp;
    int dimn = dimen;
    for(int j = 0; j < leng; j++){
      for(int l = 0; l < (lel(j) - 2); l++){
        PS.cols((l + 1) * dimn, (l + 2) * dimn - 1) = PS.cols(0, dimn - 1);
        arma::uvec b1(1); b1(0) = ind(j);
        arma::rowvec b2(dimn); b2.fill(l + 3);
        PS.submat(b1, arma::linspace<arma::uvec>((l + 1) * dimn, ((l + 2) * dimn - 1), dimn)) = b2;
      }
      dimn *= (lel(j) - 1);
    }
    P = PStrR(PS); 
    return P;
  }
}


arma::vec Prob(unsigned int cov_num, arma::vec level_num, arma::mat PStr, arma::vec pr){
  int lev_max = max(level_num);
  arma::mat p_mat(lev_max, cov_num, arma::fill::zeros);
  int sums = 0;
  for(unsigned int i = 0; i < cov_num; i++){
    arma::vec p_temp = pr.subvec(sums, (sums + level_num(i) - 1));
    p_mat.col(i).rows(0, (level_num(i) - 1)) = p_temp;
    sums += level_num(i);
  }
  arma::rowvec valv = sum(p_mat, 0);
  arma::vec temp = valv(find(valv <= 0.999999 || valv >= 1.000001));
  if(temp.n_elem > 0){
    return pr;
  }
  if(pr.n_elem != sum(level_num)){
    return pr;
  }
  else{
    int dimen = PStr.n_cols;
    arma::vec p(dimen);
    for(int i = 0; i < dimen; i++){
      arma::vec bridg = PStr.col(i);
      double prob = 1.0;
      for(unsigned int j = 0; j < cov_num; j++){
        prob *= p_mat((bridg(j) - 1), j);
      }
      p(i) = prob;
    }
    return p;
  }
}

// [[Rcpp::export]]
arma::mat Prob_S(unsigned int cov_num, arma::vec level_num, arma::vec pr){
  int lev_max = max(level_num);
  arma::mat p_mat(lev_max, cov_num, arma::fill::zeros), p_matN;
  int sums = 0;
  for(unsigned int i = 0; i < cov_num; i++){
    arma::vec p_temp = pr.subvec(sums, (sums + level_num(i) - 1));
    p_mat.col(i).rows(0, (level_num(i) - 1)) = p_temp;
    sums += level_num(i);
  }
  return p_mat;
}

arma::uvec ReturnCol(arma::mat M, arma::vec V){
  unsigned int leng = V.n_elem;
  arma::uvec Ind = find(M.row(0) == V(0));
  for(unsigned int i = 1; i < leng; i++){
    arma::uvec a(1); a(0) = i;
    arma::uvec ind = find(M.submat(a, Ind) == V(i));
    Ind = Ind(ind);
  }
  return Ind + 1;
}

// [[Rcpp::export]]
Rcpp::StringVector nameString(unsigned int cov_num, arma::vec level_num, int strt_num,
                              Rcpp::String type, Rcpp::String typeData){
  
  Rcpp::Environment base("package:base");
  Rcpp::Function paste = base["paste"];
  
  int suml = sum(level_num);
  
  if(type == "All" && typeData == "Simulation"){
    int prods = arma::prod(level_num);
    
    Rcpp::StringVector name(1 + prods + suml);
    name(0) = "overall";
    
    for(int i = 0; i < prods; i++){
      SEXP namL_temp = paste("level", i + 1, Rcpp::Named("sep", ""));
      Rcpp::String namL = Rcpp::as<Rcpp::String>(namL_temp);
      name(i + 1) = namL;
    }
    int sums = 0;
    for(unsigned int l = 0; l < cov_num; l++){
      for(int j = 0; j < level_num(l); j++){
        SEXP namM_temp = paste("margin", l + 1, j + 1, Rcpp::Named("sep", ""));
        Rcpp::String namM = Rcpp::as<Rcpp::String>(namM_temp);
        name(1 + prods + sums + j) = namM;
      }
      sums += level_num(l);
    }
    return name;
  }
  if(typeData == "Simulation" && type == "Margin"){
    Rcpp::StringVector name(cov_num);
    for(unsigned int i = 0; i < cov_num; i++){
      SEXP nameM_temp = paste("margin", i + 1, Rcpp::Named("sep", ""));
      Rcpp::String nameM = Rcpp::as<Rcpp::String>(nameM_temp);
      name(i) = nameM;
    }
    return name;
  }
  if(typeData == "Real" && type == "All"){
    Rcpp::StringVector name(1 + strt_num + suml);
    name(0) = "overall";
    for(int i = 0; i < strt_num; i++){
      SEXP nameRL_temp = paste("level", i + 1, Rcpp::Named("sep", ""));
      Rcpp::String nameRL = Rcpp::as<Rcpp::String>(nameRL_temp);
      name(i + 1) = nameRL;
    }
    int sumsR = 0;
    for(unsigned int l = 0; l < cov_num; l++){
      for(int j = 0; j < level_num(l); j++){
        SEXP nameRM_temp = paste("margin", l + 1, j + 1, Rcpp::Named("sep", ""));
        Rcpp::String nameRM = Rcpp::as<Rcpp::String>(nameRM_temp);
        name(1 + strt_num + sumsR + j) = nameRM;
      }
      sumsR += level_num(l);
    }
    return name;
  }
  else{
    Rcpp::Rcout<<"Please put in correct parameters!"<<std::endl;
    Rcpp::StringVector name;
    return name;
  }
}

// [[Rcpp::export]]
Rcpp::CharacterVector BBCDname(int n, Rcpp::String S = "pat"){
  Rcpp::Environment base("package:base");
  Rcpp::Function paste = base["paste"];
  Rcpp::CharacterVector name(n);
  for(int i = 0; i < n; i++){
    SEXP a = paste(S, i + 1, Rcpp::Named("sep", ""));
    Rcpp::String b = Rcpp::as<Rcpp::String>(a);
    name(i) = b;
  }
  return name;
}


// [[Rcpp::export]]
Rcpp::List Preprocess(Rcpp::DataFrame data){
  Rcpp::DataFrame comp;
  if(typeid(data).name() != typeid(comp).name()){
    Rcpp::Rcout<<"Error in data: data type must be Rcpp::DataFrame!"<<std::endl;
    return Rcpp::List::create(Rcpp::Named("data") = data);
  }
  else{
    Rcpp::Environment base("package:base");
    Rcpp::Function numeric = base["as.numeric"];
    Rcpp::Function sapply = base["sapply"];
    
    SEXP data_temp = sapply(data, numeric);
    Rcpp::NumericMatrix dat = Rcpp::as<Rcpp::NumericMatrix>(data_temp);
    arma::mat data_new = Rcpp::as<arma::mat>(dat);
    
    unsigned int cov_num = data_new.n_cols;
    arma::vec level_num = max(trans(data_new), 1);
    
    return Rcpp::List::create(Rcpp::Named("data") = trans(data_new),
                              Rcpp::Named("cov_num") = cov_num,
                              Rcpp::Named("level_num") = level_num);
  }
}

List Preprocess_out(DataFrame data){
  DataFrame comp;
  if(typeid(data).name() != typeid(comp).name()){
    Rcpp::Rcout<<"Error in data: data type must be dataframe!"<<std::endl;
    return List::create(Named("data") = data);
  }
  else{
    Environment base("package:base");
    Function numeric = base["as.numeric"];
    Function sapply = base["sapply"];
    Function matrix = base["apply"];
    
    SEXP data_temp = sapply(data, numeric);
    NumericMatrix dat = as<NumericMatrix>(data_temp);
    arma::mat data_new = as<arma::mat>(dat);
    
    int cov_num = data_new.n_rows-2;
    arma::vec level_num = max(data_new.rows(0,cov_num-1), 1);
    
    return List::create(Named("data") = data_new,
                        Named("cov_num") = cov_num,
                        Named("level_num") = level_num);
  }
}

arma::mat remd(arma::mat A, int bsize){
  int nr = A.n_rows, nc = A.n_cols; 
  arma::mat B(nr, nc);
  for(int i = 0; i < nr; i++){
    for(int j = 0; j < nc; j++){
      B(i, j) = std::remainder(A(i, j), bsize);
    }
  }
  B(arma::find(B < 0)) = bsize + B(arma::find(B < 0));
  return B;
}

arma::field<arma::mat> Analyze(arma::mat DIF, int row, int strt_num, int bsize, 
                               arma::mat ASSIG){
  arma::field<arma::mat> Res(2, 1); 
  
  arma::mat AD = abs(DIF);
  
  Rcpp::Environment stats("package:stats");
  Rcpp::Function quantile = stats["quantile"];
  
  arma::vec prb(3); prb(0) = 1; prb(1) = 0.95; prb(2) = 0.5;
  Rcpp::NumericMatrix QUA(row + 1, 3); QUA.fill(arma::datum::nan); 
  for(int l = 0; l < (row + 1); l++){
    SEXP q = quantile(AD.row(l), prb);
    Rcpp::NumericVector qua = Rcpp::as<Rcpp::NumericVector>(q);
    QUA.row(l) = qua;
  }
  arma::mat Q = Rcpp::as<arma::mat>(QUA);
  arma::mat R(row + 1, 4); 
  R.cols(0, 2) = Q; 
  
  arma::vec m = mean(AD, 1); 
  R.col(3) = m; 
  
  /*arma::vec SM(bsize);
  arma::mat Dstr = AD.rows(1, strt_num);
  
  arma::mat ypar = remd(SNUM, bsize);
  arma::uvec isum(bsize); isum(0) = bsize; 
  isum.subvec(1, bsize - 1) = arma::linspace<arma::uvec>(1, bsize - 1, bsize - 1);
  for(int s = 0; s < bsize; s++){
    arma::uvec ind = arma::find(ypar == s);
    arma::vec Dstr_c = Dstr(ind);
    if(Dstr_c.n_elem == 0){
      SM(isum(s) - 1) = arma::datum::nan;
    }
    else{
      double M = mean(Dstr_c);
      SM(isum(s) - 1) = M; 
    }
  }
  */
  
  Res(0, 0) = ASSIG; 
  Res(1, 0) = R; 
  //Res(2, 0) = SM; 
  
  return Res; 
}
/*
 * Output of Analyze: 
 * Rcpp::List::create(Rcpp::Named("Assign") = ASSIG, 
 * Rcpp::Named("All Imbalances") = Mix,
 * Rcpp::Named("Overall Imbalances") = Ovll,
 * Rcpp::Named("Marginal Imbalances") = MAR,
 * Rcpp::Named("Within-stratum Imbalances") = STR);
 */


/*
 * Simulate data
 */

//[[Rcpp::export]]
arma::mat genData_sim(int n, unsigned int cov_num, arma::vec level_num, 
                      arma::mat ProbS){
  arma::mat Data_S(cov_num, n);
  for(unsigned int k = 0; k < cov_num; k++){
    Data_S.row(k) = sample_vec(level_num(k),n,ProbS.col(k).head(level_num(k)));
  }
  return Data_S;
}

bool beta_check(unsigned int cov_num,arma::vec beta){
  if(beta.n_elem == cov_num){
    return TRUE;
  }
  else{
    Rcpp::Rcout<<"The length of beta must correspond to cov_num!"<<std::endl;
    Rcpp::Rcout<<"The length of beta should be:"<<cov_num<<std::endl;
    return FALSE;
  }
}

bool sigma_check(double sigma){
  if(sigma > 0){
    return TRUE;
  }
  else{
    Rcpp::Rcout<<"The error variance must be greater than zero!"<<std::endl;
    return FALSE;
  }
}


/*
 * GROUP-2: base functions
 */
// [[Rcpp::export]]
arma::field<arma::mat> HPSOne(arma::mat D, arma::mat PStr, arma::vec cov_profile, 
                              unsigned int cov_num, arma::vec level_num,
                              arma::vec omega, arma::vec strp, 
                              double p = 0.85){
  arma::field<arma::mat> result(4, 1); 
  int strt_num = PStr.n_cols; 
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  arma::uvec rvec = ReturnCol(PStr, cov_profile);
  int r = rvec(0);
  strp(r - 1) = strp(r - 1) + 1; 
  arma::uvec sub(2 + cov_num); sub(0) = 0; sub(1) = r;
  int temp = 0;
  for(unsigned int j = 0; j < cov_num; j++){
    sub(j + 2) = strt_num + temp + cov_profile(j);
    temp += level_num(j);
  }
  arma::uvec a(1); a(0) = 0;
  arma::vec arg = trans(omega) * D.submat(sub, a);
  double argval = arg(0, 0);
  
  arma::vec T_One(1); 
  if(argval > -0.000001 && argval < 0.000001){
    T_One = arma::randu<arma::vec>(1); 
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 0.5));
    result(0, 0) = strp; 
    result(1, 0) = cov_profile; 
    result(2, 0) = sum(T_One > 0.5) + 1; 
    result(3, 0) = D; 
    return result;
  }
  if(argval >= 0.000001){
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 1 - p));
    result(0, 0) = strp; 
    result(1, 0) = cov_profile; 
    result(2, 0) = sum(T_One > 1 - p) + 1; 
    result(3, 0) = D; 
    return result;
  }
  else{
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > p));
    result(0, 0) = strp; 
    result(1, 0) = cov_profile; 
    result(2, 0) = sum(T_One > p) + 1; 
    result(3, 0) = D; 
    return result;
  }
}

// [[Rcpp::export]]
arma::field<arma::mat> StrROne(arma::mat D, arma::mat PS, arma::vec cov_profile,
                               unsigned int cov_num, arma::vec level_num, int bsize, 
                               arma::mat B, arma::mat BG, arma::vec strp){
  arma::field<arma::mat> Res(4, 1); 
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  int Psize = PS.n_cols, Bsize = B.n_cols; 
  
  int r = ReturnCol(PS, cov_profile)(0);
  arma::uvec sub(2 + cov_num); sub(0) = 0; sub(1) = r;
  int temp = 0; 
  for(unsigned int j = 0; j < cov_num; j++){
    sub(2 + j) = Psize + temp + cov_profile(j);
    temp += level_num(j);
  }
  strp(r - 1) = strp(r - 1) + 1;
  arma::uvec ustrp = arma::conv_to<arma::uvec>::from(strp); 
  int pos = ustrp(r - 1) % bsize; 
  int T_STR;
  if(pos == 0){
    T_STR = BG(bsize - 1, r - 1); 
    BG.col(r - 1) = B.col(arma::randi<arma::vec>(Bsize, arma::distr_param(1, Bsize))(0) - 1);
  }
  else{
    T_STR = BG(pos - 1, r - 1);
  }
  D.rows(sub) = D.rows(sub) + brid(T_STR - 1);
  
  Res(0, 0) = strp; 
  Res(1, 0) = BG; 
  Res(2, 0) = T_STR; 
  Res(3, 0) = D; 
  return Res;
}

/*
 * input of AtBCDOne:
 * (1)D: difference between treatment 1 and treatment 2 at all three levels.
 *       The initial value is a vector with length (1 + strt_num + sum(level_num))
 *       and filled with 0.
 * (2)PStr: All involved strata.
 * (3)cov_profile: the covariate-profile of the coming patient.
 * (4)cov_num: numbers of considered covariates
 * (5)level_num: a vector whose element is the number of levels for each covariate.
 * (6)strt_num: the number of involved strata.
 * (7)num: a vector whose element is the number of patients falling in some group.
 *         The initial value is a length 2 vector filled with 0. It also the output 
 *         of BBCDOne. The output can be used in the next patient's input.
 * (8)F: a (1 + cov_num) * n matrix, the first row is always filled with 1. The initial value(i.e.
 *       when assign the 1st patient) is a (1 + cov_num) * n matrix with 1st row filled with 1. It also 
 *       the output of AtBCDOne. The output can be used in the next patient's input.
 * (13)No: the numbers of patients who are assigned to some group. The initial value is 0.
 *         It also the output of BBCDOne. The output can be used in the next patient's
 *         input.
 */

// [[Rcpp::export]]
arma::field<arma::mat> AtBCDOne(arma::vec D, arma::mat PStr, 
                                arma::vec cov_profile, unsigned int cov_num, 
                                arma::vec level_num, arma::mat F, 
                                arma::vec b, arma::vec strp, int No = 0){
  arma::field<arma::mat> Res(6, 1); 
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  int strt_num = PStr.n_cols; 
  arma::vec T_One(1); 
  if(No == 0){
    F.col(0).rows(1, cov_num) = cov_profile;
    int rf = ReturnCol(PStr, cov_profile)(0);
    strp(rf - 1) = strp(rf - 1) + 1; 
    T_One = arma::randu<arma::vec>(1); 
    int T = sum(T_One > 0.5);
    D(0) = D(0) + brid(T); D(rf) = D(rf) + brid(T);
    int temp = 0;
    for(unsigned int j = 0; j< cov_num; j++){
      int sub = strt_num + temp + cov_profile(j);
      D(sub) = D(sub) + brid(T);
      temp += level_num(j);
    }
    b(0) = brid(T);
    
    Res(0, 0) = strp; 
    Res(1, 0) = 1; 
    Res(2, 0) = F; 
    Res(3, 0) = b; 
    Res(4, 0) = T + 1; 
    Res(5, 0) = D; 
    return Res;
  }
  else{
    arma::vec bnt = F.cols(0, No - 1) * b.subvec(0, No - 1); 
    F.col(No).rows(1, cov_num) = cov_profile;
    
    arma::vec xn1 = F.col(No);
    arma::mat Fn = F.cols(0, No - 1);
    arma::mat man = Fn * trans(Fn);
    arma::mat mand = arma::pinv(man);
    arma::mat part = trans(xn1) * mand * bnt;
    double par = part(0, 0);
    double par1 = pow(1 - par, 2), par2 = pow(1 + par, 2);
    double p1 = par1 / (par1 + par2);
    
    T_One = arma::randu<arma::vec>(1); 
    int T = sum(T_One > p1);
    
    int r = ReturnCol(PStr, cov_profile)(0);
    strp(r - 1) = strp(r - 1) + 1; 
    D(0) = D(0) + brid(T); D(r) = D(r) + brid(T);
    int tp = 0;
    for(unsigned int jj = 0; jj < cov_num; jj++){
      int sub = strt_num + tp + cov_profile(jj);
      D(sub) = D(sub) + brid(T);
      tp += level_num(jj);
    }
    b(No) = brid(T);
    Res(0, 0) = strp; 
    Res(1, 0) = 1 + No; 
    Res(2, 0) = F; 
    Res(3, 0) = b; 
    Res(4, 0) = T + 1; 
    Res(5, 0) = D; 
    return Res;
  }
}

/*
 * input of AdBCDOne:
 * (1)D: difference between treatment 1 and treatment 2 at all three levels.
 *       The initial value is a vector with length (1 + strt_num + sum(level_num))
 *       and filled with 0.
 * (2)PStr: All involved strata.
 * (3)cov_profile: the covariate-profile of the coming patient.
 * (4)cov_num: numbers of considered covariates
 * (5)level_num: a vector whose element is the number of levels for each covariate.
 * (6)strt_num: the number of involved strata.
 * (7)weight: a vector with length sum(level_num). It weights all the imbalances at the
 *            each level of each covariate.
 * (8)a parameter. 
 */

// [[Rcpp::export]]
arma::field<arma::mat> AdBCDOne(arma::mat D, arma::mat PStr, arma::vec cov_profile, 
                                unsigned int cov_num, arma::vec level_num, 
                                arma::vec strp, double a = 2){
  arma::field<arma::mat> Res(3, 1); 
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  double strt_num = PStr.n_cols; 
  
  int r = ReturnCol(PStr, cov_profile)(0);
  strp(r - 1) = strp(r - 1) + 1; 
  arma::vec DAbs = abs(D); 
  double Dnd = DAbs(r, 0);
  double Dn = D(r, 0);
  arma::vec av(1); av(0) = a;
  arma::vec avec = abs(av);
  double p1, aa = avec(0);
  arma::vec T_One(1); 
  double pw = pow(Dnd, aa); 
  if(Dn >= -0.001 && Dnd <= 0.001){
    p1 = 0.5;
  }else if(Dn >= 0.99){
    p1 = 1.0 / (pw + 1.0);
  }else{
    p1 = pw / (pw + 1.0);
  }
  arma::vec pro(2); pro(0) = p1; pro(1) = 1 - p1;
  T_One = arma::randu<arma::vec>(1); 
  int T = sum(T_One > p1);
  D(0, 0) = D(0, 0) + brid(T); D(r, 0) = D(r, 0) + brid(T);
  int tps = 0;
  for(unsigned int jj = 0; jj < cov_num; jj++){
    int sub = strt_num + tps + cov_profile(jj);
    D(sub, 0) = D(sub, 0) + brid(T);
    tps += level_num(jj);
  }
  Res(0, 0) = strp; 
  Res(1, 0) = T + 1; 
  Res(2, 0) = D; 
  return Res;
}


/*
 * input of BBCDOne:
 * (1)D: difference between treatment 1 and treatment 2 at all three levels.
 *       The initial value is a vector with length (1 + strt_num + sum(level_num))
 *       and filled with 0.
 * (2)PStr: All involved strata.
 * (3)cov_profile: the covariate-profile of the coming patient.
 * (4)cov_num: numbers of considered covariates
 * (5)level_num: a vector whose element is the number of levels for each covariate.
 * (6)strt_num: the number of involved strata.
 * (7)num: a vector whose element is the number of patients falling in some group.
 *         The initial value is a length 2 vector filled with 0. It also the output 
 *         of BBCDOne. The output can be used in the next patient's input.
 * (8)numJ: a matrix whose (i, j) element represents the number of patients falling
 *          within the i-th group with class j of the category. The initial value is 
 *          a 2 * J matrix filled with 0. And it only needed when mode = "ad_hoc_J".
 *          J is the number of classes of the considered single category. It also 
 *          the output of BBCDOne. The output can be used in the next patient's input.
 * (9)gamma: a trade of parameter. The default = 1.0.
 * (10)criterion: design criterion. Default = "D-op". The alternative is "A-op".
 * (11)mode: randomization procedure. "formal", "ad_hoc", "ad_hoc_J" are welcomed.
 * (12)J: see in (8). It should be larger than 2, and is only needed when mode = "ad_hoc_J".
 * (13)No: the numbers of patients who are assigned to some group. The initial value is 0.
 *         It also the output of BBCDOne. The output can be used in the next patient's
 *         input.
 */
/*// [[Rcpp::export]]
arma::field<arma::mat> BBCDOne(arma::vec D, arma::mat PStr, arma::vec cov_profile, 
                               unsigned int cov_num, arma::vec level_num, 
                               arma::mat numJ, arma::vec strp, int J = 2, int No = 0){
  arma::field<arma::mat> Res(5, 1); 
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  int strt_num = PStr.n_cols; 
  if(No == 0){
    arma::vec Ja(2); Ja.fill(1.0 / 2);
    int TJ = arma::randi<arma::vec>(J, arma::distr_param(1, J))(0); 
    int T = arma::randi<arma::vec>(2, arma::distr_param(1, 2))(0);
    int r1 = ReturnCol(PStr, cov_profile)(0);
    strp(r1 - 1) = strp(r1 - 1) + 1; 
    D(0) = D(0) + brid(T - 1); D(r1) = D(r1) + brid(T - 1);
    int tJ = 0;
    for(unsigned int v = 0; v < cov_num; v++){
      int sJ = strt_num + tJ + cov_profile(v); 
      D(sJ) = D(sJ) + brid(T - 1); 
      tJ += level_num(v);
    }
    numJ(T - 1, TJ - 1) = numJ(T - 1, TJ - 1) + 1;
    Res(0, 0) = strp; 
    Res(1, 0) = No + 1; 
    Res(2, 0) = numJ; 
    Res(3, 0) = T; 
    Res(4, 0) = D; 
    return Res;
  }
  else{
    int TJ = arma::randi<arma::vec>(J, arma::distr_param(1, J))(0);
    arma::vec Jb(2); 
    Jb(0) = 0; Jb(1) = 1 - numJ(0, TJ - 1) * J / No;
    double p1 = Jb.max();
    arma::vec Jpro(2); Jpro(0) = p1; Jpro(1) = 1 - p1;
    arma::vec T_One = arma::randu<arma::vec>(1); 
    int T = sum(T_One > p1);
    numJ(T, TJ - 1) = numJ(T, TJ - 1) + 1;
    
    int r = ReturnCol(PStr, cov_profile)(0);
    strp(r - 1) = strp(r - 1) + 1; 
    D(0) = D(0) + brid(T); D(r) = D(r) + brid(T);
    int tJ = 0;
    for(unsigned int v = 0; v < cov_num; v++){
      int sJ = strt_num + tJ + cov_profile(v); 
      D(sJ) = D(sJ) + brid(T); 
      tJ += level_num(v);
    }
    Res(0, 0) = strp; 
    Res(1, 0) = No + 1; 
    Res(2, 0) = numJ; 
    Res(3, 0) = T + 1; 
    Res(4, 0) = D; 
    return Res;
  }
}
*/
/*
 * GROUP-3: inner functions
 */
/*
 * Groups-1: functions for simulation under 7 randomizations.
 *  (1)(2)(3) C_OneTrial: Hu and Hu's general CAR; Shao' method; Pocock and Simon's procedure
 *  (4) C_StrR: stratified permuted block randomization; 
 *  (5) C_AtkinBCD: Atkinson's optimum biased coin design; 
 *  (6) C_AdjustBCD: covariate-adaptive biased cion design;
 *  (7) C_BayseBCD: biased coin design with a baysian bias.
 */

// [[Rcpp::export]]
arma::field<arma::mat> C_HPS(int n, unsigned int cov_num, arma::vec level_num, 
                             arma::mat ProbS, arma::mat PStr,
                             arma::vec omega, double p = 0.85){
  
  arma::field<arma::mat> Res(4, 1); 
  
  arma::mat Data_S = genData_sim(n, cov_num, level_num, ProbS); 
  //arma::mat PStr = PStrR(Data_S);
  int strt_num = PStr.n_cols;
  int lnum = sum(level_num); 
  arma::vec strp(strt_num, arma::fill::zeros); 
  
  arma::mat D(1 + strt_num + lnum, 1, arma::fill::zeros);
  arma::mat CovAssig(cov_num + 1, n, arma::fill::zeros);
  CovAssig.rows(0, cov_num - 1) = Data_S;
  for(int i = 0; i < n; i++){
    arma::field<arma::mat> result = HPSOne(D, PStr, Data_S.col(i), cov_num, level_num, 
                                           omega, strp, p);
    
    CovAssig(cov_num, i) = result(2, 0)(0, 0);
    arma::vec strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp; 
    arma::mat Diff = result(3, 0);
    D.submat(0, 0, strt_num + lnum, 0) = Diff;
  }
  Res(0, 0) = strp; 
  Res(1, 0) = PStr; 
  Res(2, 0) = CovAssig;
  Res(3, 0) = D; 
  return Res;
}


// [[Rcpp::export]]
arma::field<arma::mat> C_StrR(int n, unsigned int cov_num, arma::vec level_num,
                              arma::mat ProbS, arma::mat PS,
                              int bsize, int tr_num){
  arma::field<arma::mat> Res(4, 1); 
  arma::mat Data_S = genData_sim(n, cov_num, level_num, ProbS); 
  //arma::mat PS = PStrR(Data_S);
  
  arma::mat B = Bpert(bsize, tr_num);
  
  int Psize = PS.n_cols;
  int Bsize = B.n_cols;
  int lnum = sum(level_num); 
  
  arma::vec prb(Bsize); prb.fill(1.0 / Bsize);
  Rcpp::NumericVector id = Csample(Bsize, Psize, TRUE, prb) - 1;
  arma::uvec ind = Rcpp::as<arma::uvec>(id);
  arma::mat BG = B.submat(arma::linspace<arma::uvec>(0, bsize - 1, bsize), ind);
  
  arma::vec strp(Psize, arma::fill::zeros);
  
  arma::mat CovAssig(cov_num + 1, n, arma::fill::zeros);
  CovAssig.rows(0, cov_num - 1) = Data_S;
  
  arma::mat D(1 + Psize + sum(level_num), 1, arma::fill::zeros);
  
  arma::uvec a0(1); a0(0) = 0;
  for(int i = 0; i < n; i++){
    arma::vec cov_profile = Data_S.col(i);
    
    arma::field<arma::mat> result = StrROne(D, PS, cov_profile, cov_num, level_num, 
                                            bsize, B, BG, strp); 
    int T_STR = result(2, 0)(0, 0);
    CovAssig(cov_num, i) = T_STR;
    
    arma::mat D_temp = result(3, 0); 
    D.submat(arma::linspace<arma::uvec>(0, Psize + lnum, (1 + Psize + lnum)), a0) = D_temp; 
    arma::mat BG_temp = result(1, 0);
    BG.cols(0, Psize - 1) = BG_temp; 
    arma::vec strp_temp = result(0, 0); 
    strp.subvec(0, Psize - 1) = strp_temp; 
  }
  
  Res(0, 0) = strp; 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}

// [[Rcpp::export]]
arma::field<arma::mat> C_AtkinBCD(int n, unsigned int cov_num, 
                                  arma::vec level_num, 
                                  arma::mat ProbS, arma::mat PS){
  arma::field<arma::mat> Res(4, 1);
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  arma::mat Data = genData_sim(n, cov_num, level_num, ProbS); 
  //arma::mat PS = PStrR(Data);
  int strt_num = PS.n_cols;
  int lnum = sum(level_num); 
  
  arma::mat Dus(1 + strt_num + lnum, 1, arma::fill::zeros);
  arma::mat CovAssig(cov_num + 1, n, arma::fill::zeros);
  CovAssig.rows(0, cov_num - 1) = Data;
  arma::mat F(cov_num + 1, n); F.row(0).fill(1);
  arma::vec bn(n);
  arma::vec strp(strt_num, arma::fill::zeros); 
  
  for(int l = 0; l < n; l++){
    arma::vec cov_profile = Data.col(l); 
    arma::field<arma::mat> result = AtBCDOne(Dus, PS, cov_profile, cov_num, 
                                             level_num, F, bn, strp, l);
    F.col(l).rows(1, cov_num) = cov_profile; 
    arma::mat b = result(3, 0); 
    bn.subvec(0, l) = b.submat(0, 0, l, 0); 
    arma::mat D_temp = result(5, 0); 
    Dus.submat(0, 0, strt_num + lnum, 0) = D_temp; 
    CovAssig(cov_num, l) = result(4, 0)(0, 0);
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp.col(0); 
  }
  
  Res(0, 0) = strp; 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = Dus; 
  return Res;
}

// [[Rcpp::export]]
arma::field<arma::mat> C_AdjustBCD(int n, unsigned int cov_num, 
                                   arma::vec level_num, 
                                   arma::mat ProbS, arma::mat PS, 
                                   double a = 2.0){
  arma::field<arma::mat> Res(4, 1); 
  
  arma::mat Data = genData_sim(n, cov_num, level_num, ProbS);
  //arma::mat PS = PStrR(Data);
  int strt_num = PS.n_cols;
  unsigned int lnums = sum(level_num);
  
  arma::mat D(1 + strt_num + lnums, 1, arma::fill::zeros);
  arma::mat CovAssig(1 + cov_num, n); 
  CovAssig.rows(0, cov_num - 1) = Data; 
  arma::vec strp(strt_num, arma::fill::zeros); 
  
  for(int i = 0; i < n; i++){
    arma::vec cov_profile = Data.col(i);
    arma::field<arma::mat> result = AdBCDOne(D, PS, cov_profile, cov_num, 
                                             level_num, strp, a); 
    CovAssig(cov_num, i) = result(1, 0)(0, 0);
    arma::mat D_temp = result(2, 0); 
    D.col(0) = D_temp; 
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp.col(0); 
  }
  Res(0, 0) = strp; 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}

/*
 * Group-3: functions for complete real data: 
 */
// [[Rcpp::export]]
arma::field<arma::mat> C_RHPS(arma::mat data_proc, unsigned int cov_num, 
                              arma::vec level_num, 
                              arma::vec omega, double p = 0.85){
  
  arma::field<arma::mat> Res(4, 1); 
  arma::mat P = PStrR(data_proc);
  int strt_num = P.n_cols;
  int lnum = sum(level_num); 
  
  arma::mat D(1 + strt_num + lnum, 1, arma::fill::zeros);
  
  int n = data_proc.n_cols;
  
  arma::mat CovAssig(1 + cov_num, n, arma::fill::zeros);
  CovAssig.rows(0, cov_num - 1) = data_proc;
  arma::vec strp(strt_num, arma::fill::zeros); 
  for(int i = 0; i < n; i++){
    arma::vec cov_profile = data_proc.col(i);
    arma::field<arma::mat> result = HPSOne(D, P, cov_profile, cov_num, 
                                           level_num, omega, strp, p);
    CovAssig(cov_num, i) = result(2, 0)(0, 0);
    arma::mat Diff = result(3, 0);
    D.submat(0, 0, strt_num + lnum, 0) = Diff;
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp.col(0); 
  }
  Res(0, 0) = trans(strp); 
  Res(1, 0) = P; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}

// [[Rcpp::export]]
arma::field<arma::mat> C_RStrR(arma::mat data_proc, unsigned int cov_num, 
                               arma::vec level_num, int bsize, int tr_num = 2){
  arma::field<arma::mat> Res(4, 1); 
  int n = data_proc.n_cols; 
  
  arma::mat PS = PStrR(data_proc);
  int Psize = PS.n_cols; 
  int lnum = sum(level_num);
  
  arma::mat B = Bpert(bsize, tr_num);
  int Bsize = B.n_cols;
  
  arma::vec prb(Bsize); prb.fill(1.0 / Bsize);
  Rcpp::NumericVector id = Csample(Bsize, Psize, TRUE, prb) - 1;
  arma::uvec ind = Rcpp::as<arma::uvec>(id);
  arma::mat BG = B.submat(arma::linspace<arma::uvec>(0, bsize - 1, bsize), ind);
  
  arma::vec strp(Psize, arma::fill::zeros);
  
  arma::mat CovAssig(cov_num + 1, n, arma::fill::zeros);
  CovAssig.rows(0, cov_num - 1) =data_proc;
  
  arma::mat D(1 + Psize + lnum, 1, arma::fill::zeros);
  for(int i = 0; i < n; i++){
    arma::vec cov_profile = data_proc.col(i);
    
    arma::field<arma::mat> result = StrROne(D, PS, cov_profile, cov_num, 
                                            level_num, bsize, B, BG, strp);
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, Psize - 1) = strp_temp.col(0); 
    arma::mat BG_temp = result(1, 0); 
    BG.cols(0, Psize - 1) = BG_temp; 
    arma::mat D_temp = result(3, 0); 
    D.submat(0, 0, Psize + lnum, 0) = D_temp;
    CovAssig(cov_num, i) = result(2, 0)(0, 0);
  }
  Res(0, 0) = trans(strp); 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}

// [[Rcpp::export]]
arma::field<arma::mat> C_RAtkinBCD(arma::mat data_proc,
                                   unsigned int cov_num, 
                                   arma::vec level_num){
  arma::field<arma::mat> Res(4, 1); 
  int n = data_proc.n_cols; 
  
  arma::mat PS = PStrR(data_proc);
  int strt_num = PS.n_cols;
  int lnum = sum(level_num); 
  
  arma::mat D(1 + strt_num + lnum, 1, arma::fill::zeros);
  arma::mat CovAssig(1 + cov_num, n); 
  CovAssig.rows(0, cov_num - 1) = data_proc; 
  
  arma::mat F(cov_num + 1, n); F.row(0).fill(1);
  arma::vec bn(n); 
  arma::vec strp(strt_num, arma::fill::zeros); 
  
  for(int l = 0; l < n; l++){
    arma::vec cov_profile = data_proc.col(l); 
    arma::field<arma::mat> result = AtBCDOne(D, PS, cov_profile, cov_num, 
                                             level_num, F, bn, strp, l);
    arma::mat bn_temp = result(3, 0); 
    bn.subvec(0, l) = bn_temp.submat(0, 0, l, 0); 
    F.col(l).rows(1, cov_num) = cov_profile; 
    arma::mat D_temp = result(5, 0); 
    D.submat(0, 0, strt_num + lnum, 0) = D_temp; 
    
    CovAssig(cov_num, l) = result(4, 0)(0, 0); 
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp.col(0); 
  }
  Res(0, 0) = trans(strp); 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}

// [[Rcpp::export]]
arma::field<arma::mat> C_RAdjustBCD(arma::mat data_proc, 
                                    unsigned int cov_num, 
                                    arma::vec level_num, 
                                    double a = 2.0){
  arma::field<arma::mat> Res(4, 1); 
  
  int n = data_proc.n_cols; 
  
  arma::mat PS = PStrR(data_proc);
  int strt_num = PS.n_cols;
  unsigned int lnums = sum(level_num);
  arma::mat D(1 + strt_num + lnums, 1, arma::fill::zeros);
  arma::mat CovAssig(1 + cov_num, n); 
  CovAssig.rows(0, cov_num - 1) = data_proc; 
  arma::vec strp(strt_num, arma::fill::zeros); 
  
  for(int i = 0; i < n; i++){
    arma::vec cov_profile = data_proc.col(i);
    arma::field<arma::mat> result = AdBCDOne(D, PS, cov_profile, cov_num, 
                                             level_num, strp, a); 
    arma::mat D_temp = result(2, 0); 
    D.submat(0, 0, strt_num + lnums, 0) = D_temp; 
    CovAssig(cov_num, i) = result(1, 0)(0, 0); 
    arma::mat strp_temp = result(0, 0); 
    strp.subvec(0, strt_num - 1) = strp_temp; 
  }
  Res(0, 0) = trans(strp); 
  Res(1, 0) = PS; 
  Res(2, 0) = CovAssig; 
  Res(3, 0) = D; 
  return Res;
}


/*
 * GROUP-4: inner summarization functions
 */
/*
 * Group-4: functions for assessment and comparison for randomizations: 
 *     If Replace = F, PS = PStrR(Data_S); 
 *     If Replace = T, PS = PStrGen(cov_num, level_num); 
 */

arma::field<arma::mat> SSum(arma::mat PS, arma::mat Data_S, int n, unsigned int cov_num,
                            arma::vec level_num, arma::mat ProbS, int bsize, int tr_num, arma::mat B, 
                            int Psize, bool Replace){
  arma::field<arma::mat> Res(3, 1); 
  
  int Bsize = B.n_cols;
  arma::vec prb(Bsize); prb.fill(1.0 / Bsize);
  Rcpp::NumericVector id = Csample(Bsize, Psize, TRUE, prb) - 1;
  arma::uvec ind = Rcpp::as<arma::uvec>(id);
  arma::mat BG = B.submat(arma::linspace<arma::uvec>(0, bsize - 1, bsize), ind);
  
  if(Replace == FALSE){
    arma::vec strp(Psize, arma::fill::zeros);
    arma::mat Assig(n, 1, arma::fill::zeros);
    arma::mat D(1 + Psize + sum(level_num), 1, arma::fill::zeros);
    for(int i = 0; i < n; i++){
      arma::vec cov_profile = Data_S.col(i);
      arma::field<arma::mat> res = StrROne(D, PS, cov_profile, cov_num, 
                                           level_num, bsize, B, BG, strp); 
      arma::mat Diff = res(3, 0); 
      D.col(0) = Diff; 
      arma::mat strp_temp = res(0, 0); 
      strp.subvec(0, Psize - 1) = strp_temp; 
      Assig(i, 0) = res(2, 0)(0, 0); 
    }
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res;
  }
  else{
    arma::field<arma::mat> result = C_StrR(n, cov_num, level_num, ProbS, 
                                           PS, bsize, 2); 
    arma::mat D = result(3, 0); 
    arma::mat Assig = result(2, 0).row(cov_num); 
    arma::vec strp = result(0, 0).col(0); 
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res;
  }
}

arma::field<arma::mat> OSum(arma::mat PStr, arma::mat Data_S, int n, unsigned int cov_num, 
                            arma::vec level_num, arma::mat ProbS, int strt_num, 
                            arma::vec omega, bool R, double p = 0.85){
  arma::field<arma::mat> Res(3, 1); 
  int lnum = sum(level_num); 
  if(R == FALSE){
    arma::mat D(1 + strt_num + lnum, 1, arma::fill::zeros); 
    arma::vec strp(strt_num, arma::fill::zeros); 
    arma::mat Assig(n, 1); 
    for(int i = 0; i < n; i++){
      arma::vec cov_profile = Data_S.col(i);
      arma::field<arma::mat> result = HPSOne(D, PStr, cov_profile, cov_num, 
                                             level_num, omega, strp, p); 
      arma::mat strp_temp = result(0, 0); 
      strp.subvec(0, strt_num - 1) = strp_temp; 
      arma::mat Diff = result(3, 0); 
      D.col(0) = Diff; 
      Assig(i, 0) = result(2, 0)(0, 0); 
    }
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res; 
  }
  else{
    arma::field<arma::mat> resu = C_HPS(n, cov_num, level_num, ProbS, 
                                        PStr, omega, p);
    arma::mat strp = resu(0, 0); 
    arma::mat Assig = resu(2, 0).row(cov_num); 
    arma::mat D = resu(3, 0); 
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res; 
  }
}

arma::field<arma::mat> AtBCD(arma::mat PS, arma::mat Data, int n, unsigned int cov_num, 
                             arma::vec level_num, arma::mat ProbS, 
                             int strt_num, bool Replace){
  arma::field<arma::mat> Res(3, 1); 
  
  int lnum = sum(level_num); 
  if(Replace == FALSE){
    arma::mat D(1 + strt_num + lnum, 1, arma::fill::zeros);
    arma::mat Assig(n, 1);
    arma::vec strp(strt_num, arma::fill::zeros);
    arma::mat F(1 + cov_num, n); F.row(0).fill(1); 
    arma::vec b(n); 
    for(int i = 0; i < n; i++){
      arma::vec cov_profile = Data.col(i); 
      arma::field<arma::mat> result = AtBCDOne(D, PS, cov_profile, cov_num, 
                                               level_num, F, b, strp, i);
      arma::mat strp_temp = result(0, 0); 
      strp.subvec(0, strt_num - 1) = strp_temp.col(0); 
      arma::mat b_temp = result(3, 0); 
      b.subvec(0, i) = b_temp.submat(0, 0, i, 0); 
      F.col(i).rows(1, cov_num) = cov_profile; 
      arma::mat Diff = result(5, 0); 
      D.col(0) = Diff; 
      Assig(i, 0) = result(4, 0)(0, 0); 
    }
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res;
  }
  else{
    arma::field<arma::mat> resu = C_AtkinBCD(n, cov_num, level_num,
                                             ProbS, PS); 
    arma::mat strp = resu(0, 0); 
    arma::mat Assig = resu(2, 0).row(cov_num); 
    arma::mat D = resu(3, 0); 
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res; 
  }
}

arma::field<arma::mat> AdBCD(arma::mat PS, arma::mat Data, int n, unsigned int cov_num, arma::vec level_num, 
                             arma::mat ProbS, double a, int strt_num, bool Replace){
  arma::field<arma::mat> Res(3, 1); 
  int lnums = sum(level_num);
  
  if(Replace == FALSE){
    arma::mat D(1 + strt_num + lnums, 1, arma::fill::zeros);
    arma::mat Assig(n, 1); 
    arma::vec strp(strt_num, arma::fill::zeros); 
    for(int i = 0; i < n; i++){
      arma::vec cov_profile = Data.col(i);
      arma::field<arma::mat> result = AdBCDOne(D, PS, cov_profile, cov_num, 
                                               level_num, strp, a); 
      arma::mat strp_temp = result(0, 0); 
      strp.subvec(0, strt_num - 1) = strp_temp.col(0);
      arma::mat Diff = result(2, 0); 
      D.col(0) = Diff; 
      Assig(i, 0) = result(1, 0)(0, 0); 
    }
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res;
  }
  else{
    arma::field<arma::mat> resu = C_AdjustBCD(n, cov_num, level_num, 
                                              ProbS, PS, a); 
    arma::mat strp = resu(0, 0); 
    arma::mat Assig = resu(2, 0).row(cov_num); 
    arma::mat D = resu(3, 0); 
    Res(0, 0) = strp; 
    Res(1, 0) = Assig; 
    Res(2, 0) = D; 
    return Res; 
  }
}

// [[Rcpp::export]]
arma::field<arma::mat> C_Summarize(bool Replace, unsigned int cov_num, 
                                   arma::vec level_num, arma::vec pr, 
                                   Rcpp::String method, 
                                   arma::vec omega, double p, 
                                   int bsize, double a, 
                                   int n, int N){
  arma::field<arma::mat> RR(5, 1); 
  
  arma::mat ProbS = Prob_S(cov_num, level_num, pr); 
  
  if(Replace == FALSE){
    
    arma::mat Data_S = genData_sim(n, cov_num, level_num, ProbS); 
    arma::mat PS = PStrR(Data_S);
    int strt_num = PS.n_cols;
    int lnum = sum(level_num);
    int row = strt_num + lnum;
    
    if(method == "StrPBR"){
      arma::mat B = Bpert(bsize, 2); 
      arma::mat DSA(row + 1, N), ASSIGS(n, N);
      arma::mat SNUM(strt_num, N); SNUM.fill(arma::datum::nan);
      for(int iter = 0; iter < N; iter++){
        arma::field<arma::mat> RES = SSum(PS, Data_S, n, cov_num, level_num, 
                                          ProbS, bsize, 2, B, strt_num, FALSE);
        arma::mat DS = RES(2, 0);
        DSA.col(iter) = DS;
        arma::mat AS = RES(1, 0);
        ASSIGS.col(iter) = AS;
        arma::mat SS = RES(0, 0).col(0);
        SNUM.col(iter) = SS;
      }
      arma::field<arma::mat> result = Analyze(DSA, row, strt_num, bsize, ASSIGS); 
      RR(0, 0) = result(0, 0); 
      RR(1, 0) = result(1, 0); 
      RR(2, 0) = SNUM; 
      RR(3, 0) = PS; 
      RR(4, 0) = DSA; 
      return RR;
    }
    else{
      
      arma::mat DIF(row + 1, N), ASSIG(n, N);
      arma::mat SNUM(strt_num, N);
      if(method == "HuHuCAR" || method == "PocSimMIN" || method == "StrBCD"){
        for(int i = 0; i < N; i++){
          arma::field<arma::mat> Res = OSum(PS, Data_S, n, cov_num, level_num, ProbS,
                                            strt_num, omega, FALSE, p);
          arma::mat D = Res(2, 0);
          DIF.col(i) = D;
          arma::mat A = Res(1, 0); 
          ASSIG.col(i) = A;
          arma::mat S = Res(0, 0);
          SNUM.col(i) = S;
        }
      }
      if(method == "DoptBCD"){
        for(int kk = 0; kk < N; kk++){
          arma::field<arma::mat> RES = AtBCD(PS, Data_S, n, cov_num, level_num, ProbS, 
                                             strt_num, FALSE); 
          arma::vec D = RES(2, 0);
          DIF.col(kk) = D;
          arma::mat A = RES(1, 0);
          ASSIG.col(kk) = A;
          arma::mat S = RES(0, 0);
          SNUM.col(kk) = S;
        }
      }
      if(method == "AdjBCD"){
        for(int ad = 0; ad < N; ad++){
          arma::field<arma::mat> RES = AdBCD(PS, Data_S, n, cov_num, level_num, ProbS, a,
                                             strt_num, FALSE);
          arma::vec D = RES(2, 0); 
          DIF.col(ad) = D;
          arma::mat A = RES(1, 0);
          ASSIG.col(ad) = A;
          arma::mat S = RES(0, 0);
          SNUM.col(ad) = S;
        }
      }
      arma::field<arma::mat> result = Analyze(DIF, row, strt_num, bsize, ASSIG);
      RR(0, 0) = result(0, 0); 
      RR(1, 0) = result(1, 0); 
      RR(2, 0) = SNUM; 
      RR(3, 0) = PS; 
      RR(4, 0) = DIF; 
      return RR;
    }
  }
  else{
    
    arma::mat PS = PStrGen(cov_num, level_num);
    int strt_num = PS.n_cols; 
    int lnum = sum(level_num);
    int row = strt_num + lnum;
    arma::mat Data_S(cov_num, n, arma::fill::zeros); 
    
    if(method == "StrPBR"){
      
      arma::mat BS = Bpert(bsize, 2); 
      arma::mat DSA(row + 1, N), ASSIGS(n, N); 
      arma::mat SNUM(strt_num, N); 
      for(int iter = 0; iter < N; iter++){
        
        //arma::mat Data_S(cov_num, n); 
        arma::field<arma::mat> RES = SSum(PS, Data_S, n, cov_num, level_num, ProbS, bsize, 2, BS, 
                                          strt_num, TRUE);
        arma::mat DS = RES(2, 0);
        DSA.col(iter) = DS;
        arma::mat AS = RES(1, 0).t();
        ASSIGS.col(iter) = AS;
        arma::mat SS = RES(0, 0);
        SNUM.col(iter) = SS;
      }
      arma::field<arma::mat> result = Analyze(DSA, row, strt_num, bsize, ASSIGS);
      RR(0, 0) = result(0, 0); 
      RR(1, 0) = result(1, 0); 
      RR(2, 0) = SNUM; 
      RR(3, 0) = PS; 
      RR(4, 0) = DSA; 
      return RR;
    }
    else{
      arma::mat DIF(row + 1, N), ASSIG(n, N);
      arma::mat SNUM(strt_num, N); SNUM.fill(arma::datum::nan);
      
      if(method == "HuHuCAR" || method == "PocSimMIN" || method == "StrBCD"){
        for(int i = 0; i < N; i++){
          
          arma::mat Data_S(cov_num, n); 
          arma::field<arma::mat> RES = OSum(PS, Data_S, n, cov_num, level_num, ProbS,
                                            strt_num, omega, TRUE, p);
          arma::mat D = RES(2, 0);
          DIF.col(i) = D;
          arma::mat A = RES(1, 0).t(); 
          ASSIG.col(i) = A;
          arma::mat S = RES(0, 0); 
          SNUM.col(i) = S;
        }
      }
      if(method == "DoptBCD"){
        for(int Ak = 0; Ak < N; Ak++){
          //arma::mat Data_S(cov_num, n);
          arma::field<arma::mat> RES = AtBCD(PS, Data_S, n, cov_num, level_num, ProbS, strt_num, TRUE);
          arma::vec D = RES(2, 0);
          DIF.col(Ak) = D;
          arma::vec A = RES(1, 0).t();
          ASSIG.col(Ak) = A;
          arma::mat S = RES(0, 0);
          SNUM.col(Ak) = S;
        }
      }
      if(method == "AdjBCD"){
        for(int Dk = 0; Dk < N; Dk++){
          arma::mat Data_S(cov_num, n); 
          arma::field<arma::mat> RES = AdBCD(PS, Data_S, n, cov_num, level_num, ProbS, a, strt_num, TRUE);
          arma::mat D = RES(2, 0);
          DIF.col(Dk) = D;
          arma::mat A = RES(1, 0).t(); 
          ASSIG.col(Dk) = A;
          arma::mat S = RES(0, 0); 
          SNUM.col(Dk) = S;
        }
      }
      arma::field<arma::mat> result = Analyze(DIF, row, strt_num, bsize, ASSIG);
      RR(0, 0) = result(0, 0); 
      RR(1, 0) = result(1, 0); 
      RR(2, 0) = SNUM; 
      RR(3, 0) = PS; 
      RR(4, 0) = DIF; 
      return RR;
    }
  }
}

// [[Rcpp::export]]
arma::field<arma::mat> C_RSummarize(arma::mat data_proc, unsigned int cov_num, arma::vec level_num, 
                                    Rcpp::String method, arma::vec omega, double p, int bsize, 
                                    double a, int N){
  arma::field<arma::mat> RR(5, 1); 
  arma::mat PS = PStrR(data_proc);
  int strt_num = PS.n_cols;
  int lnum = sum(level_num);
  int row = strt_num + lnum;
  int n = data_proc.n_cols; 
  arma::mat ProbS; 
  
  if(method == "StrPBR"){
    arma::mat B = Bpert(bsize, 2); 
    arma::mat DSA(row + 1, N), ASSIGS(n, N);
    arma::mat SNUM(strt_num, N); SNUM.fill(arma::datum::nan);
    for(int iter = 0; iter < N; iter++){
      arma::field<arma::mat> RES = SSum(PS, data_proc, n, cov_num, level_num, 
                                        ProbS, bsize, 2, B, strt_num, FALSE);
      arma::mat DS = RES(2, 0);
      DSA.col(iter) = DS;
      arma::mat AS = RES(1, 0);
      ASSIGS.col(iter) = AS;
      arma::mat SS = RES(0, 0).col(0);
      SNUM.col(iter) = SS;
    }
    arma::field<arma::mat> result = Analyze(DSA, row, strt_num, bsize, ASSIGS);
    RR(0, 0) = result(0, 0); 
    RR(1, 0) = result(1, 0); 
    RR(2, 0) = SNUM; 
    RR(3, 0) = PS; 
    RR(4, 0) = DSA; 
    return RR;
  }
  else{
    
    arma::mat DIF(row + 1, N), ASSIG(n, N);
    arma::mat SNUM(strt_num, N);
    if(method == "HuHuCAR" || method == "StrBCD" || method == "PocSimMIN"){
      
      for(int i = 0; i < N; i++){
        arma::field<arma::mat> Res = OSum(PS, data_proc, n, cov_num, level_num, ProbS,
                                          strt_num, omega, FALSE, p);
        arma::mat D = Res(2, 0);
        DIF.col(i) = D;
        arma::mat A = Res(1, 0); 
        ASSIG.col(i) = A;
        arma::mat S = Res(0, 0);
        SNUM.col(i) = S;
      }
    }
    if(method == "DoptBCD"){
      for(int kk = 0; kk < N; kk++){
        arma::field<arma::mat> RES = AtBCD(PS, data_proc, n, cov_num, level_num, ProbS, 
                                           strt_num, FALSE); 
        arma::vec D = RES(2, 0);
        DIF.col(kk) = D;
        arma::mat A = RES(1, 0);
        ASSIG.col(kk) = A;
        arma::mat S = RES(0, 0);
        SNUM.col(kk) = S;
      }
    }
    if(method == "AdjBCD"){
      for(int ad = 0; ad < N; ad++){
        arma::field<arma::mat> RES = AdBCD(PS, data_proc, n, cov_num, level_num, ProbS, a,
                                           strt_num, FALSE);
        arma::vec D = RES(2, 0); 
        DIF.col(ad) = D;
        arma::mat A = RES(1, 0);
        ASSIG.col(ad) = A;
        arma::mat S = RES(0, 0);
        SNUM.col(ad) = S;
      }
    }
    arma::field<arma::mat> result = Analyze(DIF, row, strt_num, bsize, ASSIG);
    RR(0, 0) = result(0, 0); 
    RR(1, 0) = result(1, 0); 
    RR(2, 0) = SNUM; 
    RR(3, 0) = PS; 
    RR(4, 0) = DIF; 
    return RR;
  }
}





/*
 * Functions about testing
 */


arma::mat HuHuCAR_In(arma::mat P, arma::mat D, arma::vec cov_profile, unsigned int cov_num,
                     arma::vec level_num, arma::vec omega, double p = 0.85){
  int strt_num = P.n_cols; 
  
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  arma::uvec rvec = ReturnCol(P, cov_profile);
  int r = rvec(0);
  arma::uvec sub(2 + cov_num); sub(0) = 0; sub(1) = r;
  int temp = 0;
  for(unsigned int j = 0; j < cov_num; j++){
    sub(j + 2) = strt_num + temp + cov_profile(j);
    temp += level_num(j);
  }
  arma::uvec a(1); a(0) = 0;
  arma::vec arg = trans(omega) * D.submat(sub, a);
  double argval = arg(0, 0);
  
  arma::vec T_One(1); 
  if(argval > -0.000001 && argval < 0.000001){
    T_One = arma::randu<arma::vec>(1); 
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 0.5));
    D(D.n_rows - 1, 0) = sum(T_One > 0.5) + 1; 
    return D;
  }
  if(argval >= 0.000001){
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 1 - p));
    D(D.n_rows - 1, 0) = sum(T_One > 1 - p) + 1; 
    return D;
  }
  else{
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > p));
    D(D.n_rows - 1, 0) = sum(T_One > p) + 1;
    return D;
  }
}

arma::rowvec Assign(arma::mat data,arma::mat D,arma::mat P,int n,int cov_num,int strt_num,arma::vec level_num,arma::vec omega, double p){
  arma::rowvec assignew(n);
  for(int i = 0; i< n; i++){
    D = HuHuCAR_In(P, D, data.col(i).head(cov_num), cov_num, level_num, omega, p);
    assignew(i) = D(D.n_rows-1,0);
  }
  return assignew;
}

arma::rowvec AssignB(arma::mat data,arma::mat D,arma::mat P,int n,int cov_num,int strt_num,arma::vec level_num,arma::vec omega,arma::uvec shuffle, double p){
  arma::rowvec assignew(n);
  arma::mat datanew = data.cols(shuffle);
  for(int i = 0; i< n; i++){
    D = HuHuCAR_In(P, D, datanew.col(i).head(cov_num), cov_num, level_num, omega, p);
    assignew(i) = D(D.n_rows-1,0);
  }
  return assignew;
}

//[[Rcpp::export]]
arma::mat HuHuCAR_getData(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          double mu1,double mu2,double sigma,arma::vec omega,
                          double p){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    arma::mat ProbS = Prob_S(cov_num,level_num,pr);
    Tdata.rows(0,cov_num-1) = genData_sim(n,cov_num,level_num,ProbS);
    arma::vec level_num_data = max(Tdata.rows(0,cov_num-1), 1);
    arma::mat P = PStrR(Tdata.rows(0,cov_num-1));
    int strt_num = P.n_cols;
    arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
    Tdata.row(cov_num) = Assign(Tdata,D,P,n,cov_num,strt_num,level_num_data,omega,p);
    arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
    if(type == "logit"){
      arma::vec m = exp(yita)/(1+exp(yita));
      arma::vec y = arma::randu(n);
      y.elem(find(y>=m)).fill(2);
      y.elem(find(y<m)).fill(1);
      y.replace(2,0);
      Tdata.insert_rows(cov_num+1,y.t());
    }
    if(type == "linear"){
      if(sigma_check(sigma) == TRUE){
        arma::vec epsilon(n,arma::fill::randn);
        epsilon = epsilon*sigma;
        arma::vec y = yita + epsilon;
        Tdata.insert_rows(cov_num+1,y.t());
      }
      else{
        Tdata.zeros();
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}


//[[Rcpp::export]]
Rcpp::List HuHuCAR_RT(DataFrame data,double Reps,arma::vec omega,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  for(int i=0; i < Reps; i++){
    D.zeros();
    assignew = Assign(data_proc,D,P,n,cov_num,strt_num,level_num,omega,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = Reps,
                              Named("pval") = 0,
                              Named("estimate") = diff_data);
  }
  else{
    double index = min(find(diff>=diff_data));
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = index,
                              Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                     Named("estimate") = diff_data);
  }
}

//[[Rcpp::export]]
arma::vec HuHuCAR_BT(DataFrame data,double B,arma::vec omega,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data_proc.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data_proc.row(cov_num+1);
  arma::vec D(2 + strt_num + sum(level_num));
  arma::rowvec assignew(n);
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    assignew = AssignB(data_proc,D,P,n,cov_num,strt_num,level_num,omega,shuffle,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec r(4);
  r(0) = s;
  r(1) = vb;
  r(2) = s/vb;
  r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  return r;
}

//[[Rcpp::export]]
double HuHuCAR_RT_In(arma::mat data,double Reps,arma::vec omega,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data.n_cols;
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  arma::vec D(2 + strt_num + sum(level_num));
  arma::rowvec assignew(n);
  for(int i=0; i < Reps; i++){
    D.zeros();
    assignew = Assign(data,D,P,n,cov_num,strt_num,level_num,omega,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = arma::find(diff>=diff_data);
  if(fi.is_empty()){
    return 0;
  }
  else{
    double index = min(find(diff>=diff_data));
    return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
  }
}

//[[Rcpp::export]]
double HuHuCAR_BT_In(arma::mat data,double B,arma::vec omega,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data.row(cov_num+1);
  arma::vec D(2 + strt_num + sum(level_num));
  arma::rowvec assignew(n,arma::fill::zeros);
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    assignew = AssignB(data,D,P,n,cov_num,strt_num,level_num,omega,shuffle,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
}

// [[Rcpp::export]]
arma::vec HuHuCAR_RT_power(int n,unsigned int cov_num,arma::vec level_num,arma::vec pr,std::string type,arma::vec beta,arma::vec mu1,arma::vec mu2,double sigma,double Iternum,double sl,arma::vec omega,double p,double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = HuHuCAR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,omega,p);
        double pval = HuHuCAR_RT_In(data,Reps,omega,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec HuHuCAR_BT_power(int n,unsigned int cov_num,arma::vec level_num,arma::vec pr,std::string type,arma::vec beta,arma::vec mu1,arma::vec mu2,double sigma,double Iternum,double sl,arma::vec omega,double p,double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = HuHuCAR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,omega,p);
        double pval = HuHuCAR_BT_In(data,B,omega,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier   
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

//[[Rcpp::export]]
arma::mat PocSimMIN_getData(int n,unsigned int cov_num,arma::vec level_num,
                            arma::vec pr,std::string type,arma::vec beta,
                            double mu1,double mu2,double sigma,
                            arma::vec weight,double p){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    if(weight.n_elem != cov_num && weight.n_elem != 0){
      Rcpp::Rcout<<
        "Error in weight: length of weight should correspond to number of covariates for PocSimMIN!"
        <<std::endl;
        Rcpp::Rcout<<"Required lenth of weight :"<<cov_num<<std::endl;
        Tdata.zeros();
    }
    else{
      arma::vec omeganew(2 + cov_num, arma::fill::zeros);
      if(weight.n_elem == 0){
        omeganew.subvec(2, 1 + cov_num).fill(1 / cov_num);
      }
      else{
        arma::vec w = abs(weight) / sum(abs(weight));
        omeganew.subvec(2, 1 + cov_num) = w;
      }
      arma::mat ProbS = Prob_S(cov_num,level_num,pr);
      Tdata.rows(0,cov_num-1) = genData_sim(n,cov_num,level_num,ProbS);
      arma::vec level_num_data = max(Tdata.rows(0,cov_num-1), 1);
      arma::mat P = PStrR(Tdata.rows(0,cov_num-1));
      int strt_num = P.n_cols;
      arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
      Tdata.row(cov_num) = Assign(Tdata,D,P,n,cov_num,strt_num,level_num_data,omeganew,p);
      arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
      if(type == "logit"){
        arma::vec m = exp(yita)/(1+exp(yita));
        arma::vec y = arma::randu(n);
        y.elem(find(y>=m)).fill(2);
        y.elem(find(y<m)).fill(1);
        y.replace(2,0);
        Tdata.insert_rows(cov_num+1,y.t());
      }
      if(type == "linear"){
        if(sigma_check(sigma) == TRUE){
          arma::vec epsilon(n,arma::fill::randn);
          epsilon = epsilon*sigma;
          arma::vec y = yita + epsilon;
          Tdata.insert_rows(cov_num+1,y.t());
        }
        else{
          Tdata.zeros();
        }
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}

//[[Rcpp::export]]
Rcpp::List PocSimMIN_RT(DataFrame data,double Reps,arma::vec weight,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  Rcpp::List result;
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  if(weight.n_elem != cov_num && weight.n_elem != 0){
    Rcpp::Rcout<<
      "Error in weight: length of weight should correspond to number of covariates for PocSimMIN!"
      <<std::endl;
      Rcpp::Rcout<<"Required lenth of weight :"<<cov_num<<std::endl;
      return Rcpp::List::create(Named("Randata") = 0,
                                Named("index") = 0,
                                Named("pval") = 0,
                                Named("estimate") = 0);
  }
  else{
    arma::vec omeganew(2 + cov_num, arma::fill::zeros);
    if(weight.n_elem == 0){
      omeganew.subvec(2, 1 + cov_num).fill(1 / cov_num);
    }
    else{
      arma::vec w = abs(weight) / sum(abs(weight));
      omeganew.subvec(2, 1 + cov_num) = w;
    }
    for(int i=0; i < Reps; i++){
      D.zeros();
      assignew = Assign(data_proc,D,P,n,cov_num,strt_num,level_num,omeganew,p);
      n1 = -sum(assignew-2);
      n0 = n - n1;
      diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
    }
    diff = sort(diff);
    arma::uvec fi = find(diff>=diff_data);
    if(fi.is_empty()){
      return Rcpp::List::create(Named("Randata") = diff,
                                Named("index") = Reps,
                                Named("pval") = 0,
                                Named("estimate") = diff_data);
    }
    else{
      double index = min(find(diff>=diff_data));
      return Rcpp::List::create(Named("Randata") = diff,
                                Named("index") = index,
                                Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                       Named("estimate") = diff_data);
    }
  }
}

//[[Rcpp::export]]
arma::vec PocSimMIN_BT(DataFrame data,double B,arma::vec weight,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  if(weight.n_elem != cov_num && weight.n_elem != 0){
    Rcpp::Rcout<<
      "Error in weight: length of weight should correspond to number of covariates for PocSimMIN!"
      <<std::endl;
      Rcpp::Rcout<<"Required lenth of weight :"<<cov_num<<std::endl;
      return level_num.zeros();
  }
  else{
    arma::vec omeganew(2 + cov_num, arma::fill::zeros);
    if(weight.n_elem == 0){
      omeganew.subvec(2, 1 + cov_num).fill(1 / cov_num);
    }
    else{
      arma::vec w = abs(weight) / sum(abs(weight));
      omeganew.subvec(2, 1 + cov_num) = w;
    }
    arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
    int strt_num = P.n_cols;
    arma::uvec a0(1,arma::fill::zeros);
    int n = data_proc.n_cols;
    arma::uvec shuffle(n,arma::fill::zeros);
    arma::vec theta(B,arma::fill::zeros);
    Rcpp::List result;
    double n0,n1;
    arma::rowvec ynew(n);
    arma::rowvec y = data_proc.row(cov_num+1);
    arma::vec D(1 + strt_num + sum(level_num));
    arma::rowvec assignew(n,arma::fill::zeros);
    for(int i=0; i < B; i++){
      shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
      ynew = trans(y(shuffle));
      D.zeros();
      assignew = AssignB(data_proc,D,P,n,cov_num,strt_num,level_num,omeganew,shuffle,p);
      n1 = -sum(assignew-2);
      n0 = n - n1;
      theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
    }
    n1 = -sum(data_proc.row(cov_num)-2);
    n0 = n - n1;
    double vb = stddev(theta);
    double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
    arma::vec r(4);
    r(0) = s;
    r(1) = vb;
    r(2) = s/vb;
    r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
    return r;
  }
}

//[[Rcpp::export]]
double PocSimMIN_RT_In(arma::mat data,double Reps,arma::vec weight,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  Rcpp::List result;
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  if(weight.n_elem != cov_num && weight.n_elem != 0){
    Rcpp::Rcout<<
      "Error in weight: length of weight should correspond to number of covariates for PocSimMIN!"
      <<std::endl;
      Rcpp::Rcout<<"Required lenth of weight :"<<cov_num<<std::endl;
      return 0;
  }
  else{
    arma::vec omeganew(2 + cov_num, arma::fill::zeros);
    if(weight.n_elem == 0){
      omeganew.subvec(2, 1 + cov_num).fill(1 / cov_num);
    }
    else{
      arma::vec w = abs(weight) / sum(abs(weight));
      omeganew.subvec(2, 1 + cov_num) = w;
    }
    for(int i=0; i < Reps; i++){
      D.zeros();
      assignew = Assign(data,D,P,n,cov_num,strt_num,level_num,omeganew,p);
      n1 = -sum(assignew-2);
      n0 = n - n1;
      diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
    }
    diff = sort(diff);
    arma::uvec fi = find(diff>=diff_data);
    if(fi.is_empty()){
      return 0;
    }
    else{
      double index = min(find(diff>=diff_data));
      return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
    }
  }
}


//[[Rcpp::export]]
double PocSimMIN_BT_In(arma::mat data,double B,arma::vec weight,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  if(weight.n_elem != cov_num && weight.n_elem != 0){
    Rcpp::Rcout<<
      "Error in weight: length of weight should correspond to number of covariates for PocSimMIN!"
      <<std::endl;
      Rcpp::Rcout<<"Required lenth of weight :"<<cov_num<<std::endl;
      return 0;
  }
  else{
    arma::vec omeganew(2 + cov_num, arma::fill::zeros);
    if(weight.n_elem == 0){
      omeganew.subvec(2, 1 + cov_num).fill(1 / cov_num);
    }
    else{
      arma::vec w = abs(weight) / sum(abs(weight));
      omeganew.subvec(2, 1 + cov_num) = w;
    }
    arma::mat P = PStrR(data.rows(0,cov_num-1));
    int strt_num = P.n_cols;
    arma::uvec a0(1,arma::fill::zeros);
    int n = data.n_cols;
    arma::uvec shuffle(n,arma::fill::zeros);
    arma::vec theta(B,arma::fill::zeros);
    Rcpp::List result;
    double n0,n1;
    arma::rowvec ynew(n);
    arma::rowvec y = data.row(cov_num+1);
    arma::vec D(1 + strt_num + sum(level_num));
    arma::rowvec assignew(n,arma::fill::zeros);
    for(int i=0; i < B; i++){
      shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
      ynew = trans(y(shuffle));
      D.zeros();
      assignew = AssignB(data,D,P,n,cov_num,strt_num,level_num,omeganew,shuffle,p);
      n1 = -sum(assignew-2);
      n0 = n - n1;
      theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
    }
    n1 = -sum(data.row(cov_num)-2);
    n0 = n - n1;
    double vb = stddev(theta);
    double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
    return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  }
}

// [[Rcpp::export]]
arma::vec PocSimMIN_RT_power(int n,unsigned int cov_num,arma::vec level_num,arma::vec pr,std::string type,arma::vec beta,arma::vec mu1,arma::vec mu2,double sigma,double Iternum,double sl,arma::vec weight,double p,double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = PocSimMIN_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,weight,p);
        double pval = PocSimMIN_RT_In(data,Reps,weight,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec PocSimMIN_BT_power(int n,unsigned int cov_num,arma::vec level_num,
                             arma::vec pr,std::string type,arma::vec beta,
                             arma::vec mu1,arma::vec mu2,double sigma,
                             double Iternum,double sl,arma::vec weight,
                             double p,double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = PocSimMIN_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,weight,p);
        double pval = PocSimMIN_BT_In(data,B,weight,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

//[[Rcpp::export]]
arma::mat StrBCD_getData(int n,unsigned int cov_num,arma::vec level_num,
                         arma::vec pr,std::string type,arma::vec beta,double mu1,
                         double mu2,double sigma,double p){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
    arma::mat ProbS = Prob_S(cov_num,level_num,pr);
    Tdata.rows(0,cov_num-1) = genData_sim(n,cov_num,level_num,ProbS);
    arma::vec level_num_data = max(Tdata.rows(0,cov_num-1), 1);
    arma::mat P = PStrR(Tdata.rows(0,cov_num-1));
    int strt_num = P.n_cols;
    arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
    Tdata.row(cov_num) = Assign(Tdata,D,P,n,cov_num,strt_num,level_num_data,omeganew,p);
    arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
    if(type == "logit"){
      arma::vec m = exp(yita)/(1+exp(yita));
      arma::vec y = arma::randu(n);
      y.elem(find(y>=m)).fill(2);
      y.elem(find(y<m)).fill(1);
      y.replace(2,0);
      Tdata.insert_rows(cov_num+1,y.t());
    }
    if(type == "linear"){
      if(sigma_check(sigma) == TRUE){
        arma::vec epsilon(n,arma::fill::randn);
        epsilon = epsilon*sigma;
        arma::vec y = yita + epsilon;
        Tdata.insert_rows(cov_num+1,y.t());
      }
      else{
        Tdata.zeros();
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}

//[[Rcpp::export]]
Rcpp::List StrBCD_RT(DataFrame data,double Reps,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
  for(int i=0; i < Reps; i++){
    D.zeros();
    assignew = Assign(data_proc,D,P,n,cov_num,strt_num,level_num,omeganew,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = Reps,
                              Named("pval") = 0,
                              Named("estimate") = diff_data);
  }
  else{
    double index = min(find(diff>=diff_data));
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = index,
                              Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                     Named("estimate") = diff_data);
  }
}

//[[Rcpp::export]]
arma::vec StrBCD_BT(DataFrame data,double B,double p){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data_proc.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data_proc.row(cov_num+1);
  arma::vec D(1 + strt_num + sum(level_num));
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    assignew = AssignB(data_proc,D,P,n,cov_num,strt_num,level_num,omeganew,shuffle,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec r(4);
  r(0) = s;
  r(1) = vb;
  r(2) = s/vb;
  r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  return r;
}

//[[Rcpp::export]]
double StrBCD_RT_In(arma::mat data,double Reps,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
  for(int i=0; i < Reps; i++){
    D.zeros();
    assignew = Assign(data,D,P,n,cov_num,strt_num,level_num,omeganew,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return 0;
  }
  else{
    double index = min(find(diff>=diff_data));
    return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
  }
}

//[[Rcpp::export]]
double StrBCD_BT_In(arma::mat data,double B,double p){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data.row(cov_num+1);
  arma::vec D(1 + strt_num + sum(level_num));
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    assignew = AssignB(data,D,P,n,cov_num,strt_num,level_num,omeganew,shuffle,p);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
}


// [[Rcpp::export]]
arma::vec StrBCD_RT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double p,double Iternum,double sl,
                          double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,p);
        double pval = StrBCD_RT_In(data,Reps,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec StrBCD_BT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,
                          double sigma,double Iternum,double sl,
                          double p,double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,p);
        double pval = StrBCD_BT_In(data,B,p);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::mat DoptBCD_In(int n, unsigned int cov_num,arma::vec level_num, arma::vec pr){
  
  arma::mat ProbS = Prob_S(cov_num, level_num, pr);
  arma::mat Data = genData_sim(n,cov_num,level_num,ProbS);
  arma::mat P = PStrR(Data);
  
  arma::rowvec assig(n);
  arma::mat Tdata(cov_num+1,n);
  arma::vec T_One(1);
  if(level_num.n_elem != cov_num){
    Rcpp::Rcout<<
      "Error in level_num: default of level_num is (2,  2). Length of level_num should match cov_num!"
      <<std::endl;
      return Tdata.zeros();
  }
  arma::vec val = level_num(find(level_num <= 1.999999));
  if(val.n_elem > 0){
    Rcpp::Rcout<<
      "Error in level_num: number of levels of each covariate must be larger than 2!"
      <<std::endl;
      return Tdata.zeros();
  }
  if(pr.n_elem != sum(level_num)){
    Rcpp::Rcout<<
      "Error in pr: default of pr is (0.5, 0.5, 0.5, 0.5). Length of pr should correspond to number of all levels for all margin!"
      <<std::endl;
      return Tdata.zeros();
  }
  arma::rowvec valv = sum(ProbS, 0);
  arma::vec valvCh = valv(find(valv <=0.999999 || valv >= 1.000001));
  if(valvCh.n_elem > 0){
    Rcpp::Rcout<<"Error in pr: probabilities of each margin must sum up to 1!"<<std::endl;
    return Tdata.zeros();
  }
  else{
    Tdata.rows(0,cov_num-1) = Data;
    arma::rowvec In(n,arma::fill::ones);
    Data.insert_rows(0,In);
    arma::vec T_One(1);
    arma::mat M(cov_num+1,cov_num+1);
    double pp;
    T_One = arma::randi<arma::vec>(1,arma::distr_param(0,1));
    assig(0) = T_One(0);
    assig(1) = 1 - T_One(0);
    arma::mat d(1,1);
    arma::vec x(cov_num);
    arma::vec b = (2*assig(0)-1)*Data.col(0)+(2*assig(1)-1)*Data.col(1);
    for(int i =2; i< n;i++){
      T_One = arma::randu<arma::vec>(1);
      M = pinv(Data.cols(0,i-1)*(Data.cols(0,i-1).t()));
      d = Data.col(i).t()*M*b;
      pp = pow(1.0-d(0,0),2)/(pow(1.0-d(0,0),2)+pow(1.0+d(0,0),2));
      assig(i) = sum(T_One<pp);
      b = b + (2*assig(i)-1)*Data.col(i);
    }
    Tdata.row(cov_num) = assig+1;
    return Tdata;
  }
}


//[[Rcpp::export]]
arma::mat DoptBCD_getData(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          double mu1,double mu2,double sigma){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    Tdata = DoptBCD_In(n,cov_num,level_num,pr);
    arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
    if(type == "logit"){
      arma::vec m = exp(yita)/(1+exp(yita));
      arma::vec y = arma::randu(n);
      y.elem(find(y>=m)).fill(2);
      y.elem(find(y<m)).fill(1);
      y.replace(2,0);
      Tdata.insert_rows(cov_num+1,y.t());
    }
    if(type == "linear"){
      if(sigma_check(sigma) == TRUE){
        arma::vec epsilon(n,arma::fill::randn);
        epsilon = epsilon*sigma;
        arma::vec y = yita + epsilon;
        Tdata.insert_rows(cov_num+1,y.t());
      }
      else{
        Tdata.zeros();
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}

arma::rowvec DoptBCDOne(arma::mat Data,int n,int cov_num){
  arma::rowvec assig(n);
  arma::rowvec In(n,arma::fill::ones);
  Data.insert_rows(0,In);
  arma::vec T_One(1);
  arma::mat M(cov_num+1,cov_num+1);
  double pp;
  T_One = arma::randi<arma::vec>(1,arma::distr_param(0,1));
  assig(0) = T_One(0);
  assig(1) = 1 - T_One(0);
  arma::mat d(1,1);
  arma::vec x(cov_num);
  arma::vec b = (2*assig(0)-1)*Data.col(0)+(2*assig(1)-1)*Data.col(1);
  for(int i =2; i< n;i++){
    T_One = arma::randu<arma::vec>(1);
    M = pinv(Data.cols(0,i-1)*(Data.cols(0,i-1).t()));
    d = Data.col(i).t()*M*b;
    pp = pow(1.0-d(0,0),2)/(pow(1.0-d(0,0),2)+pow(1.0+d(0,0),2));
    assig(i) = sum(T_One<pp);
    b = b + (2*assig(i)-1)*Data.col(i);
  }
  return assig+1;
}

//[[Rcpp::export]]
Rcpp::List DoptBCD_RT(DataFrame data,double Reps){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  arma::vec diff(Reps);
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  for(int i = 0; i < Reps; i++){
    assignew = DoptBCDOne(data_proc.rows(0,cov_num-1),n,cov_num);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = Reps,
                              Named("pval") = 0,
                              Named("estimate") = diff_data);
  }
  else{
    double index = min(find(diff>=diff_data));
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = index,
                              Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                     Named("estimate") = diff_data);
  }
}

//[[Rcpp::export]]
arma::vec DoptBCD_BT(DataFrame data,double B){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  int n = data_proc.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  arma::rowvec assignew(n,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::mat dtemp = data_proc.rows(0,cov_num-1);
  arma::rowvec y = data_proc.row(cov_num+1);
  for(int i = 0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    assignew = DoptBCDOne(dtemp.cols(shuffle),n,cov_num);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec r(4);
  r(0) = s;
  r(1) = vb;
  r(2) = s/vb;
  r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  return r;
}


//[[Rcpp::export]]
double DoptBCD_RT_In(arma::mat data,double Reps){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  int n = data.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  arma::vec diff(Reps);
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  for(int i = 0; i < Reps; i++){
    assignew = DoptBCDOne(data.rows(0,cov_num-1),n,cov_num);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return 0;
  }
  else{
    double index = min(find(diff>=diff_data));
    return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
  }
}

//[[Rcpp::export]]
double DoptBCD_BT_In(arma::mat data,double B){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  int n = data.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  arma::rowvec assignew(n,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::mat dtemp = data.rows(0,cov_num-1);
  arma::rowvec y = data.row(cov_num+1);
  for(int i = 0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    assignew = DoptBCDOne(dtemp.cols(shuffle),n,cov_num);
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
}

// [[Rcpp::export]]
arma::vec DoptBCD_RT_power(int n,unsigned int cov_num,arma::vec level_num,
                           arma::vec pr,std::string type,arma::vec beta,
                           arma::vec mu1,arma::vec mu2,double sigma,
                           double Iternum,double sl,double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = DoptBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma);
        double pval = DoptBCD_RT_In(data,Reps);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}


// [[Rcpp::export]]
arma::vec DoptBCD_BT_power(int n,unsigned int cov_num,arma::vec level_num,
                           arma::vec pr,std::string type,arma::vec beta,
                           arma::vec mu1,arma::vec mu2,double sigma,
                           double Iternum,double sl,double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = DoptBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma);
        double pval = DoptBCD_BT_In(data,B);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

double ppc(int x,double a){
  if(x>-0.1&& x<0.1){
    return 0.5;
  }
  else if(x > 0.1){
    return 1.0/(pow(x,a)+1.0);
  }
  else{
    return 1.0-1.0/(pow(x,a)+1.0);
  }
}

// [[Rcpp::export]]
arma::mat AdjBCD_In(int n, unsigned int cov_num, arma::vec level_num, arma::vec pr, 
                    double a = 2){
  arma::mat Tdata(cov_num+1,n);
  if(level_num.n_elem != cov_num){
    Rcpp::Rcout<<
      "Error in level_num: default of level_num is (2,  2). Length of level_num should match cov_num!"
      <<std::endl;
      return Tdata.zeros();
  }
  arma::vec val = level_num(find(level_num <= 1.999999));
  if(val.n_elem > 0){
    Rcpp::Rcout<<
      "Error in level_num: number of levels of each covariate must be larger than 2!"
      <<std::endl;
      return Tdata.zeros();
  }
  if(pr.n_elem != sum(level_num)){
    Rcpp::Rcout<<
      "Error in pr: default of pr is (0.5, 0.5, 0.5, 0.5). Length of pr should correspond to number of all levels for all margin!"
      <<std::endl;
      return Tdata.zeros();
  }
  arma::mat ProbS = Prob_S(cov_num, level_num, pr);
  arma::rowvec valv = sum(ProbS, 0);
  arma::vec valvCh = valv(find(valv <=0.999999 || valv >= 1.000001));
  if(valvCh.n_elem > 0){
    Rcpp::Rcout<<"Error in pr: probabilities of each margin must sum up to 1!"<<std::endl;
    return Tdata.zeros();
  }
  else{
    arma::mat Data = genData_sim(n,cov_num,level_num,ProbS);
    arma::mat P = PStrR(Data);
    int strt_num = P.n_cols;
    
    arma::vec D(strt_num,arma::fill::zeros);
    arma::rowvec assig(n); 
    int r;
    double pp;
    arma::vec T_One(1);
    arma::vec brid(2); brid(0) = 1; brid(1) = -1;
    for(int i=0;i<n;i++){
      r = ReturnCol(P,Data.col(i))(0);
      pp = ppc(D(r-1),a);
      T_One = arma::randu<arma::vec>(1);
      assig(i) = sum(T_One>pp) + 1;
      D(r-1) = D(r-1)+brid(sum(T_One>pp));
    }
    Tdata.rows(0,cov_num-1) = Data;
    Tdata.row(cov_num) = assig;
    return Tdata;
  }
}

//[[Rcpp::export]]
arma::mat AdjBCD_getData(int n,unsigned int cov_num,arma::vec level_num,
                         arma::vec pr,std::string type,arma::vec beta,
                         double mu1,double mu2,double sigma,double a){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    Tdata = AdjBCD_In(n,cov_num,level_num,pr,a);
    arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
    if(type == "logit"){
      arma::vec m = exp(yita)/(1+exp(yita));
      arma::vec y = arma::randu(n);
      y.elem(find(y>=m)).fill(2);
      y.elem(find(y<m)).fill(1);
      y.replace(2,0);
      Tdata.insert_rows(cov_num+1,y.t());
    }
    if(type == "linear"){
      if(sigma_check(sigma) == TRUE){
        arma::vec epsilon(n,arma::fill::randn);
        epsilon = epsilon*sigma;
        arma::vec y = yita + epsilon;
        Tdata.insert_rows(cov_num+1,y.t());
      }
      else{
        Tdata.zeros();
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}

arma::vec AdjBCDOne(arma::vec D,arma::mat P,arma::vec cov_profile,double a){
  arma::vec brid(2); brid(0) = 1; brid(1) = -1;
  int r = ReturnCol(P,cov_profile)(0);
  double pp = ppc(D(r-1),a);
  arma::vec T_One = arma::randu<arma::vec>(1);
  D(r-1) = D(r-1)+brid(sum(T_One>pp));
  D(D.n_elem-1) = sum(T_One>pp) + 1;
  return D;
}

//[[Rcpp::export]]
Rcpp::List AdjBCD_RT(DataFrame data,double Reps,double a){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec D(1 + strt_num + sum(level_num));
  for(int i = 0; i < Reps; i++){
    D.zeros();
    for(int j = 0; j < n; j++){
      D = AdjBCDOne(D,P,data_proc.col(j).head(cov_num),a);
      assignew(j) = D(D.n_elem-1);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = Reps,
                              Named("pval") = 0,
                              Named("estimate") = diff_data,
                              Named("D") = D);
  }
  else{
    double index = min(find(diff>=diff_data));
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = index,
                              Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                     Named("estimate") = diff_data,
                                                     Named("D") = D);
  }
}


//[[Rcpp::export]]
arma::vec AdjBCD_BT(DataFrame data,double B,double a){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  double n = data_proc.n_cols;
  int m = data_proc.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data_proc.row(cov_num+1);
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec D(2 + strt_num + sum(level_num));
  for(int i = 0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(m,arma::distr_param(0,m-1));
    ynew = trans(y(shuffle));
    D.zeros();
    for(int j = 0; j < n; j++){
      D = AdjBCDOne(D,P,data_proc.col(shuffle(j)).head(cov_num),a);
      assignew(j) = D(D.n_elem-1);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec r(4);
  r(0) = s;
  r(1) = vb;
  r(2) = s/vb;
  r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  return r;
}

//[[Rcpp::export]]
double AdjBCD_RT_In(arma::mat data,double Reps,double a){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  arma::vec D(2 + strt_num + sum(level_num));
  for(int i = 0; i < Reps; i++){
    D.zeros();
    for(int j = 0; j < n; j++){
      D = AdjBCDOne(D,P,data.col(j).head(cov_num),a);
      assignew(j) = D(D.n_elem-1);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return 0;
  }
  else{
    double index = min(find(diff>=diff_data));
    return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
  }
}

//[[Rcpp::export]]
double AdjBCD_BT_In(arma::mat data,double B,double a){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  double n = data.n_cols;
  int m = data.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data.row(cov_num+1);
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec D(2 + strt_num + sum(level_num));
  for(int i = 0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(m,arma::distr_param(0,m-1));
    ynew = trans(y(shuffle));
    D.zeros();
    for(int j = 0; j < n; j++){
      D = AdjBCDOne(D,P,data.col(shuffle(j)).head(cov_num),a);
      assignew(j) = D(D.n_elem-1);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  arma::vec r(4);
  return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
}

// [[Rcpp::export]]
arma::vec AdjBCD_RT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double Iternum,double sl,double a,
                          double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)  
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = AdjBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,a);
        double pval = AdjBCD_RT_In(data,Reps,a);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec AdjBCD_BT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double Iternum,double sl,double a,double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)   
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = AdjBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,a);
        double pval = AdjBCD_RT_In(data,B,a);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

//[[Rcpp::export]]
arma::mat StrPBR_getData(int n,unsigned int cov_num,arma::vec level_num,
                         arma::vec pr,std::string type,arma::vec beta,
                         double mu1,double mu2,double sigma,int bsize){
  bool check_beta = beta_check(cov_num,beta);
  arma::mat Tdata(cov_num+1,n);
  if(check_beta == TRUE){
    arma::mat ProbS = Prob_S(cov_num,level_num,pr);
    Tdata.rows(0,cov_num-1) = genData_sim(n,cov_num,level_num,ProbS);
    arma::vec level_num_data = max(Tdata.rows(0,cov_num-1), 1);
    arma::mat P = PStrR(Tdata.rows(0,cov_num-1));
    int strt_num = P.n_cols;
    arma::vec D(1 + strt_num + sum(level_num_data),arma::fill::zeros);
    arma::mat B = Bpert(bsize,2);
    int Bsize = B.n_cols;
    arma::uvec shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
    arma::mat BG = B.cols(shuffle);
    arma::vec strp(strt_num,arma::fill::zeros);
    arma::field<arma::mat> Res(4,1);
    for(int i=0; i < n; i++){
      Res = StrROne(D,P,Tdata.col(i).head(cov_num),cov_num,level_num_data,bsize,B,BG,strp);
      strp = Res(0,0);
      BG = Res(1,0);
      Tdata(cov_num,i) = Res(2,0)(0,0);
      D = Res(3,0);
    }
    arma::vec yita = (Tdata.rows(0,cov_num-1)).t()*beta+(Tdata.row(cov_num)).t()*(mu2-mu1)+2*mu1-mu2;
    if(type == "logit"){
      arma::vec m = exp(yita)/(1+exp(yita));
      arma::vec y = arma::randu(n);
      y.elem(find(y>=m)).fill(2);
      y.elem(find(y<m)).fill(1);
      y.replace(2,0);
      Tdata.insert_rows(cov_num+1,y.t());
    }
    if(type == "linear"){
      if(sigma_check(sigma) == TRUE){
        arma::vec epsilon(n,arma::fill::randn);
        epsilon = epsilon*sigma;
        arma::vec y = yita + epsilon;
        Tdata.insert_rows(cov_num+1,y.t());
      }
      else{
        Tdata.zeros();
      }
    }
  }
  else{
    Tdata.zeros();
  }
  return Tdata;
}

//[[Rcpp::export]]
Rcpp::List StrPBR_RT(DataFrame data,double Reps,int bsize){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data_proc.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(1 + strt_num + sum(level_num),arma::fill::zeros);
  arma::mat B = Bpert(bsize,2);
  int Bsize = B.n_cols;
  arma::uvec shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
  arma::mat BG = B.cols(shuffle);
  arma::vec strp(strt_num,arma::fill::zeros);
  arma::field<arma::mat> Res(4,1);
  double diff_data = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  for(int i=0; i < Reps; i++){
    D.zeros();
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,B.n_cols-1));
    BG = B.cols(shuffle);
    strp.zeros();
    for(int j=0; j < n; j++){
      Res = StrROne(D,P,data_proc.col(j).head(cov_num),cov_num,level_num,bsize,B,BG,strp);
      strp = Res(0,0);
      BG = Res(1,0);
      assignew(j) = Res(2,0)(0,0);
      D = Res(3,0);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data_proc.row(cov_num+1)%(assignew-2))/n1-sum(data_proc.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = Reps,
                              Named("pval") = 0,
                              Named("estimate") = diff_data);
  }
  else{
    double index = min(find(diff>=diff_data));
    return Rcpp::List::create(Named("Randata") = diff,
                              Named("index") = index,
                              Named("pval") = index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps),
                                                     Named("estimate") = diff_data);
  }
}

//[[Rcpp::export]]
arma::vec StrPBR_BT(DataFrame data,double B,int bsize){
  Rcpp::List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  arma::vec level_num = resu["level_num"];
  arma::mat P = PStrR(data_proc.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  int n = data_proc.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data_proc.row(cov_num+1);
  arma::vec D(1 + strt_num + sum(level_num),arma::fill::zeros);
  arma::mat Bp = Bpert(bsize,2);
  int Bsize = Bp.n_cols;
  arma::uvec shuffleB = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
  arma::mat BG = Bp.cols(shuffleB);
  arma::vec strp(strt_num,arma::fill::zeros);
  arma::field<arma::mat> Res(4,1);
  arma::rowvec assignew(n,arma::fill::zeros);
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    shuffleB = arma::randi<arma::uvec>(n,arma::distr_param(0,Bp.n_cols-1));
    BG = Bp.cols(shuffleB);
    strp.zeros();
    for(int j=0; j < n; j++){
      Res = StrROne(D,P,data_proc.col(shuffle(j)).head(cov_num),cov_num,level_num,bsize,Bp,BG,strp);
      strp = Res(0,0);
      BG = Res(1,0);
      assignew(j) = Res(2,0)(0,0);
      D = Res(3,0);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data_proc.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-2))/n1-sum(data_proc.row(cov_num+1)%(data_proc.row(cov_num)-1))/n0;
  arma::vec r(4);
  r(0) = s;
  r(1) = vb;
  r(2) = s/vb;
  r(3) = arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
  return r;
}

//[[Rcpp::export]]
double StrPBR_RT_In(arma::mat data,double Reps,int bsize){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data.n_cols;
  arma::rowvec assignew(n,arma::fill::zeros);
  arma::vec diff(Reps,arma::fill::zeros);
  double n0,n1;
  n1 = -sum(data.row(cov_num)-2);
  n0 = n-n1;
  arma::vec D(1 + strt_num + sum(level_num),arma::fill::zeros);
  arma::mat B = Bpert(bsize,2);
  int Bsize = B.n_cols;
  arma::uvec shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
  arma::mat BG = B.cols(shuffle);
  arma::vec strp(strt_num,arma::fill::zeros);
  arma::field<arma::mat> Res(4,1);
  double diff_data = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  for(int i=0; i < Reps; i++){
    D.zeros();
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,B.n_cols-1));
    BG = B.cols(shuffle);
    strp.zeros();
    for(int j=0; j < n; j++){
      Res = StrROne(D,P,data.col(j).head(cov_num),cov_num,level_num,bsize,B,BG,strp);
      strp = Res(0,0);
      BG = Res(1,0);
      assignew(j) = Res(2,0)(0,0);
      D = Res(3,0);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    diff(i) = -sum(data.row(cov_num+1)%(assignew-2))/n1-sum(data.row(cov_num+1)%(assignew-1))/n0;
  }
  diff = sort(diff);
  arma::uvec fi = find(diff>=diff_data);
  if(fi.is_empty()){
    return 0;
  }
  else{
    double index = min(find(diff>=diff_data));
    return index/Reps<(1-index/Reps)?index/Reps:(1-index/Reps);
  }
}

//[[Rcpp::export]]
double StrPBR_BT_In(arma::mat data,double B,int bsize){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  arma::mat P = PStrR(data.rows(0,cov_num-1));
  int strt_num = P.n_cols;
  arma::uvec a0(1,arma::fill::zeros);
  int n = data.n_cols;
  arma::uvec shuffle(n,arma::fill::zeros);
  arma::vec theta(B,arma::fill::zeros);
  double n0,n1;
  arma::rowvec ynew(n);
  arma::rowvec y = data.row(cov_num+1);
  arma::vec D(1 + strt_num + sum(level_num),arma::fill::zeros);
  arma::mat Bp = Bpert(bsize,2);
  int Bsize = Bp.n_cols;
  arma::uvec shuffleB = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
  arma::mat BG = Bp.cols(shuffleB);
  arma::vec strp(strt_num,arma::fill::zeros);
  arma::field<arma::mat> Res(4,1);
  arma::rowvec assignew(n,arma::fill::zeros);
  for(int i=0; i < B; i++){
    shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,n-1));
    ynew = trans(y(shuffle));
    D.zeros();
    shuffleB = arma::randi<arma::uvec>(n,arma::distr_param(0,Bp.n_cols-1));
    BG = Bp.cols(shuffleB);
    strp.zeros();
    for(int j=0; j < n; j++){
      Res = StrROne(D,P,data.col(shuffle(j)).head(cov_num),cov_num,level_num,bsize,Bp,BG,strp);
      strp = Res(0,0);
      BG = Res(1,0);
      assignew(j) = Res(2,0)(0,0);
      D = Res(3,0);
    }
    n1 = -sum(assignew-2);
    n0 = n - n1;
    theta(i)= -sum(ynew%(assignew-2))/n1-sum(ynew%(assignew-1))/n0;
  }
  n1 = -sum(data.row(cov_num)-2);
  n0 = n - n1;
  double vb = stddev(theta);
  double s = -sum(data.row(cov_num+1)%(data.row(cov_num)-2))/n1-sum(data.row(cov_num+1)%(data.row(cov_num)-1))/n0;
  return arma::normcdf(s/vb)<(1-arma::normcdf(s/vb))?arma::normcdf(s/vb):(1-arma::normcdf(s/vb));
}


// [[Rcpp::export]]
arma::vec StrPBR_RT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double Iternum,double sl,int bsize,
                          double Reps,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)   
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrPBR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,bsize);
        double pval = StrPBR_RT_In(data,Reps,bsize);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec StrPBR_BT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double Iternum,double sl,int bsize,
                          double B,int nthreads){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    
#pragma omp parallel for shared(n,N,mu1,mu2) private(i,j) num_threads(nthreads) schedule(static,N1*N/nthreads) collapse(2)    
    for(i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrPBR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,bsize);
        double pval = StrPBR_BT_In(data,B,bsize);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
#pragma omp barrier
    
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

//[[Rcpp::export]]
arma::vec CTT(DataFrame data){
  List resu = Preprocess_out(data);
  arma::mat data_proc = resu["data"];
  unsigned int cov_num = resu["cov_num"];
  int n = data_proc.n_cols;
  arma::rowvec T = data_proc.row(cov_num);
  T.elem(find(T!=1)).zeros();
  arma::mat G(cov_num+2,n);
  G.row(0) = T;
  G.row(1) = 1-T;
  G.rows(2,cov_num+1) = data_proc.rows(0,cov_num-1);
  arma::mat Gi = pinv(G*G.t());
  arma::vec thetahat = Gi*G*data_proc.row(cov_num+1).t();
  arma::vec L(cov_num+2,arma::fill::zeros);
  L(0) = 1;
  L(1) = -1;
  arma::mat a = L.t()*Gi*L;
  arma::vec r(4);
  r(0) = sum(L%thetahat);
  r(1) = sqrt(sum(pow((data_proc.row(cov_num+1).t()-G.t()*thetahat),2))/(n-cov_num-2)*a(0,0));
  r(2) = r(0)/r(1);
  r(3) = r(2)>0?(1-arma::normcdf(r(2))):arma::normcdf(r(2));
  return r;
}

//[[Rcpp::export]]
double CTT_In(arma::mat data){
  unsigned int cov_num = data.n_rows-2;
  arma::vec level_num = max(data.rows(0,cov_num-1), 1);
  int n = data.n_cols;
  arma::rowvec T = data.row(cov_num);
  T.elem(find(T!=1)).zeros();
  arma::mat G(cov_num+2,n);
  G.row(0) = T;
  G.row(1) = 1-T;
  G.rows(2,cov_num+1) = data.rows(0,cov_num-1);
  arma::mat Gi = pinv(G*G.t());
  arma::vec thetahat = Gi*G*data.row(cov_num+1).t();
  arma::vec L(cov_num+2,arma::fill::zeros);
  L(0) = 1;
  L(1) = -1;
  arma::mat a = L.t()*Gi*L;
  double t = sum(L%thetahat)/(sqrt(sum(pow((data.row(cov_num+1).t()-G.t()*thetahat),2))/(n-cov_num-2)*a(0,0)));
  return t>0?(1-arma::normcdf(t)):arma::normcdf(t);
}

// [[Rcpp::export]]
arma::vec HuHuCAR_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                           arma::vec pr,std::string type,arma::vec beta,
                           arma::vec mu1, arma::vec mu2,double sigma,
                           double Iternum,double sl,arma::vec omega,double p){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = HuHuCAR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,omega,p);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec PocSimMIN_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                             arma::vec pr,std::string type,arma::vec beta,
                             arma::vec mu1, arma::vec mu2,double sigma,
                             double Iternum,double sl,arma::vec weight,double p){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = PocSimMIN_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,weight,p);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec StrBCD_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1, arma::vec mu2,double sigma,
                          double Iternum,double sl,double p){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,p);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec DoptBCD_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                           arma::vec pr,std::string type,arma::vec beta,
                           arma::vec mu1, arma::vec mu2,double sigma,
                           int Iternum,double sl){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = DoptBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec AdjBCD_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          int Iternum,double sl,double a){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = AdjBCD_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,a);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}

// [[Rcpp::export]]
arma::vec StrPBR_CT_power(int n,unsigned int cov_num,arma::vec level_num,
                          arma::vec pr,std::string type,arma::vec beta,
                          arma::vec mu1,arma::vec mu2,double sigma,
                          double Iternum,double sl,int bsize){
  if(mu1.n_elem != mu2.n_elem){
    arma::vec result(mu1.n_elem);
    Rcpp::Rcout<<"The length of two mu's must match!"<<std::endl;
    return result;
  }
  else{
    unsigned int N = Iternum;
    unsigned int N1 = mu1.n_elem; 
    unsigned int i,j;
    arma::mat p_all(N1,N,arma::fill::zeros);
    for (i = 0; i < N1; i++){
      for(j = 0; j < N; j++){
        arma::mat data = StrPBR_getData(n,cov_num,level_num,pr,type,beta,mu1(i),mu2(i),sigma,bsize);
        double pval = CTT_In(data);
        p_all(i,j) = pval<sl/2?1:0;
      }
    }
    arma::vec result(2*N1);
    for(i = 0; i < N1; i++){
      result(i) = sum(p_all.row(i))/Iternum;
      result(i+N1) = stddev(p_all.row(i));
    }
    return result;
  }
}
