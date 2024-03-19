// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <R.h>
#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <RcppEigen.h>
#include <random>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export]]
std::normal_distribution<double> normalDistribution(0.0, 1.0);
std::random_device rd{};
std::mt19937 gen{rd()};

// Create matrix of given dimensions with normally distributed entries
// [[Rcpp::export]]
MatrixXd normalRandomMatrix(int nrows, int ncols) {
  MatrixXd temp = MatrixXd::Zero(nrows,ncols);
  for(int i=0;i<nrows;i++) {
    for(int j=0;j<ncols;j++) {
      temp(i,j) = normalDistribution(gen);
    }
  }
  return temp;
}

// Genomic relationship matrix
// [[Rcpp::export]]
MatrixXd grmCpp(MatrixXd X) {
  VectorXd pv = 0.5*X.rowwise().mean();
  VectorXd qv = 1.0-pv.array();
  VectorXd rv = 2.0*pv.cwiseProduct(qv);
  X = X.colwise()-2.0*pv;
  double n = 1.0/rv.sum();
  return n*X.transpose()*X;
}

// fast SVD algorithm of Halko et al. (2009)
// [[Rcpp::export]]
MatrixXd fastSVDCpp(MatrixXd X, int k, int q=2) {
  
  int nrows = X.rows();
  int ncols = k*q;
  MatrixXd Y = normalRandomMatrix(nrows, ncols);
  Y = X*Y;
  for(int i = 1; i <= q; i++){
    Y = X*X.transpose()*Y;
  }
  
  HouseholderQR<MatrixXd> qr(Y);
  MatrixXd Q = qr.householderQ();
  
  MatrixXd B = Q.transpose()*X;
  
  JacobiSVD<MatrixXd> svd(B,ComputeFullU);
  
  MatrixXd thinU = MatrixXd::Identity(Y.rows(), k);
  MatrixXd U = svd.matrixU()*thinU;
  
  return Q*U;
  
}

// fast GRM eigenvectors
// [[Rcpp::export]]
MatrixXd fastGRMCpp(MatrixXd X, int k, int q=2) {
  
  X = grmCpp(X);
  
  MatrixXd V = fastSVDCpp(X, k, q);
  
  return V;
  
}



// [[Rcpp::export]]
NumericVector Rcpp_seq(double start, double end, double bins) {
  
  NumericVector sequence(bins + 1);
  
  int by = round((end - start + 1)/bins);
  
    
  int element = start;
  int i = 0;
  while( element < end ){
    sequence[i] = element;
    element = element + by;
    i = i + 1;
  }

  sequence[i] = end + 1;

  return sequence;
}

// // [[Rcpp::export]]
// DataFrame trimVCF(DataFrame vcf, double w = 1, double n = 100){
// 
//   NumericVector s = Rcpp_seq(1, vcf.nrow(), n);
//   int low = int(s[w - 1]);
//   int high = int(s[w])-1;
//   DataFrame w_gt_table;
//   for(int col = 8; col < vcf.size(); col++){
//     CharacterVector v = vcf[col];
//     w_gt_table.push_back(v[Rcpp::Range(low - 1, high - 1)]);
//   }
// 
//   return w_gt_table;
// 
// }
// 
// DataFrame get_allReaddepth(DataFrame vcf){
// 
//   DataFrame allele_depth;
// 
//   for(int variant = 0; variant < nrow(vcf); variant++) {
//     gt_pos = grep('GT',strsplit(vcf[variant,1], ':')[[1]]);
//     ad_pos = grep('AD',strsplit(vcf[variant,1], ':')[[1]]);
//     temp_gts = vcf[variant,-1]
// 
//     for(sample in 1:length(temp_gts)){
//       gt = strsplit(strsplit(as.character(temp_gts[sample]),':')[[1]][gt_pos], '/')[[1]]
//       ad = strsplit(strsplit(as.character(temp_gts[sample]),':')[[1]][ad_pos], ',')[[1]]
// 
//       gt = unique(gt)
// 
//       if(gt[1] == '.'){
//         ad_df = data.frame(vcf_pos = variant + low - 1,
//                            sample = names(temp_gts)[sample],
//                                                    allele = NA,
//                                                    allele_depth = NA,
//                                                    ap_to_major = NA,
//                                                    typeof_gt = NA,
//                                                    typeof_all = NA)
//       }else{
// 
//         ad_df = data.frame(vcf_pos = variant + low - 1,
//                            sample = names(temp_gts)[sample],
//                                                    allele = as.character(0:(length(ad)-1)),
//                                                    allele_depth = as.numeric(ad))
// 
//         ad_df = ad_df[ad_df$allele_depth != 0,]
// 
//         if(nrow(ad_df) == 0){
// 
//           ad_df = data.frame(vcf_pos = variant + low - 1,
//                              sample = names(temp_gts)[sample],
//                                                      allele = 'uninformative alleles',
//                                                      allele_depth = 0,
//                                                      ap_to_major = NA,
//                                                      typeof_gt = NA,
//                                                      typeof_all = NA)
// 
//         }else{
// 
//           ad_df = ad_df[order(ad_df$allele_depth, decreasing = T),]
// 
//           ad_df$ap_to_major = ad_df$allele_depth/max(ad_df$allele_depth)
// 
//           ad_df$typeof_gt = ifelse(length(gt) == 1, 'Homozygous', 'Heterozygous')
// 
//           ad_df$typeof_all = NA
// 
//           if(length(gt) == 1){
// 
//             ad_df[1, ][['typeof_all']] = 'Homozygous'
// 
//             if(nrow(ad_df[-1,]) > 0 ){
// 
//               ad_df[-1,][['typeof_all']] = 'Excluded minor allele'
// 
//             }
// 
//           }else if(length(gt) == 2){
// 
//             ad_df[1, ][['typeof_all']] = 'Major allele'
//             ad_df[2, ][['typeof_all']] = 'Minor allele'
// 
//             if(nrow(ad_df[-1:-2,]) > 0 ){
// 
//               ad_df[-1:-2,][['typeof_all']] = 'Excluded minor allele'
// 
//             }
//           }
// 
//         }
// 
//       }
// 
//       allele_depth = rbind(allele_depth, ad_df)
// 
//     }
//   }
// 
//   return(allele_depth)
// 
// }