#include "Data2.h"
#include "globals.h"


arma::field<arma::umat> Data2::get_zt() const
{
  arma::field<arma::umat> czt(zt.size());
  for(int i = 0; i < zt.size(); i++) {
    czt(i) = Rcpp::as<arma::umat>(zt[i]);
  }
  return czt;
}

arma::umat Data2::get_zy() const
{
  arma::umat mzy = Rcpp::as<arma::umat>(zy);
  return mzy;
}

arma::uvec Data2::get_e() const
{
  arma::uvec ve = Rcpp::as<arma::uvec>(E);
  return ve;
}

arma::vec Data2::get_Y() const
{
  arma::vec vy = Rcpp::as<arma::vec>(Y);
  return vy;
}
#include "Data.h"

Data::Data(
  Rcpp::IntegerMatrix mat1Z,
  Rcpp::NumericMatrix mat1f,
  Rcpp::List mat2Z,
  Rcpp::IntegerMatrix range0)
{
  this->mat1Z = mat1Z;
  this->mat1f = mat1f;
  // List mat2(10);
  // for(int i = 0; i< 10; i++)
  // {
  //   mat2[i] = as<arma::umat>(mat2Z[i]);
  // }
  this->mat2Z = mat2Z;
  this->range0 = range0;
}

// [[Rcpp::export]]
arma::field<arma::mat> make_mat2(arma::vec& tk0,
				 arma::vec& Y0,
				 arma::uvec& id0,
				 arma::mat& z0) {
  int K = tk0.size();
  arma::field<arma::mat> out(K);
  for(int k = 0; k < K; k++) {
    arma::uvec ind = arma::find( Y0 >= tk0[k] );
    arma::uvec idAR = id0(ind);
    arma::mat subX = z0.rows(ind);
    arma::uvec ind2 = arma::find(Rcpp::as<arma::uvec>(Rcpp::duplicated(Rcpp::IntegerVector(idAR.begin(), idAR.end()) )) == 0);
    out[k] = subX.rows(ind2);
  }
  return out;
}

// [[Rcpp::export]]
arma::field<arma::mat> make_mat2_t(arma::vec& tk0,
				   arma::vec& Y0,
				   arma::uvec& id0,
				   arma::mat& z0) {
  int K = tk0.size();
  arma::field<arma::mat> out(K);
  for(int k = 0; k < K; k++) {
    arma::uvec ind = arma::find( Y0 >= tk0[k] );
    arma::uvec idAR = id0(ind);
    arma::mat subX = z0.rows(ind);
    arma::uvec ind2 = arma::find(Rcpp::as<arma::uvec>(Rcpp::duplicated(Rcpp::IntegerVector(idAR.begin(), idAR.end()) )) == 0);
    out[k] = subX.rows(ind2).t();
  }
  return out;
}


Rcpp::IntegerMatrix Data::get_mat1Z() const {
  return mat1Z;
}

Rcpp::NumericMatrix Data::get_mat1f() const {
  return mat1f;
}

Rcpp::IntegerMatrix Data::get_range0() const {
  return range0;
}

Rcpp::List Data::get_mat2ZK() const {
  return mat2Z;
}

Rcpp::IntegerMatrix Data::get_mat2Zk( int k ) const {
  return mat2Z[k];
}

arma::uword Data::get_len_mat2k( int k ) const {
  arma::umat m = mat2Z[k];
  return m.n_rows;
}
#include "Forest.h"
#include "globals.h"

int Forest::trainRF(std::vector<std::shared_ptr<Tree> >& trees,
                    const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf,
                    const arma::umat& range0,
                    const arma::umat& ids,
                    const arma::uvec& e) {
  int n = mat1Z.n_rows;
  for(size_t i = 0; i != NUM_TREE; i++) {
    arma::field<arma::umat>  mat2Zsub(K);
    for(size_t k = 0; k != K; k++) {
      int nr = mat2Zf(k).n_rows;
      mat2Zsub(k) = mat2Zf(k).rows( ids( ids.n_rows*i + find(ids.col(i) >= n-nr) ) + nr - n );
    }
    //arma::umat mat1Zsub = mat1Z.rows( ids.col(i) );
    //arma::mat mat1fsub = mat1f.cols( ids.col(i) );
    trees.push_back( train( mat1Z.rows( ids.col(i) ),
                            mat1f.cols( ids.col(i) ),
                            mat2Zsub, range0,
                            e( ids.col(i) )) );
  }
  return 1;
}

std::shared_ptr<Tree> Forest::train(const arma::umat& mat1Z,
                                    const arma::mat& mat1f,
                                    const arma::field<arma::umat>& mat2Zf,
                                    const arma::umat& range0,
                                    const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = range0.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  if(spCriterion == 1) {
    size_t ndcount = 0;
    size_t countsp = 0;
    int end = 0;
    while(end == 0) {
      end = split_DICON(mat1Z, mat1f, mat2Zf,
                        left_childs, right_childs,
                        split_vars, split_values, isLeaf,
                        ranges, nodeSampleY, nodeSample,
                        countsp, ndcount, e);
      if(ndcount >= MAX_NODE - 2) {
        isLeaf( arma::find(left_childs == 0) ).ones();
        break;
      }
    }
    arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
    std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty), right_childs(nonEmpty),
                                      split_vars(nonEmpty), split_values(nonEmpty), isLeaf(nonEmpty)));
    return tr;
  } else {
    arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
    size_t ndcount = 0;
    size_t countsp = 0;
    int end = 0;
    while(end == 0) {
      end = split_ICON(mat1Z, mat1f, mat2Zf,
                       left_childs, right_childs,
                       split_vars, split_values, isLeaf,
                       fmat, Smat, 
                       ranges, nodeSampleY, nodeSample,
                       countsp, ndcount, e);
      if(ndcount >= MAX_NODE - 2) {
        isLeaf( arma::find(left_childs == 0) ).ones();
        break;
      }
    }
    arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
    //Rcpp::Rcout << "2";
    std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                      right_childs(nonEmpty),
                                      split_vars(nonEmpty),
                                      split_values(nonEmpty),
                                      isLeaf(nonEmpty)));
    return tr;
  }
}


int Forest::split_DICON(const arma::umat& mat1Z,
                        const arma::mat& mat1f,
                        const arma::field<arma::umat>& mat2Zf, // dat
                        arma::uvec& left_childs,
                        arma::uvec& right_childs,
                        arma::uvec& split_vars,
                        arma::uvec& split_values,
                        arma::uvec& isLeaf, // tree
                        arma::ucube& ranges,
                        arma::field<arma::uvec>& nodeSampleY,
                        arma::field<arma::uvec>& nodeSample,
                        size_t& countsp,
                        size_t& ndcount,
                        const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  while(varsp == -1 && countsp <= ndcount) {
    nd = countsp;
    arma::ivec bestSp = find_split_DICON(nd,
                                         mat1Z, mat1f, mat2Zf,
                                         ranges, nodeSampleY, nodeSample, e);
    varsp = bestSp(1);
    cutsp = bestSp(2);
    if(varsp == -1) {
      isLeaf(nd) = 1;
      while(countsp <= ndcount) {
        countsp++;
        if(isLeaf(countsp) == 0) break;
      }
    }
  }
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <=cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub >cutsp) );
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub<=cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub>cutsp) );
      // if(nodeSample(k, ndc1).size() < MIN_SPLIT1)
      // {isLeaf(ndc1) = 1;}
      // if(nodeSample(k, ndc2).size() < MIN_SPLIT1)
      // {isLeaf(ndc2) = 1;}
    }
    //if( arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 )
    //if( nodeSampleY(ndc1).size() < MIN_SPLIT2 && isLeaf(ndc1) == 1)
    //if( arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 &&  isLeaf(ndc1) )
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) {
      isLeaf(ndc1) = 1;
    } else {
      isLeaf(ndc1) = 0;
    }
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) {
      isLeaf(ndc2) = 1;
    } else {
      isLeaf(ndc2) = 0;
    }
    // if(nodeSample(0, ndc1).size() < MIN_SPLIT1)
    // {isLeaf(ndc1) = 1;}
    // if(nodeSample(0, ndc2).size() < MIN_SPLIT1)
    // {isLeaf(ndc2) = 1;}
    // if( isL1 == 1 && isLeaf(ndc1) == 1)
    // {
    //   isLeaf(ndc1) = 1;
    // }else{
    //   isLeaf(ndc1) = 0;
    // }
    //
    // if( isL2 == 1 && isLeaf(ndc2) == 1)
    // {
    //   isLeaf(ndc2) = 1;
    // }else{
    //   isLeaf(ndc2) = 0;
    // }
    //arma::sum( e(nodeSampleY(ndc1))) < MIN_SPLIT2 &&
    // if( nodeSample(0, ndc1).size() < MIN_SPLIT1 )
    // {
    //   isLeaf(ndc1) = 1;
    // }
    // // arma::sum( e(nodeSampleY(ndc2))) < MIN_SPLIT2 &&
    // if( nodeSample(0, ndc2).size() < MIN_SPLIT1)
    // {
    //   isLeaf(ndc2) = 1;
    // }
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount += 2;
    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0) break;
    }
  } else {
    end = 1;
  }
  return end;
}

arma::ivec Forest::find_split_DICON(size_t nd,
                                    const arma::umat& mat1Z,
                                    const arma::mat& mat1f,
                                    const arma::field<arma::umat>& mat2Zf, // dat
                                    const arma::ucube& ranges,
                                    const arma::field<arma::uvec>& nodeSampleY,
                                    const arma::field<arma::uvec>& nodeSample,
                                    const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::uvec spSet = arma::shuffle( arma::regspace<arma::uvec>(0,P-1) );
  for(auto p : spSet.head(mtry)) {
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    int nel = 0;
    // int nelr = arma::sum(e( indY ) );
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0),
                                                     ranges(nd, p, 1));
     //arma::vec den = (fLSum + fRSum)%(SLSum + SRSum);
     //den( arma::find(den == 0) ).ones();
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        if( mat1Z(indYj, p) == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }
      // if(  (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) )
      //if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) && ( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) ) //|| (SLSum(0) > 9*SRSum(0) || SRSum(0) > 9*SLSum(0))
      //if( (( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) && (SLSum.min() < MIN_NODE1 || SRSum.min() < MIN_NODE1) )  )
      if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        //arma::vec hL = fLSum / SLSum;
        //arma::vec hR = fRSum / SRSum;
        //hL(arma::find(SLSum == 0)).zeros();
        //hR(arma::find(SRSum == 0)).zeros();
        //dICONTemp = abs(arma::sum( (hL - hR) ));
        dICONTemp = arma::sum( abs(fLSum%SRSum - fRSum%SLSum) );
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}


int Forest::split_ICON(const arma::umat& mat1Z,
                       const arma::mat& mat1f,
                       const arma::field<arma::umat>& mat2Zf, // dat
                       arma::uvec& left_childs,
                       arma::uvec& right_childs,
                       arma::uvec& split_vars,
                       arma::uvec& split_values,
                       arma::uvec& isLeaf,
                       arma::mat& fmat,
                       arma::umat& Smat,// tree
                       arma::ucube& ranges,
                       arma::field<arma::uvec>& nodeSampleY,
                       arma::field<arma::uvec>& nodeSample,
                       size_t& countsp,
                       size_t& ndcount,
                       const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  while(varsp == -1 && countsp <= ndcount) {
    nd = countsp;
    arma::ivec bestSp(3);
    bestSp = find_split_ICON(nd,
                             mat1Z, mat1f, mat2Zf, isLeaf,
                             ranges, nodeSampleY, nodeSample,
                             fmat, Smat, ndcount, e);
    varsp = bestSp(1);
    cutsp = bestSp(2);
    if(varsp == -1) {
      isLeaf(nd) = 1;
      while(countsp <= ndcount) {
        countsp++;
        if(isLeaf(countsp) == 0) break;
      }
    }
  }
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <=cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub >cutsp) );
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub<=cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub>cutsp) );
    }
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) isLeaf(ndc1) = 1;
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) isLeaf(ndc2) = 1;
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount = ndcount+2;
    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0) break;
    }
  } else {
    end = 1;
  }
  return end;
}

arma::ivec Forest::find_split_ICON(size_t nd,
                                   const arma::umat& mat1Z,
                                   const arma::mat& mat1f,
                                   const arma::field<arma::umat>& mat2Zf, // dat
                                   const arma::uvec& isLeaf,
                                   const  arma::ucube& ranges,
                                   const arma::field<arma::uvec>& nodeSampleY,
                                   const arma::field<arma::uvec>& nodeSample,
                                   arma::mat& fmat,
                                   arma::umat& Smat,
                                   int ndcount,
                                   const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::mat fmatTerm = fmat.cols(arma::find(isLeaf == 1));
  arma::umat SmatTerm = Smat.cols(arma::find(isLeaf == 1));
  arma::uvec spSet = arma::shuffle( arma::regspace<arma::uvec>(0,P-1) );
  for(auto p : spSet.head(mtry)) {
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0),
                                                     ranges(nd, p, 1));
    int nel = 0;
    // int nelr = arma::sum(e( indY ) );
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        size_t z = mat1Z(indYj, p);
        if(z == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }
      if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum))/2;
        for(size_t i = 0; i < fmatTerm.n_cols; i++) {
          if(i != nd ) {
            for(size_t k = 0; k < K; k++) {
              if( (fLSum(k)+fRSum(k))*SmatTerm(k,i) <=  fmatTerm(k,i)*(SLSum(k)+SRSum(k)) ) {
                if( fmatTerm(k,i)*SLSum(k) <  fLSum(k)*SmatTerm(k,i) ) 
                {// lambdaR < lambdand <= lambdai < lambdaL
                  dICONTemp = dICONTemp +  fLSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SLSum(k) ;
                }
                if( fmatTerm(k,i)*SRSum(k) <  fRSum(k)*SmatTerm(k,i))
                {// lambdaL < lambdand <= lambdai < lambdaR
                  dICONTemp = dICONTemp + fRSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SRSum(k) ;
                }
              }else{
                if( fmatTerm(k,i)*SLSum(k) >  fLSum(k)*SmatTerm(k,i) )
                {// lambdaL < lambdai < lambdand < lambdaR
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SLSum(k) - fLSum(k)*SmatTerm(k,i) ;
                }
                if( fmatTerm(k,i)*SRSum(k) >  fRSum(k)*SmatTerm(k,i) )
                {// lambdaR < lambdai < lambdand < lambdaL
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SRSum(k) - fRSum(k)*SmatTerm(k,i) ;
                }
              }
            }
          }
        }
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
        fmat.col(ndcount + 1) = fLSum;
        fmat.col(ndcount + 2) = fRSum;
        Smat.col(ndcount + 1) = SLSum;
        Smat.col(ndcount + 2) = SRSum;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}
#include "ForestPrediction.h"
#include "globals.h"
#include <algorithm>


void ForestPrediction::transformZ(const arma::mat& z,
                                  arma::umat& z2,
                                  const arma::mat& matX,
                                  const arma::uvec& e,
                                  const arma::vec& breaks,
                                  const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = matX.col(p);
      // BE CAREFUL
      int j = 0;
      for(int i = 0; i < n; i++) {
        if(e(i) == 1) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          double z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) ) - zpref.begin())/(n-i+0.0);
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1 ;
          j++;
        }
        ind = ind + 1;
        if(ind.n_elem > 0) ind.shed_row(0);
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

void ForestPrediction::transformZH(const arma::mat& z,  // z on tg
                                   const arma::vec& tg, //grid of time point
                                   arma::umat& z2,
                                   const arma::mat& dat, // dat is predictor matrix
                                   const arma::vec& y,
                                   const arma::uvec& e,
                                   const arma::vec& breaks,
                                   const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  int NG = tg.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = dat.col(p);
      // BE CAREFUL
      int j = 0;
      double z2ecdf;
      for(int i = 0; i < n-1; i++) {
        if( y(i) <= tg(j) && tg(j) < y(i+1) ) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,j) )- zpref.begin())/(n-i+0.0)  ;
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1;
          j++;
          if(j == NG) break;
        }
        ind = ind + 1;
        if(ind.n_elem > 0) ind.shed_row(0);
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

// The new version
arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
                                        const arma::vec& y,
                                        const arma::uvec& e,
                                        const arma::field<arma::uvec>&& nodeSizeB0,
                                        const arma::umat&& nodeLabelB0,
                                        const arma::field<arma::uvec>&& tnd3B0,
                                        const arma::umat&& ids,
                                        const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  arma::vec w0 = arma::zeros<arma::vec>(NE);
  arma::vec w1 = arma::zeros<arma::vec>(NE);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = (tnd3B0(b));
    arma::uvec zY;
    // take the ids with uncensored survival times
    // IDs have been sorted when sampling is done; need to sort if not sorted before
    arma::uvec idb = ids.col(b);
    //arma::uvec idb = arma::sort(ids.col(b));
    arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
    int j = 0;
    int nidbe = idbe.n_elem;
    for(size_t i = 0; i < NE; i++) {
      int count = 0;
      while(j < nidbe) {
        if( idbe(j) == idY2(i) ) {
          count++;
          j++;
        } else {
          break;
        }
      }
      zY = zt2.col(i);
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
	varsp = vars(k);
	cutsp = values(k);
	if(zY(varsp) > cutsp) {
	  k = rcs(k);
	} else {
            k = lcs(k);
	}
	isl = il(k);
      }
      arma::uvec nbi = nodeSizeB0(b + NUM_TREE * i);
      if(nodeLabelB0(b,i) == k) {
	w0(i) += count ;
      }
      w1(i) += nbi(tnd3(k));
    }   
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w(e1) = w0/w1;
  w.replace(arma::datum::nan, 0);
  return exp(-cumsum(w));
}



arma::vec ForestPrediction::getHazard(const arma::umat& ztvec,
                                      const arma::vec& tg,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::mat& Kmat,
                                      const double h,
                                      const arma::field<arma::uvec>&& nodeSizeB0,
                                      const arma::umat&& nodeLabelB0,
                                      const arma::field<arma::uvec>&& tnd3B0,
                                      const arma::umat&& ids,
                                      const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  uint NG = tg.n_elem;
  arma::mat v0 = arma::zeros<arma::mat>(NE, NG);
  arma::mat v1 = arma::zeros<arma::mat>(NE, NG);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = tnd3B0(b);
    // take the ids with uncensored survival times
    arma::uvec idb = ids.col(b); // IDs have been sorted when sampling is done; need to sort if not sorted before
    //arma::uvec idb = arma::sort(ids.col(b));
    arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
    int nidbe = idbe.n_elem;
    for(size_t l = 0; l < NG; l++) {
      arma::uvec zl = ztvec.col(l);
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zl(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      int j = 0;
      for(size_t i = 0; i < NE; i++) {
        int count = 0;
        while(j < nidbe) {
          if( idbe(j) == idY2(i) ) {
            count++;
            j++;
          } else {
            break;
          }
        }
	if( y2(i) >= tg(l) - h && y2(i) <= tg(l) + h && nodeLabelB0(b,i) == k) {
          v0(i,l) += count;
        }
        arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
        v1(i,l) += nbi(tnd3(k));
      }
    }
  }
  arma::vec hz = arma::zeros<arma::vec>(NG);
  v1(arma::find(v1 == 0)).ones();
  for(size_t l = 0; l < NG; l++) {
    hz(l) = arma::sum(Kmat.col(l) % v0.col(l) / v1.col(l)) ;
  }
  return hz;
}


// The old version
arma::vec ForestPrediction::getSurvival2(const arma::umat& zt2,
                                        const arma::vec& y,
                                        const arma::uvec& e,
                                        const arma::field<arma::uvec>&& nodeSizeB0,
                                        const arma::umat&& nodeLabelB0,
                                        const arma::field<arma::uvec>&& tnd3B0,
                                        const arma::umat&& ids,
                                        const arma::field<arma::umat>&& trees)
{
  arma::uvec e1 = arma::find(e == 1);
  arma::vec y2 = y( e1 );
  arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
  arma::uvec idY2 = idAll( e1 );
  uint NE = y2.n_elem;
  uint NUM_TREE = trees.size();
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(size_t b = 0; b != NUM_TREE; b++) {
    // take the b th tree
    arma::umat tt = (trees(b));
    arma::uvec vars =  tt.col(0);
    arma::uvec values = tt.col(1);
    arma::uvec lcs = tt.col(2);
    arma::uvec rcs = tt.col(3);
    arma::uvec il = tt.col(4);
    arma::uvec tnd3 = (tnd3B0(b));
    arma::uvec zY;
    for(size_t i = 0; i < NE; i++) {
      // need to add this back if using sample-spliting
      //if(std::find(ids.begin_col(b), ids.end_col(b), idY2(i) ) != ids.end_col(b))
      zY = zt2.col(i);
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
	varsp = vars(k);
	cutsp = values(k);
	if(zY(varsp) > cutsp) {
	  k = rcs(k);
	} else {
	  k = lcs(k);
	}
	isl = il(k);
      }
      arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
      double den = nbi(tnd3(k));
      if(nodeLabelB0(b,i) == k) {
	w2(i) += (1.0/NUM_TREE)/den;
      }
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( e1 ) = w2;
  //return w;
  return exp(-cumsum(w));
}

ForestPrediction::ForestPrediction(const arma::umat& zy,
                                   const arma::field<arma::umat>& zt,
                                   const arma::umat& ids,
                                   const std::vector<std::shared_ptr<Tree> >& trees,
                                   arma::uword n)

{
  size_t nT = zt.size();
  uint NUM_TREE = trees.size();
  arma::umat ndy(NUM_TREE, nT);
  arma::field<arma::uvec> ndsz(NUM_TREE, nT);
  arma::field<arma::uvec> tnd3B(NUM_TREE);
  std::vector<std::shared_ptr<Tree> >::const_iterator it;
  int i = 0;
  for(it = trees.begin(); it != trees.end(); it++, i++) {
    // take the i th tree
    std::shared_ptr<Tree> tt = *it;
    arma::uvec vars =  tt->get_split_vars();
    arma::uvec lcs = tt->get_left_childs();
    arma::uvec rcs = tt->get_right_childs();
    arma::uvec values = tt->get_split_values();
    arma::uvec il = tt->get_isLeaf();
    int nNd = arma::accu(il);
    //arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
    arma::uvec tnd2 = (arma::find(il == 1));
    arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
    tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
    for(size_t j = 0; j != nT; j++) {
      arma::uvec zyj = zy.col( j );
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      size_t k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zyj(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndy(i,j) = k;
    }
    arma::uvec idi = ids.col(i);
    for(arma::uword c=0; c < nT; c++) {
      arma::uvec ndszic = arma::zeros<arma::uvec>(nNd);
      arma::umat m = zt(c);
      //arma::umat m(zt(c).memptr(), zt(c).n_rows, zt(c).n_cols,false);
      //arma::ivec idc = arma::conv_to<arma::ivec>::from(idi + m.n_cols - n);
      //arma::ivec idcp = idc( find(idc>=0) );
      arma::uvec idcp = idi( arma::find(idi >= n-m.n_cols) ) + m.n_cols - n;
      int j_end = idcp.n_elem;
      if(j_end > 0) {
        for(int j = 0; j != j_end; j++) {
          arma::uvec zti = m.col( idcp(j) );
          int isl = 0;
          int varsp = 0;
          size_t cutsp = 0;
          int k = 0;
          while(isl == 0) {
            varsp = vars(k);
            cutsp = values(k);
            if(zti(varsp) > cutsp) {
              k = rcs(k);
            } else {
              k = lcs(k);
            }
            isl = il(k);
          }
          //Rcpp::Rcout << k;
          ndszic(tnd3(k))++;
        }
      }
      //Rcpp::Rcout << ndszic;
      ndsz(i,c) = ndszic;
    }
    tnd3B(i) = tnd3;
  }
  this->nodeLabelB = ndy;
  this->nodeSizeB = ndsz;
  this->tnd3B = tnd3B;
}

// arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
//                                         const Data2* dat2,
//                                         const arma::umat& ids,
//                                         const std::vector<std::shared_ptr<Tree>>& trees)
// {
//   arma::vec y = dat2->get_Y();
//   arma::uvec e = dat2->get_e();
//   arma::vec y2 = y( arma::find(e == 1));
//   arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
//   arma::uvec idY2 = idAll( arma::find(e == 1));
//   uint NE = y2.n_elem;
//   uint NUM_TREE = trees.size();
//   arma::vec w2 = arma::zeros<arma::vec>(NE);
//
//
//   std::vector<std::shared_ptr<Tree>>::const_iterator it;
//   int b = 0;
//   for(it = trees.begin(); it != trees.end(); it++)
//   {
//     // take the i th tree
//     std::shared_ptr<Tree> tt = *it;
//     const arma::uvec& vars =  tt->get_split_vars();
//     const arma::uvec& lcs = tt->get_left_childs();
//     const arma::uvec& rcs = tt->get_right_childs();
//     const arma::uvec& values = tt->get_split_values();
//     const arma::uvec& il = tt->get_isLeaf();
//
//     arma::uvec tnd3 = tnd3B(b);
//
//     for(int i = 0; i < NE; i++)
//     {
//       arma::uvec zY = (zt2.col(i));
//
//       int isl = 0;
//       int varsp = 0;
//       int cutsp = 0;
//       int k = 0;
//
//       while(isl == 0)
//       {
//         varsp = vars(k);
//         cutsp = values(k);
//         if(zY(varsp) > cutsp)
//         {
//           k = rcs(k);
//         }else{
//           k = lcs(k);
//         }
//         isl = il(k);
//       }
//
//       int nd = k;
//       arma::uvec nbi = nodeSizeB(b,i);
//       arma::uword den = nbi(tnd3(nd));
//       if(nodeLabelB(b,i) == nd && den > 0 && arma::sum(idY2(i) == ids.col(b))>0 )//&&
//       {
//         w2(i) = w2(i) + (1.0/NUM_TREE)/den;
//       }
//     }
//     b++;
//   }
//
//
//   arma::vec w = arma::zeros<arma::vec>(y.n_elem);
//   w( arma::find(e == 1)) = w2;
//
//   //return w;
//   return exp(-cumsum(w));
// }
//
//
//
// ForestPrediction::ForestPrediction(const Data2* dat2,
//                                    const arma::umat& ids,
//                                    const std::vector<std::shared_ptr<Tree>>& trees,
//                                    arma::uword n)
//
// {
//   arma::umat zy = dat2->get_zy();
//   arma::field<arma::umat> zt = dat2->get_zt();
//
//   int nT = zt.size();
//   uint NUM_TREE = trees.size();
//   arma::umat ndy(NUM_TREE, nT);
//   arma::field<arma::uvec> ndsz(NUM_TREE, nT);
//   arma::field<arma::uvec> tnd3B(NUM_TREE);
//
//   std::vector<std::shared_ptr<Tree>>::const_iterator it;
//   int i = 0;
//   for(it = trees.begin(); it != trees.end(); it++)
//   {
//     // take the i th tree
//     std::shared_ptr<Tree> tt = *it;
//     arma::uvec vars =  tt->get_split_vars();
//     arma::uvec lcs = tt->get_left_childs();
//     arma::uvec rcs = tt->get_right_childs();
//     arma::uvec values = tt->get_split_values();
//     arma::uvec il = tt->get_isLeaf();
//     int nNd = arma::accu(il);
//     arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
//     arma::uvec tnd2 = tnd(arma::find(il == 1));
//     arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
//     tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
//
//     for(int j = 0; j != nT; j++)
//     {
//       arma::uvec zyj = (zy.col( j  ));
//
//       int isl = 0;
//       int varsp = 0;
//       int cutsp = 0;
//       int k = 0;
//       while(isl == 0)
//       {
//         varsp = vars(k);
//         cutsp = values(k);
//         if(zyj(varsp) > cutsp)
//         {
//           k = rcs(k);
//         }else{
//           k = lcs(k);
//         }
//         isl = il(k);
//       }
//       ndy(i,j) = k;
//     }
//
//     arma::uvec idi = ids.col(i);
//
//     for(arma::uword c=0; c < nT; c++)
//     {
//       arma::uvec ndszic = arma::zeros<arma::uvec>(nNd);
//       arma::umat m = zt(c);
//       arma::ivec idc = arma::conv_to<arma::ivec>::from(idi + m.n_cols - n);
//       arma::ivec idcp = idc( find(idc>=0) );
//       int j_end = idcp.n_elem;
//
//
//       if(j_end > 0)
//       {
//         for(int j = 0; j != j_end; j++)
//         {
//           arma::uvec zti = (m.col( idcp(j) ));
//
//           int isl = 0;
//           int varsp = 0;
//           int cutsp = 0;
//           int k = 0;
//           while(isl == 0)
//           {
//             varsp = vars(k);
//             cutsp = values(k);
//             if(zti(varsp) > cutsp)
//             {
//               k = rcs(k);
//             }else{
//               k = lcs(k);
//             }
//             isl = il(k);
//           }
//
//           //Rcpp::Rcout << k;
//           ndszic(tnd3(k))++;
//         }
//       }
//       //Rcpp::Rcout << ndszic;
//       ndsz(i,c) = ndszic;
//     }
//
//
//     tnd3B(i) = tnd3;
//
//     i++;
//
//   }
//
//
//
//   this->nodeLabelB = ndy;
//   this->nodeSizeB = ndsz;
//   this->tnd3B = tnd3B;
//
// }

// Time-independent way
// void ForestPrediction::transformZ0(const arma::mat& z,
//                                    arma::umat& z2,
//                                    const arma::mat& matX,
//                                    const arma::uvec& e,
//                                    const arma::vec& breaks,
//                                    const arma::uvec& disc)
// {
//   int P = z.n_rows;
//   int n = e.n_elem;
//   arma::vec::const_iterator bb = breaks.begin();
//   arma::vec::const_iterator be = breaks.end();
//
//   for(int p = 0; p < P; p++)
//   {
//     if(disc(p) == 0)
//     {
//       arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
//       arma::vec zp = matX.col( p );
//       arma::vec zpref = arma::sort(zp( ind ));
//       // BE CAREFUL
//       int j = 0;
//       for(int i = 0; i < n; i++)
//       {
//         if(e(i) == 1)
//         {
//           double z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) )- zpref.begin())/(n+0.0)  ;
//           z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1 ;
//           j++;
//         }
//       }
//     }
//     else{
//       z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
//     }
//   }
// }

// Averaging CHF
// arma::vec ForestPrediction::getSurvival(const arma::umat& zt2,
//                                          const arma::vec& y,
//                                          const arma::uvec& e,
//                                          const arma::field<arma::uvec>&& nodeSizeB0,
//                                          const arma::umat&& nodeLabelB0,
//                                          const arma::field<arma::uvec>&& tnd3B0,
//                                          const arma::umat&& ids,
//                                          const arma::field<arma::umat>&& trees)
// {
//   arma::uvec e1 = arma::find(e == 1);
//   arma::vec y2 = y( e1 );
//   arma::uvec idAll = arma::regspace<arma::uvec>(0, y.n_elem-1);
//   arma::uvec idY2 = idAll( e1 );
//   uint NE = y2.n_elem;
//   uint NUM_TREE = trees.size();
//
//   arma::vec Hes = arma::zeros<arma::vec>(y.n_elem);
//
//   for(int b = 0; b != NUM_TREE; b++)
//   {
//     // take the b th tree
//     arma::umat tt = (trees(b));
//     arma::uvec vars =  tt.col(0);
//     arma::uvec values = tt.col(1);
//     arma::uvec lcs = tt.col(2);
//     arma::uvec rcs = tt.col(3);
//     arma::uvec il = tt.col(4);
//     arma::uvec tnd3 = (tnd3B0(b));
//
//     arma::uvec zY;
//
//     // take the ids with uncensored survival times
//     arma::uvec idb = ids.col(b); // IDs have been sorted when sampling is done; need to sort if not sorted before
//     //arma::uvec idb = arma::sort(ids.col(b));
//     arma::uvec idbe = idb( arma::find( e(idb) == 1 )  );
//
//     int j = 0;
//     int nidbe = idbe.n_elem;
//
//     arma::vec wb0 = arma::zeros<arma::vec>(NE);
//     arma::vec wb1 = arma::zeros<arma::vec>(NE);
//
//     for(int i = 0; i < NE; i++)
//     {
//       int count = 0;
//       while(j < nidbe)
//       {
//         if( idbe(j) == idY2(i) )
//         {
//           count++;
//           j++;
//         }else{
//           break;
//         }
//       }
//
//
//       {
//         zY = zt2.col(i);
//
//         int isl = 0;
//         int varsp = 0;
//         int cutsp = 0;
//         int k = 0;
//
//         while(isl == 0)
//         {
//           varsp = vars(k);
//           cutsp = values(k);
//           if(zY(varsp) > cutsp)
//           {
//             k = rcs(k);
//           }else{
//             k = lcs(k);
//           }
//           isl = il(k);
//         }
//
//         arma::uvec nbi = nodeSizeB0(b+NUM_TREE*i);
//
//         if(nodeLabelB0(b,i) == k)
//         {
//           wb0(i) += count ;
//         }
//         wb1(i) += nbi(tnd3(k));
//       }
//
//     }
//     arma::vec wb = arma::zeros<arma::vec>(y.n_elem);
//     wb(e1) = wb0/wb1;
//     wb.replace(arma::datum::nan, 0);
//     Hes += arma::cumsum(wb)/NUM_TREE;
//   }
//
//   return exp(-Hes);
// }
// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// make_mat2
arma::field<arma::mat> make_mat2(arma::vec& tk0, arma::vec& Y0, arma::uvec& id0, arma::mat& z0);
RcppExport SEXP _rocTree_make_mat2(SEXP tk0SEXP, SEXP Y0SEXP, SEXP id0SEXP, SEXP z0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type tk0(tk0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type id0(id0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z0(z0SEXP);
    rcpp_result_gen = Rcpp::wrap(make_mat2(tk0, Y0, id0, z0));
    return rcpp_result_gen;
END_RCPP
}
// make_mat2_t
arma::field<arma::mat> make_mat2_t(arma::vec& tk0, arma::vec& Y0, arma::uvec& id0, arma::mat& z0);
RcppExport SEXP _rocTree_make_mat2_t(SEXP tk0SEXP, SEXP Y0SEXP, SEXP id0SEXP, SEXP z0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type tk0(tk0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type id0(id0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z0(z0SEXP);
    rcpp_result_gen = Rcpp::wrap(make_mat2_t(tk0, Y0, id0, z0));
    return rcpp_result_gen;
END_RCPP
}
// rocForest_C
SEXP rocForest_C(const arma::mat& mat1f0, const arma::umat& mat1Z0, const arma::field<arma::umat>& mat2Zf0, const arma::umat& r0, const arma::field<arma::umat>& zt0, const arma::umat& zy0, const arma::uvec& e0, int spCriterion, int numTree, int minNode1, int minSplit1, int maxNode, int mtry);
RcppExport SEXP _rocTree_rocForest_C(SEXP mat1f0SEXP, SEXP mat1Z0SEXP, SEXP mat2Zf0SEXP, SEXP r0SEXP, SEXP zt0SEXP, SEXP zy0SEXP, SEXP e0SEXP, SEXP spCriterionSEXP, SEXP numTreeSEXP, SEXP minNode1SEXP, SEXP minSplit1SEXP, SEXP maxNodeSEXP, SEXP mtrySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat1f0(mat1f0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type mat1Z0(mat1Z0SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::umat>& >::type mat2Zf0(mat2Zf0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::umat>& >::type zt0(zt0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type zy0(zy0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< int >::type spCriterion(spCriterionSEXP);
    Rcpp::traits::input_parameter< int >::type numTree(numTreeSEXP);
    Rcpp::traits::input_parameter< int >::type minNode1(minNode1SEXP);
    Rcpp::traits::input_parameter< int >::type minSplit1(minSplit1SEXP);
    Rcpp::traits::input_parameter< int >::type maxNode(maxNodeSEXP);
    Rcpp::traits::input_parameter< int >::type mtry(mtrySEXP);
    rcpp_result_gen = Rcpp::wrap(rocForest_C(mat1f0, mat1Z0, mat2Zf0, r0, zt0, zy0, e0, spCriterion, numTree, minNode1, minSplit1, maxNode, mtry));
    return rcpp_result_gen;
END_RCPP
}
// predict_rocForest_C
SEXP predict_rocForest_C(const arma::mat& zraw0, const arma::vec& y0, const arma::uvec& e0, const Rcpp::List& forestobj, const arma::mat& matX, const arma::uvec& disc, const arma::vec& breaks);
RcppExport SEXP _rocTree_predict_rocForest_C(SEXP zraw0SEXP, SEXP y0SEXP, SEXP e0SEXP, SEXP forestobjSEXP, SEXP matXSEXP, SEXP discSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zraw0(zraw0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type forestobj(forestobjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type matX(matXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type disc(discSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_rocForest_C(zraw0, y0, e0, forestobj, matX, disc, breaks));
    return rcpp_result_gen;
END_RCPP
}
// rocTree_C
SEXP rocTree_C(const arma::mat& mat1f0, const arma::umat& mat1Z0, const arma::field<arma::umat>& mat2Zf0, const arma::umat& r0, const arma::field<arma::umat>& zt0, const arma::umat& zy0, const arma::uvec& e0, int spCriterion, int numFold, int minNode1, int minSplit1, int maxNode);
RcppExport SEXP _rocTree_rocTree_C(SEXP mat1f0SEXP, SEXP mat1Z0SEXP, SEXP mat2Zf0SEXP, SEXP r0SEXP, SEXP zt0SEXP, SEXP zy0SEXP, SEXP e0SEXP, SEXP spCriterionSEXP, SEXP numFoldSEXP, SEXP minNode1SEXP, SEXP minSplit1SEXP, SEXP maxNodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat1f0(mat1f0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type mat1Z0(mat1Z0SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::umat>& >::type mat2Zf0(mat2Zf0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::umat>& >::type zt0(zt0SEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type zy0(zy0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< int >::type spCriterion(spCriterionSEXP);
    Rcpp::traits::input_parameter< int >::type numFold(numFoldSEXP);
    Rcpp::traits::input_parameter< int >::type minNode1(minNode1SEXP);
    Rcpp::traits::input_parameter< int >::type minSplit1(minSplit1SEXP);
    Rcpp::traits::input_parameter< int >::type maxNode(maxNodeSEXP);
    rcpp_result_gen = Rcpp::wrap(rocTree_C(mat1f0, mat1Z0, mat2Zf0, r0, zt0, zy0, e0, spCriterion, numFold, minNode1, minSplit1, maxNode));
    return rcpp_result_gen;
END_RCPP
}
// predict_rocTree_C
arma::vec predict_rocTree_C(const arma::mat& zraw0, const arma::vec& y0, const arma::uvec& e0, const Rcpp::List& treeobj, const arma::mat& matX, const arma::uvec& disc, const arma::vec& breaks);
RcppExport SEXP _rocTree_predict_rocTree_C(SEXP zraw0SEXP, SEXP y0SEXP, SEXP e0SEXP, SEXP treeobjSEXP, SEXP matXSEXP, SEXP discSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zraw0(zraw0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type treeobj(treeobjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type matX(matXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type disc(discSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_rocTree_C(zraw0, y0, e0, treeobj, matX, disc, breaks));
    return rcpp_result_gen;
END_RCPP
}
// predict_rocTreeHZ_C
arma::vec predict_rocTreeHZ_C(const arma::mat& zraw0, const arma::vec& tg0, const arma::vec& y0, const arma::uvec& e0, const arma::mat& fy20, const double h0, const Rcpp::List& treeobj, const arma::mat& matX, const arma::uvec& disc, const arma::vec& breaks);
RcppExport SEXP _rocTree_predict_rocTreeHZ_C(SEXP zraw0SEXP, SEXP tg0SEXP, SEXP y0SEXP, SEXP e0SEXP, SEXP fy20SEXP, SEXP h0SEXP, SEXP treeobjSEXP, SEXP matXSEXP, SEXP discSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zraw0(zraw0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tg0(tg0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fy20(fy20SEXP);
    Rcpp::traits::input_parameter< const double >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type treeobj(treeobjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type matX(matXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type disc(discSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_rocTreeHZ_C(zraw0, tg0, y0, e0, fy20, h0, treeobj, matX, disc, breaks));
    return rcpp_result_gen;
END_RCPP
}
// predict_rocForestHZ_C
arma::vec predict_rocForestHZ_C(const arma::mat& zraw0, const arma::vec& tg0, const arma::vec& y0, const arma::uvec& e0, const arma::mat& Kmat0, const double h0, const Rcpp::List& forestobj, const arma::mat& matX, const arma::uvec& disc, const arma::vec& breaks);
RcppExport SEXP _rocTree_predict_rocForestHZ_C(SEXP zraw0SEXP, SEXP tg0SEXP, SEXP y0SEXP, SEXP e0SEXP, SEXP Kmat0SEXP, SEXP h0SEXP, SEXP forestobjSEXP, SEXP matXSEXP, SEXP discSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type zraw0(zraw0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tg0(tg0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Kmat0(Kmat0SEXP);
    Rcpp::traits::input_parameter< const double >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type forestobj(forestobjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type matX(matXSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type disc(discSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_rocForestHZ_C(zraw0, tg0, y0, e0, Kmat0, h0, forestobj, matX, disc, breaks));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rocTree_make_mat2", (DL_FUNC) &_rocTree_make_mat2, 4},
    {"_rocTree_make_mat2_t", (DL_FUNC) &_rocTree_make_mat2_t, 4},
    {"_rocTree_rocForest_C", (DL_FUNC) &_rocTree_rocForest_C, 13},
    {"_rocTree_predict_rocForest_C", (DL_FUNC) &_rocTree_predict_rocForest_C, 7},
    {"_rocTree_rocTree_C", (DL_FUNC) &_rocTree_rocTree_C, 12},
    {"_rocTree_predict_rocTree_C", (DL_FUNC) &_rocTree_predict_rocTree_C, 7},
    {"_rocTree_predict_rocTreeHZ_C", (DL_FUNC) &_rocTree_predict_rocTreeHZ_C, 10},
    {"_rocTree_predict_rocForestHZ_C", (DL_FUNC) &_rocTree_predict_rocForestHZ_C, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_rocTree(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "Data.h"
#include "Data2.h"
#include "globals.h"
#include "Tree.h"
#include "Forest.h"
#include "ForestPrediction.h"
#include "TreeGrow.h"
#include "TreePrediction.h"

#ifndef OLD_WIN_R_BUILD
#include <thread>
#include <chrono>
#endif

//' @noRd
// [[Rcpp::export]]
SEXP rocForest_C(const arma::mat& mat1f0,
		 const arma::umat& mat1Z0,
		 const arma::field<arma::umat>& mat2Zf0,
		 const arma::umat& r0,
		 const arma::field<arma::umat>& zt0,
		 const arma::umat& zy0,
		 const arma::uvec& e0,
		 int spCriterion, // dICON if 1, ICON otherwise
		 int numTree,
		 int minNode1,
		 int minSplit1,
		 int maxNode,
		 int mtry) {
  int K = mat1f0.n_rows;
  int n = mat1Z0.n_rows;
  std::vector<std::shared_ptr<Tree> > trees;
  trees.reserve(numTree);
  Forest ff(spCriterion, numTree, maxNode, minNode1, minSplit1, K, mtry);
  //arma::umat ids(n, numTree);
  arma::umat ids0(n, numTree);
  //ff.sampleWithoutReplacementSplit(n, n, ids);
  ff.sampleWithReplacementSplit(n, n, ids0);
  // Bootstrapping
  arma::umat id1 = ids0;//ids.rows( arma::regspace<arma::uvec>(0, s-1)  );
  ff.trainRF(trees, mat1Z0, mat1f0, mat2Zf0, r0, id1, e0);
  arma::field<arma::umat> treeList(numTree);
  std::vector<std::shared_ptr<Tree> >::const_iterator it;
  int i = 0;
  for(it = trees.begin(); it != trees.end(); it++, i++) {
    std::shared_ptr<Tree> tt = *it;
    const arma::uvec& vars0 = tt->get_split_vars();
    arma::umat treeMat(vars0.n_elem,5);
    treeMat.col(0) = vars0;
    treeMat.col(1) = tt->get_split_values();
    treeMat.col(2) = tt->get_left_childs();
    treeMat.col(3) = tt->get_right_childs();
    treeMat.col(4) = tt->get_isLeaf();
    treeList(i) = treeMat;
  }
  // use bootstrap observations
  arma::umat id2 = ids0;
  //ids.rows( arma::regspace<arma::uvec>(0, n-1)  );
  ForestPrediction fp(zy0, zt0, id2, trees, n);
  return Rcpp::List::create(Rcpp::Named("trees") = treeList,
                            Rcpp::Named("nodeLabel") = fp.get_nodeLabel(),
                            Rcpp::Named("nodeSize") = fp.get_nodeSize(),
                            Rcpp::Named("nodeMap") = fp.get_nodeMap(),
                            Rcpp::Named("boot.id") = id2);
}

// [[Rcpp::export]]
SEXP predict_rocForest_C(const arma::mat& zraw0,
			 const arma::vec& y0,
			 const arma::uvec& e0,
			 const Rcpp::List& forestobj,
			 const arma::mat& matX,
			 const arma::uvec& disc,
			 const arma::vec& breaks) {
  arma::umat z0(zraw0.n_rows, arma::sum(e0));
  ForestPrediction::transformZ(zraw0, z0, matX, e0, breaks, disc);
  arma::vec sy = ForestPrediction::getSurvival(z0,
                                               y0,
                                               e0,
                                               forestobj[2],
                                               forestobj[1],
                                               forestobj[3],
                                               forestobj[4],
                                               forestobj[0]);
  return Rcpp::wrap(sy);
}

//' Main function for tree.
//' Data are prepared outside of the function.
//' 
//' @param mat1f0 is a K x n matrix; K is the number of cutpoint
//'               for i = 1, ..., n, k = 1, ..., K,
//'               the transpose of it has the (i, j)th element of f_i(t_j) 
//' @param mat1Z0 is a n x p matrix; p is the number of covariates
//'               instead of covaraites values, this gives the 'order' that indicates which
//'               time intervals the covaraite values lie in.
//' @param mat2Z0 is a list of length K. The lth list consists of a r x p matrix that
//'               specified the order of covariates for subjects at risk at tk;
//'               r is the number of subjects at risk at t_l.
//' @param r0     is a 2 x p matrix. Each column specified the smallest and the largest intervals
//'               the pth covariates lie in.
//' @param zt0    is a list of length equal to the numbers of events (\sum\Delta).
//'               This is created similar to `mat2Z` but it stores the order of covariates for
//'               subjects at risk at the event times (Y0)
//' @param zy0    is a transpose of a subset of `mat1Z0` on these event times,
//'               e.g., .zy <- t(.mat1Z[.E0 == 1,])
//' @param e0     censoring indicator; 1 = event and 0 = control
//'
//'
//' Tunning parameters
//' @param spCriterion specifies the splitting criterion; 1 = DICON and 0 = ICON
//' @param numFold is the number of fold used in cross-validation
//' @param minNode1 is the minimum number of baseline observations in each terminal node; default = 15
//'                 we used time-invariant tree, so the node size criterion is on the baseline obs.
//' @param minSplit1 is the minimum number of baseline observations in each splitable node; default 30
//' @param maxNode is the maximum number of terminal nodes
//'
//' Returns a list with the following elements
//'   *treeMat - a xx by 5 matrix that consists of 5 columns;
//'     the number of rows depends on how depth the tree is
//'      Column1: Which covariate to split? Index corresponds to the columns in .X
//'      Column2: What value to split? (this is the ranked value; corresponding to which interval.
//'      Column3: Left child node (<= the cut value)
//'      Column4: Right child node (> the cut value)
//'      Column5: Node indicator; 1 = terminal node, 0 = non-terminal node
//'
//'   *nodeLabel - a n-dimensional vector, where n is the number of unique uncensored survival times.
//'
//'   *nodeSize - a xx by n matrix;
//'     the number of rows = number of terminal nodes
//'     the number of columns = number of unique ordered uncensored survival times  (Y0)
//'       This gives node size at each of the terminal nodes at each unique Y0
//' 
//'   *nodeMap - a xx-dimensional vector, where xx is the number of rows in treeMat.
//'     The elements indicates the element in 'nodeSize' to be used as risk set size in survival pred.
//'     Only the elements correspond to terminal nodes (treeMat[,5] == 1) will be used.
//' 
//' @noRd
// [[Rcpp::export]]
SEXP rocTree_C(const arma::mat& mat1f0,
	       const arma::umat& mat1Z0,
	       const arma::field<arma::umat>& mat2Zf0,
	       const arma::umat& r0,
	       const arma::field<arma::umat>& zt0,
	       const arma::umat& zy0,
	       const arma::uvec& e0,
	       int spCriterion, 
	       int numFold,
	       int minNode1,
	       int minSplit1,
	       int maxNode) {
  int K = mat1f0.n_rows;
  TreeGrow tg(numFold, K, spCriterion, maxNode, minNode1, minSplit1);
  std::shared_ptr<Tree> tr2 = tg.trainCV(mat1Z0, mat1f0, mat2Zf0, r0, e0);
  const arma::uvec& vars0 = tr2->get_split_vars();
  arma::umat treeMat(vars0.n_elem,5);
  treeMat.col(0) = vars0;
  treeMat.col(1) = tr2->get_split_values();
  treeMat.col(2) = tr2->get_left_childs();
  treeMat.col(3) = tr2->get_right_childs();
  treeMat.col(4) = tr2->get_isLeaf();
  TreePrediction tp( zy0, zt0, treeMat.col(0), treeMat.col(1),
		     treeMat.col(2), treeMat.col(3), treeMat.col(4));
  return Rcpp::List::create(Rcpp::Named("treeMat") = treeMat,
                            Rcpp::Named("nodeLabel") = tp.get_nodeLabel(),
                            Rcpp::Named("nodeSize") = tp.get_nodeSize(),
                            Rcpp::Named("nodeMap") = tp.get_nodeMap());
}


//' Main function for survival prediction from a rocTree object
//'
//' @param zraw0 is a matrix consists of new covariates.
//'             The number of columns equal to the total number of unique survival times (n).
//'             The number of rows equal to the number of covariates
//' @param y0 is a n-dimensional vector consists of all survival times 
//' @param e0 is the censoring indicator
//' @param treeobj is the rocTree object
//' @param matX is the original covariate matrix,
//'             the number of rows is equal to the total observed survival times
//' @param disc is a vector with length equal to the number of covariates.
//'             This indicates whether a covariate is continuous (0) or not (1)
//' @param breaks is the cutoff points, e.g., cutoff <- (1:nc) / (nc+1)
//'             This gives the boundaries of intervals
//'
//' @noRd
// [[Rcpp::export]]
arma::vec predict_rocTree_C(const arma::mat& zraw0,
			    const arma::vec& y0,
			    const arma::uvec& e0,
			    const Rcpp::List& treeobj,
			    const arma::mat& matX, 
			    const arma::uvec& disc,
			    const arma::vec& breaks) {
  arma::umat z0(zraw0.n_rows, arma::sum(e0));
  ForestPrediction::transformZ(zraw0, z0, matX, e0, breaks, disc);
  arma::umat nodeSize2 = treeobj[2];
  arma::uvec nodeLabel2 = treeobj[1];
  arma::uvec tnd32 = treeobj[3];
  arma::umat treeMat2 = treeobj[0];
  arma::vec sy = TreePrediction::getSurvival(z0, y0, e0, nodeSize2, nodeLabel2, tnd32, treeMat2);
  return sy;
}

//' Main function for hazard prediction from a rocTree object
//'
//' Because of smoothing, the hazard estimate is based on fewer time points
//'
//' @param zraw0 is a matrix consists of new covariates.
//'              The number of columns equal to the length of `tg0`
//'              The number of rows equal to the number of covariates
//' @param tg0 is a vector that gives the time points where the prediction takes place.
//'              This corresponds to zraw0
//' @param y0 is a n-dimensional vector consists of all survival times
//' @param e0 is the censoring indicator
//' @param fy20 is a new version of `.mat1f` based on `zraw0` and `tg0`
//' @param h0 is the smoothing parameter
//' @param treeobj is the rocTree object
//' @param matX is the original covariate matrix,
//'             the number of rows is equal to the total observed survival times
//' @param disc is a vector with length equal to the number of covariates.
//'             This indicates whether a covariate is continuous (0) or not (1)
//' @param breaks is the cutoff points, e.g., cutoff <- (1:nc) / (nc+1)
//'             This gives the boundaries of intervals
//'
//' @noRd
// [[Rcpp::export]]
arma::vec predict_rocTreeHZ_C(const arma::mat& zraw0,
			      const arma::vec& tg0,
			      const arma::vec& y0,
			      const arma::uvec& e0,
			      const arma::mat& fy20,
			      const double h0,
			      const Rcpp::List& treeobj,
			      const arma::mat& matX, // the columns of covariates in dat from simu()
			      const arma::uvec& disc,
			      const arma::vec& breaks) {
  arma::umat z0(zraw0.n_rows, zraw0.n_cols);
  TreePrediction::transformZH(zraw0, tg0, z0, matX, y0, e0, breaks, disc);
  arma::umat nodeSize2 = treeobj[2];
  arma::uvec nodeLabel2 = treeobj[1];
  arma::uvec tnd32 = treeobj[3];
  arma::umat treeMat2 = treeobj[0];
  arma::vec hztg = TreePrediction::getHazard(z0, tg0,
                                             y0, e0,
                                             fy20, h0,
                                             nodeSize2, nodeLabel2, tnd32, treeMat2);
  return hztg;
}

//' @noRd
// [[Rcpp::export]]
arma::vec predict_rocForestHZ_C(const arma::mat& zraw0,
				const arma::vec& tg0,
				const arma::vec& y0,
				const arma::uvec& e0,
				const arma::mat& Kmat0,
				const double h0,
				const Rcpp::List& forestobj,
				const arma::mat& matX, // the columns of covariates in dat from simu()
				const arma::uvec& disc,
				const arma::vec& breaks) {
  arma::umat z0(zraw0.n_rows, zraw0.n_cols);
  ForestPrediction::transformZH(zraw0, tg0, z0, matX, y0, e0, breaks, disc);
  arma::vec hztg = ForestPrediction::getHazard(z0, tg0,
                                               y0, e0,
                                               Kmat0, h0,
                                               forestobj[2],
                                               forestobj[1],
                                               forestobj[3],
                                               forestobj[4],
                                               forestobj[0]);
  return hztg;
}


#include "Tree.h"
#include "globals.h"

const arma::uvec& Tree::get_split_vars() const
{
  return split_vars;
}

const arma::uvec& Tree::get_split_values() const
{
  return split_values;
}

const arma::uvec& Tree::get_left_childs() const
{
  return left_childs;
}

const arma::uvec& Tree::get_right_childs() const
{
  return right_childs;
}

const arma::uvec& Tree::get_isLeaf() const
{
  return isLeaf;
}

const arma::uvec& Tree::get_parents() const
{
  return parents;
}

void Tree::setzero(uint i, uint ndcount) {
  uint lid = left_childs(i);
  uint rid = right_childs(i);
  if(lid <= ndcount && isLeaf(lid) == 0) {
    setzero(lid,ndcount);
  }
  if(rid <= ndcount && isLeaf(rid) == 0) {
    setzero(rid,ndcount);
  }
  left_childs(i) = 0;
  right_childs(i) = 0;
  split_values(i) = 0;
  split_vars(i) = 0;
}

void Tree::cut(arma::uvec& nodeTerminal)
{
  arma::uvec isLeaf2 = isLeaf;
  isLeaf2.zeros();
  isLeaf2( nodeTerminal ).ones();
  int ndcount = isLeaf.n_elem-1;
  int ndcount2 = nodeTerminal(nodeTerminal.n_elem-1);
  for(int i = 0; i <= ndcount2; i++) {
    if(isLeaf2(i) == 1 && isLeaf(i) == 0) {
      setzero(i, ndcount);
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount2);
  right_childs = right_childs.elem(nonEmpty);
  left_childs = left_childs.elem(nonEmpty);
  split_values = split_values.elem(nonEmpty);
  split_vars = split_vars(nonEmpty);
  isLeaf = isLeaf2(nonEmpty);
}

double Tree::get_ICONTrain(const arma::uvec& isLeafTemp,
                           const arma::mat& fmat,
                           const arma::umat& Smat)
{
  arma::mat fmatTemp = fmat.cols( arma::find(isLeafTemp == 1) );
  arma::umat SmatTemp = Smat.cols( arma::find(isLeafTemp == 1) );
  int numLeafTemp = arma::sum(isLeafTemp == 1);
  uint K = fmat.n_rows;
  //SmatTemp.print("a");
  arma::vec icon = arma::zeros<arma::vec>(K);
  for(size_t k = 0; k != K; k++) {
    for(int i = 0; i < numLeafTemp; i++) {
      for(int j = 0; j < i; j++) {
        double fi = fmatTemp(k,i);
        double fj = fmatTemp(k,j);
        int Sj = SmatTemp(k,j);
        int Si = SmatTemp(k,i);
        if(fi*Sj > fj*Si) {
          icon(k) = icon(k) + fi*Sj;
        } else {
          icon(k) = icon(k) + fj*Si;
        }
      }
      icon(k) = icon(k) + fmatTemp(k,i)*SmatTemp(k,i)*0.5;
    }
  }
  icon = icon / (arma::sum(fmatTemp, 1) % arma::sum(SmatTemp, 1));
  return arma::sum(icon)/K;
}

//
// void Tree::findOptimalSizekSubtree(arma::mat& fmat, arma::umat& Smat,
//                                    arma::vec& iconAll, arma::field<arma::uvec>& nodeSetList, uint numLeaf)
// {
//   arma::uvec nodeID = arma::regspace<arma::uvec>(0, isLeaf.n_elem-1);
//   arma::uvec isLeafTemp = isLeaf;
//   arma::uvec isLeafTemp2 = isLeafTemp;
//   arma::uvec isLeafTemp3= isLeafTemp;
//   iconAll(numLeaf-1) = get_ICONTrain(isLeaf, fmat, Smat);
//   nodeSetList(numLeaf-1) = nodeID(arma::find(isLeafTemp == 1) );
//   while(numLeaf >= 2)
//   {
//     int cutNd;
//     double iconMax = 0;
//     double iconl = 0;
//     arma::uvec nodeTermTemp = nodeID(arma::find(isLeafTemp == 1) );
//
//     isLeafTemp2 = isLeafTemp;
//     isLeafTemp3 = isLeafTemp;
//
//     if(nodeTermTemp.n_elem > 1)
//     {
//       for(int l = 0; l < nodeTermTemp.n_elem; l++)
//       {
//         isLeafTemp2 = isLeafTemp3;
//         cutNd = nodeTermTemp(l);
//
//         uint pid = parents(cutNd);
//         if( arma::sum( nodeTermTemp == right_childs(pid) )>0 && arma::sum( nodeTermTemp == left_childs(pid))>0 )
//         {
//           isLeafTemp2( pid ) = 1;
//           isLeafTemp2( right_childs(pid) ) = 0;
//           isLeafTemp2( left_childs(pid) ) = 0;
//
//           iconl = get_ICONTrain(isLeafTemp2, fmat, Smat);
//           if(iconl >= iconMax) // = is used for the case with only one terminal nodes
//           {
//             iconMax = iconl;
//             nodeSetList(numLeaf-2) = nodeID(arma::find(isLeafTemp2 == 1) );
//             isLeafTemp = isLeafTemp2;
//           }
//         }
//
//       }
//     }else{
//     }
//
//     iconAll(numLeaf-2) = iconMax;
//     numLeaf--;
//   }
//
// }



void Tree::findOptimalSizekSubtree(arma::mat& fmat, arma::umat& Smat,
                                   arma::vec& iconAll, arma::field<arma::uvec>& nodeSetList,
				   uint numLeaf)
{
  arma::uvec nodeID = arma::regspace<arma::uvec>(0, isLeaf.n_elem-1);
  arma::uvec isLeafTemp = arma::zeros<arma::uvec>(isLeaf.n_elem);
  isLeafTemp(0) = 1;
  iconAll(0) = get_ICONTrain(isLeafTemp, fmat, Smat);
  nodeSetList(0) = nodeID(arma::find(isLeafTemp == 1) );
  arma::uvec isLeafTemp2 = isLeafTemp;
  arma::uvec isLeafTemp3 = isLeafTemp;
  size_t i = 1;
  while(i < numLeaf - 1) {
    double iconMax = 0;
    double iconl = 0;
    isLeafTemp2 = isLeafTemp;
    isLeafTemp3 = isLeafTemp;
    arma::uvec nodeTermTemp = nodeID(arma::find(isLeafTemp == 1) );
    for(size_t l = 0; l < nodeTermTemp.n_elem; l++ ) {
      int spNd = nodeTermTemp(l);
      if(isLeaf(spNd) == 0) {
        isLeafTemp2 = isLeafTemp3;
        int lid = left_childs(spNd);
        int rid = right_childs(spNd);
        isLeafTemp2(spNd) = 0;
        isLeafTemp2(lid) = 1;
        isLeafTemp2(rid) = 1;
        iconl = get_ICONTrain(isLeafTemp2, fmat, Smat);
        if(iconl > iconMax) {
          iconMax = iconl;
          nodeSetList(i) = nodeID(arma::find(isLeafTemp2 == 1) );
          isLeafTemp = isLeafTemp2;
        }
      }
    }
    iconAll(i) = iconMax;
    i++;
  }
  iconAll(i) = get_ICONTrain(isLeaf, fmat, Smat);
  nodeSetList(i) = nodeID(arma::find(isLeaf == 1) );
}

void Tree::findBeta(arma::vec& iconAll, arma::vec& beta, arma::uvec& sizeTree)
{
  arma::vec alpha(iconAll.n_elem);
  int L = iconAll.n_elem;
  size_t q = 1;
  alpha(0) = 0;
  sizeTree(0) = L;
  while( L > 1 ) {
    arma::vec iconSmallerTree = iconAll.head( L-1 );
    //arma::vec alphaTT = (iconAll(L-1) - iconSmallerTree)/(  arma::regspace<arma::uvec>(L-1,-1,1)  ); // L - 1 to L - (L-1)
    arma::vec alphaTT = (iconAll(L-1) - iconSmallerTree)/(pow(double(L),1) - arma::pow(arma::regspace<arma::vec>(1,L-1),1)  ); // L - 1 to L - (L-1)
    alpha(q) = alphaTT.min();
    sizeTree(q) = alphaTT.index_min() + 1;
    L = sizeTree(q);
    q++;
  }
  if(q < iconAll.n_elem) {
    sizeTree.shed_rows(q, iconAll.n_elem-1);
    alpha.shed_rows(q, iconAll.n_elem-1);
    beta.shed_rows(q, iconAll.n_elem-1);
  }
  for(size_t i = 0; i < alpha.n_elem; i++) {
    if(i < alpha.n_elem - 1) {
      beta(i) = sqrt( alpha(i)*alpha(i+1)  );
    } else {
      beta(i) = alpha(i);
    }
  }
}
#include <memory>
#include <RcppArmadillo.h>

#include "TreeGrow.h"
#include "globals.h"


double TreeGrow::get_ICONTest(const arma::uvec& isLeafTemp,
                              const arma::mat& fmat, const arma::umat& Smat,
                              const arma::mat& fmat2, const arma::umat& Smat2) const
{
  arma::uvec leafTemp = arma::find(isLeafTemp==1);
  int numLeafTemp = leafTemp.n_elem;
  arma::vec icon = arma::zeros<arma::vec>(K);
  arma::vec fSum = arma::zeros<arma::vec>(K);
  arma::vec SSum = arma::zeros<arma::vec>(K);
  for(size_t k = 0; k != K; k++) {
    for(int i = 0; i < numLeafTemp; i++) {
      int li = leafTemp(i);
      double fi = fmat(k,li);
      double Si = Smat(k,li);
      double fi2 = fmat2(k,li);
      double Si2 = Smat2(k,li);
      for(int j = 0; j < i; j++) {
        int lj = leafTemp(j);
        if(fi * Smat(k,lj) > fmat(k,lj) * Si) {
          icon(k) += fi2 * Smat2(k,lj);
        } else {
          icon(k) += fmat2(k,lj) * Si2;
        }
      }
      fSum(k) += fi2;
      SSum(k) += Si2;
      icon(k) += fi2 * Si2 / 2;
    }
  }
  icon = icon / (fSum % SSum);
  arma::uvec zeroDen = arma::find(fSum == 0 || SSum == 0);
  icon ( zeroDen ).ones();
  icon ( zeroDen ) = icon ( zeroDen )/2; // set icon(t_k) to 0.5 if denominator is 0
  return arma::sum(icon) / K;
}

std::shared_ptr<Tree> TreeGrow::trainCV(const arma::umat& mat1Z,
                                        const arma::mat& mat1f,
                                        const arma::field<arma::umat>& mat2Zf,
                                        const arma::umat& range0,
                                        const arma::uvec& e) const
{
  arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
  arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
  std::shared_ptr<Tree> tr = grow(mat1Z, mat1f, mat2Zf, range0, fmat, Smat, e);
  if (NUM_FOLD > 1) {
    const arma::uvec& isl = tr->get_isLeaf();
    uint numLeaf = arma::sum(isl);
    arma::field<arma::uvec> nodeSetList(numLeaf);
    arma::vec iconAll(numLeaf);
    tr->findOptimalSizekSubtree(fmat, Smat, iconAll, nodeSetList, numLeaf);
    // Initialization
    arma::uvec sizeTree = arma::regspace<arma::uvec>(1,iconAll.n_elem);
    arma::vec beta(iconAll.n_elem);
    Tree::findBeta(iconAll, beta, sizeTree); 
    arma::vec iconBeta = prune(beta, mat1Z, mat1f, mat2Zf, range0, e);
    uint qo = iconBeta.index_max();
    arma::uvec nodeSetFinal = nodeSetList(sizeTree(qo)-1);
    tr->cut(nodeSetFinal);
  }
  return tr;
}

std::shared_ptr<Tree> TreeGrow::grow(const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf,
                                     const arma::umat& range0,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = mat1Z.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec parents = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K, MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  size_t ndcount = 0;
  size_t countsp = 0;
  int end = 0;
  while(end == 0) {
    end = split(mat1Z, mat1f, mat2Zf,
                left_childs, right_childs,
                split_vars, split_values, isLeaf,
                parents, fmat, Smat, 
                ranges, nodeSampleY, nodeSample,
                countsp, ndcount, e);
    if(ndcount >= MAX_NODE - 2) {
      isLeaf( arma::find(left_childs == 0) ).ones();
      break;
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
  std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                    right_childs(nonEmpty),
                                    split_vars(nonEmpty),
                                    split_values(nonEmpty),
                                    isLeaf(nonEmpty),
                                    parents(nonEmpty)  ));
  if( ndcount + 1 <=  MAX_NODE-1 ) {
    arma::uvec Empty = arma::regspace<arma::uvec>(ndcount+1, MAX_NODE-1);
    fmat.shed_cols(Empty);
    Smat.shed_cols(Empty);
  }
  if(ndcount > 1) {
    fmat.col(0) = fmat.col(1) + fmat.col(2);
    Smat.col(0) = Smat.col(1) + Smat.col(2);
  } else {
    // no need to prune if there is only one node
  }
  return tr;
}



std::shared_ptr<Tree> TreeGrow::grow(const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf,
                                     const arma::umat& mat1ZVal,
                                     const arma::mat& mat1fVal,
                                     const arma::field<arma::umat>& mat2ZfVal,
                                     const arma::umat& range0,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     arma::mat& fmat2,
                                     arma::umat& Smat2,
                                     const arma::uvec& e) const
{
  int n = mat1Z.n_rows;
  int P = mat1Z.n_cols;
  arma::ucube ranges = arma::zeros<arma::ucube>(MAX_NODE, P, 2);
  arma::uvec left_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec right_childs = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_vars = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec split_values = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec isLeaf = arma::zeros<arma::uvec>(MAX_NODE);
  arma::uvec parents = arma::zeros<arma::uvec>(MAX_NODE);
  ranges.row(0) = range0.t();
  arma::field<arma::uvec> nodeSampleY(MAX_NODE);
  nodeSampleY(0) = arma::regspace<arma::uvec>(0, n-1);
  arma::field<arma::uvec> nodeSample(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSample(k,0) = arma::regspace<arma::uvec>(0, mat2Zf(k).n_rows-1);
  }
  int nVal = mat1ZVal.n_rows;
  arma::field<arma::uvec> nodeSampleYVal(MAX_NODE);
  nodeSampleYVal(0) = arma::regspace<arma::uvec>(0, nVal-1);
  arma::field<arma::uvec> nodeSampleVal(K,MAX_NODE);
  for(size_t k = 0; k < K; k++) {
    nodeSampleVal(k,0) = arma::regspace<arma::uvec>(0, mat2ZfVal(k).n_rows-1);
  }
  size_t ndcount = 0;
  size_t countsp = 0;
  int end = 0;
  while(end == 0) {
    end = split(mat1Z, mat1f, mat2Zf,mat1ZVal, mat1fVal, mat2ZfVal,
                left_childs, right_childs,
                split_vars, split_values, isLeaf,
                parents, fmat, Smat, fmat2, Smat2, 
                ranges, nodeSampleY, nodeSample, nodeSampleYVal, nodeSampleVal,
                countsp, ndcount, e);
    if(ndcount >= MAX_NODE - 2) {
      isLeaf( arma::find(left_childs == 0) ).ones();
      break;
    }
  }
  arma::uvec nonEmpty = arma::regspace<arma::uvec>(0, ndcount);
  std::shared_ptr<Tree> tr(new Tree(left_childs(nonEmpty),
                                    right_childs(nonEmpty),
                                    split_vars(nonEmpty),
                                    split_values(nonEmpty),
                                    isLeaf(nonEmpty),
                                    parents(nonEmpty)  ));
  if( ndcount + 1 <=  MAX_NODE-1 ) {
    arma::uvec Empty = arma::regspace<arma::uvec>(ndcount+1, MAX_NODE-1);
    fmat.shed_cols(Empty);
    Smat.shed_cols(Empty);
    fmat2.shed_cols(Empty);
    Smat2.shed_cols(Empty);
  }


  if(ndcount > 1)
  {
    fmat.col(0) = fmat.col(1) + fmat.col(2);
    Smat.col(0) = Smat.col(1) + Smat.col(2);
    fmat2.col(0) = fmat2.col(1) + fmat2.col(2);
    Smat2.col(0) = Smat2.col(1) + Smat2.col(2);
  }else{
    // no need to prune if there is only one node
  }
  return tr;
}

arma::vec TreeGrow::prune(arma::vec& beta,
                          const arma::umat& mat1Z,
                          const arma::mat& mat1f,
                          const arma::field<arma::umat>& mat2Zf,
                          const arma::umat& range0,
                          const arma::uvec& e) const
{

  // create folds
  int n = mat1Z.n_rows;
  arma::uvec s = arma::shuffle(arma::regspace<arma::uvec>(0,n-1));
  int quotient = n / NUM_FOLD;
  int remainder = n % NUM_FOLD;
  arma::uvec foldstart = arma::regspace<arma::uvec>(0, quotient, n-1);
  arma::uvec foldend = arma::regspace<arma::uvec>(quotient-1, quotient, n-1);
  arma::uvec::iterator it = foldstart.end();
  arma::uvec::iterator it2 = foldend.end();
  it--;
  it2--;
  while(remainder > 0 ) {
    (*it) = *it + remainder;
    (*it2) = *it2 + remainder;
    it--;
    it2--;
    remainder--;
  }
  // create folds end
  arma::mat iconVal(beta.n_elem, NUM_FOLD);
  for(size_t l = 0; l != NUM_FOLD; l++) {
    // get training/validation set
    arma::uvec valid = s(arma::regspace<arma::uvec>(foldstart(l), foldend(l)));
    arma::uvec trainid;
    if(l == 0) {
      trainid = s(arma::regspace<arma::uvec>(foldend(l) + 1, n - 1));
    }
    else if(l == NUM_FOLD-1) {
      trainid = s(arma::regspace<arma::uvec>(0, foldstart(l) - 1));
    } else {
      arma::uvec trainid1 = arma::regspace<arma::uvec>(0, foldstart(l) - 1);
      arma::uvec trainid2 = arma::regspace<arma::uvec>(foldend(l) + 1, n - 1);
      trainid = s(join_cols(trainid1, trainid2));
    }
    arma::field<arma::umat>  mat2Ztrain(K);
    arma::field<arma::umat>  mat2Zval(K);
    for(size_t k = 0; k != K; k++) {
      arma::umat mat2Zfk = mat2Zf(k);
      arma::ivec idk = arma::conv_to<arma::ivec>::from(trainid + mat2Zfk.n_rows - n);
      arma::ivec idkp = idk( find(idk>=0) );
      arma::uvec idkp2 = arma::conv_to<arma::uvec>::from(idkp);
      mat2Ztrain(k) = mat2Zfk.rows( idkp2  );
      arma::ivec idkv = arma::conv_to<arma::ivec>::from(valid + mat2Zfk.n_rows - n);
      arma::ivec idkpv = idkv( find(idkv>=0) );
      arma::uvec idkpv2 = arma::conv_to<arma::uvec>::from(idkpv);
      mat2Zval(k) = mat2Zfk.rows( idkpv2  );
    }
    // arma::umat mat1Ztrain = mat1Z.rows( trainid );
    // arma::mat mat1ftrain = mat1f.cols( trainid );
    // arma::umat mat1Zval = mat1Z.rows( valid );
    // arma::mat mat1fval = mat1f.cols( valid );
    // arma::uvec etrain = e(trainid);
    // get training/validation set end
    arma::mat fmat = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat = arma::zeros<arma::umat>(K, MAX_NODE);
    arma::mat fmat2 = arma::zeros<arma::mat>(K, MAX_NODE);
    arma::umat Smat2 = arma::zeros<arma::umat>(K, MAX_NODE);
    TreeGrow tg(K, spCriterion, MAX_NODE, MIN_NODE1, MIN_SPLIT1);
    std::shared_ptr<Tree> trl = tg.grow(mat1Z.rows( trainid ),
                                        mat1f.cols( trainid ),
                                        mat2Ztrain,
                                        mat1Z.rows( valid ),
                                        mat1f.cols( valid ),
                                        mat2Zval,
                                        range0,
                                        fmat, Smat, fmat2, Smat2, e(trainid));
    const arma::uvec& il = trl->get_isLeaf();
    uint numLeaf = arma::sum(il);
    arma::field<arma::uvec> nodeSetList(numLeaf);
    arma::vec iconAll(numLeaf);
    if (numLeaf > 1) 
      trl->findOptimalSizekSubtree(fmat, Smat, iconAll, nodeSetList, numLeaf);
    arma::vec sizeTree = arma::regspace<arma::vec>(1,iconAll.n_elem);
    arma::uvec isLeafTemp(il.n_elem);
    for(size_t j = 0; j < beta.n_elem; j++) {
      arma::vec iconbetajAll = iconAll -  beta(j) * arma::pow(sizeTree,1);
      uint opt = iconbetajAll.index_max();
      isLeafTemp.zeros();
      isLeafTemp( nodeSetList(opt) ).ones();
      iconVal(j, l) = TreeGrow::get_ICONTest(isLeafTemp, fmat, Smat, fmat2, Smat2);
    }
  }
  return arma::sum(iconVal, 1)/NUM_FOLD;
}



int TreeGrow::split(const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf, // dat
                    arma::uvec& left_childs,
                    arma::uvec& right_childs,
                    arma::uvec& split_vars,
                    arma::uvec& split_values,
                    arma::uvec& isLeaf,
                    arma::uvec& parents,
                    arma::mat& fmat,
                    arma::umat& Smat,// tree
                    arma::ucube& ranges,
                    arma::field<arma::uvec>& nodeSampleY,
                    arma::field<arma::uvec>& nodeSample,
                    size_t& countsp,
                    size_t& ndcount,
                    const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  while(varsp == -1 && countsp <= ndcount) {
    nd = countsp;
    arma::ivec bestSp(3);
    if( spCriterion == 1 ) {
      bestSp = find_split_DICON(nd,
                                mat1Z, mat1f, mat2Zf,
                                ranges, nodeSampleY, nodeSample,
                                fmat, Smat, ndcount, e);
    } else {
      bestSp = find_split_ICON(nd,
                               mat1Z, mat1f, mat2Zf, isLeaf,
                               ranges, nodeSampleY, nodeSample,
                               fmat, Smat, ndcount, e);
    }
    varsp = bestSp(1);
    cutsp = bestSp(2);
    if(varsp == -1) {
      isLeaf(nd) = 1;
      while(countsp <= ndcount) {
        countsp++;
        if(isLeaf(countsp) == 0) break;
      }
    }
  }
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    parents(ndc1) = nd;
    parents(ndc2) = nd;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*mat1Z.n_rows + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <=cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub >cutsp) );
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub<=cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub>cutsp) );
    }
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) isLeaf(ndc1) = 1;
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) isLeaf(ndc2) = 1;
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount += 2;
    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0) break;
    }
  } else {
    end = 1;
  }
  return end;
}

int TreeGrow::split(const arma::umat& mat1Z,
                    const arma::mat& mat1f,
                    const arma::field<arma::umat>& mat2Zf,
                    const arma::umat& mat1ZVal,
                    const arma::mat& mat1fVal,
                    const arma::field<arma::umat>& mat2ZfVal,// dat
                    arma::uvec& left_childs,
                    arma::uvec& right_childs,
                    arma::uvec& split_vars,
                    arma::uvec& split_values,
                    arma::uvec& isLeaf,
                    arma::uvec& parents,
                    arma::mat& fmat,
                    arma::umat& Smat,
                    arma::mat& fmat2,
                    arma::umat& Smat2,// tree
                    arma::ucube& ranges,
                    arma::field<arma::uvec>& nodeSampleY,
                    arma::field<arma::uvec>& nodeSample,
                    arma::field<arma::uvec>& nodeSampleYVal,
                    arma::field<arma::uvec>& nodeSampleVal,
                    size_t& countsp,
                    size_t& ndcount,
                    const arma::uvec& e) const {
  int end = 0;
  int varsp = -1;
  int cutsp = 0;
  size_t nd = countsp;
  if( spCriterion == 1) {
    while(varsp == -1 && countsp <= ndcount) {
      nd = countsp;
      arma::ivec bestSp(3);
      bestSp = find_split_DICON(nd,
				mat1Z, mat1f, mat2Zf,
				ranges, nodeSampleY, nodeSample,
				fmat, Smat, ndcount, e);
      varsp = bestSp(1);
      cutsp = bestSp(2);
      if(varsp == -1) {
	isLeaf(nd) = 1;
	while(countsp <= ndcount) {
	  countsp++;
	  if(isLeaf(countsp) == 0) break;
	}
      }
    }
  }
  if( spCriterion != 1) {
    while(varsp == -1 && countsp <= ndcount) {
      nd = countsp;
      arma::ivec bestSp(3);
      bestSp = find_split_ICON(nd,
			       mat1Z, mat1f, mat2Zf, isLeaf,
			       ranges, nodeSampleY, nodeSample,
			       fmat, Smat, ndcount, e);
      varsp = bestSp(1);
      cutsp = bestSp(2);
      if(varsp == -1) {
	isLeaf(nd) = 1;
	while(countsp <= ndcount) {
	  countsp++;
	  if(isLeaf(countsp) == 0) break;
	}
      }
    }
  }         
  int n = mat1Z.n_rows;
  if(varsp != -1) {
    split_vars(nd) = varsp;
    split_values(nd) = cutsp;
    arma::uword ndc1 = ndcount + 1;
    arma::uword ndc2 = ndcount + 2;
    left_childs(nd) = ndc1;
    right_childs(nd) = ndc2;
    parents(ndc1) = nd;
    parents(ndc2) = nd;
    arma::uvec nodeSampleYnd = std::move(nodeSampleY(nd));
    arma::uvec zvarspsub = mat1Z( varsp*n + nodeSampleYnd );
    nodeSampleY(ndc1) = nodeSampleYnd( arma::find(zvarspsub <= cutsp) );
    nodeSampleY(ndc2) = nodeSampleYnd( arma::find(zvarspsub > cutsp) );
    arma::uvec nodeSampleYndVal = std::move(nodeSampleYVal(nd));
    arma::uvec zvarspsubVal = mat1ZVal( varsp*mat1ZVal.n_rows + nodeSampleYndVal );
    nodeSampleYVal(ndc1) = nodeSampleYndVal( arma::find(zvarspsubVal <=cutsp) );
    nodeSampleYVal(ndc2) = nodeSampleYndVal( arma::find(zvarspsubVal >cutsp) );
    fmat2.col(ndc1) = arma::sum( mat1fVal.cols( nodeSampleYVal(ndc1) ), 1);
    fmat2.col(ndc2) = arma::sum( mat1fVal.cols( nodeSampleYVal(ndc2) ), 1);
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSampleknd = std::move(nodeSample(k,nd));
      arma::uvec zvarspsub = mat2Zf(k)( varsp*mat2Zf(k).n_rows + nodeSampleknd );
      nodeSample(k, ndc1) = nodeSampleknd( arma::find(zvarspsub <= cutsp) );
      nodeSample(k, ndc2) = nodeSampleknd( arma::find(zvarspsub > cutsp) );
    }
    if(nodeSample(0, ndc1).size() < MIN_SPLIT1) isLeaf(ndc1) = 1;
    if(nodeSample(0, ndc2).size() < MIN_SPLIT1) isLeaf(ndc2) = 1;
    for(size_t k = 0; k < K; k++) {
      arma::uvec nodeSamplekndVal = std::move(nodeSampleVal(k,nd));
      arma::uvec zvarspsubVal = mat2ZfVal(k)( varsp*mat2ZfVal(k).n_rows + nodeSamplekndVal );
      nodeSampleVal(k, ndc1) = nodeSamplekndVal( arma::find(zvarspsubVal<=cutsp) );
      nodeSampleVal(k, ndc2) = nodeSamplekndVal( arma::find(zvarspsubVal>cutsp) );
      Smat2(k,ndc1) = nodeSampleVal(k, ndc1).n_elem;
      Smat2(k,ndc2) = nodeSampleVal(k, ndc2).n_elem;
    }
    ranges.row(ndc1) = ranges.row(nd);
    ranges.row(ndc2) = ranges.row(nd);
    ranges(ndc2,varsp,0) = cutsp+1;
    ranges(ndc1,varsp,1) = cutsp;
    ndcount += 2;

    while(countsp <= ndcount) {
      countsp++;
      if(isLeaf(countsp) == 0)
      {break;}
    }
  } else {
    end = 1;
  }
  return end;
}



arma::ivec TreeGrow::find_split_DICON(size_t nd,
                                      const arma::umat& mat1Z,
                                      const arma::mat& mat1f,
                                      const arma::field<arma::umat>& mat2Zf, // dat
                                      const arma::ucube& ranges,
                                      const arma::field<arma::uvec>& nodeSampleY,
                                      const arma::field<arma::uvec>& nodeSample,
                                      arma::mat& fmat,
                                      arma::umat& Smat,
                                      size_t ndcount,
                                      const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  for(int p = 0; p < P; p++) {
    // arma::uvec zp = mat1Z.col(p);
    // arma::uvec zpsub = zp(nodeSampleY(nd));
    // arma::uvec ind = arma::regspace<arma::uvec>(0,n); // should be n-1
    // arma::uvec indsub = ind(nodeSampleY(nd));
    // arma::uvec indY = indsub(sort_index(zpsub));
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    //indY.print();
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      // arma::umat mat2zp = mat2Zf(k);
      // arma::uvec zp = mat2zp.col(p);
      // arma::uvec zpsub = zp(nodeSample(k,nd));
      // arma::uvec ind = arma::regspace<arma::uvec>(0,zp.size() );
      // arma::uvec indsub = ind(nodeSample(k,nd));
      // indp(k) = indsub(sort_index(zpsub));
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = arma::sum(mat1f.cols(indY),1);
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0), ranges(nd, p, 1));
    int nel = 0;
    // int nelr = arma::sum(e( indY ) );
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        size_t z = mat1Z(indYj, p);
        if(z == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }

     // if( (SLSum.min() < MIN_NODE1 || SRSum.min() < MIN_NODE1) && ( nel < MIN_NODE2 || nelr - nel < MIN_NODE2) )
      if( (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
        //Rcpp::Rcout << "small 1" << SLSum.min() << "a"<<SRSum.min() << "a"<<j <<"a"<< indY.size()-j << "\n";
      } else {
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum) );
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
        fmat.col(ndcount + 1) = fLSum;
        fmat.col(ndcount + 2) = fRSum;
        Smat.col(ndcount + 1) = SLSum;
        Smat.col(ndcount + 2) = SRSum;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}

arma::ivec TreeGrow::find_split_ICON(size_t nd,
                                     const arma::umat& mat1Z,
                                     const arma::mat& mat1f,
                                     const arma::field<arma::umat>& mat2Zf, // dat
                                     const arma::uvec& isLeaf,
                                     arma::ucube& ranges,
                                     arma::field<arma::uvec>& nodeSampleY,
                                     arma::field<arma::uvec>& nodeSample,
                                     arma::mat& fmat,
                                     arma::umat& Smat,
                                     size_t ndcount,
                                     const arma::uvec& e) const {
  int P = mat1Z.n_cols;
  int n = mat1Z.n_rows;
  int varsp = -1;
  int cutsp = 0;
  double dICONmax = 0;
  double dICONTemp = 0;
  arma::mat fmatTerm = fmat.cols(arma::find(isLeaf == 1));
  arma::umat SmatTerm = Smat.cols(arma::find(isLeaf == 1));
  for(int p = 0; p < P; p++) {
    arma::uvec indY = nodeSampleY(nd)( sort_index( mat1Z(p*n + nodeSampleY(nd)) ));
    arma::field<arma::uvec> indp(K);
    arma::uvec SRSum = arma::zeros<arma::uvec>(K);
    for(size_t k = 0; k < K; k++) {
      arma::uvec zpsub = mat2Zf(k)( p*mat2Zf(k).n_rows + nodeSample(k,nd) );
      indp(k) = nodeSample(k,nd)(sort_index(zpsub));
      SRSum(k) = zpsub.size();
    }
    arma::vec fLSum = arma::zeros<arma::vec>(K);
    arma::uvec SLSum = arma::zeros<arma::uvec>(K);
    arma::vec fRSum = (sum(mat1f.cols(indY),1));
    int j = 0;
    arma::uvec jv = arma::zeros<arma::uvec>(K);
    int nj = indY.size();
    int nel = 0;
    // int nelr = arma::sum( e(indY) );
    arma::uvec rangeCut = arma::regspace<arma::uvec>(ranges(nd, p, 0), ranges(nd, p, 1));
    for(auto cu : rangeCut) {
      while(j < nj) {
        int indYj = indY(j);
        size_t z = mat1Z(indYj, p);
        if(z == cu) {
          arma::vec df = mat1f.col( indYj );
          fLSum = fLSum + df;
          fRSum = fRSum - df;
          nel += e( indYj );
          j++;
        } else {
          break;
        }
      }
      for(size_t k = 0; k < K; k++) {
        arma::uvec indpk = indp(k);
        while(jv(k) < indpk.size()) {
          if(mat2Zf(k)( indpk( jv(k) ) , p) == cu) {
            SLSum(k)++;
            SRSum(k)--;
            jv(k)++;
          } else {
            break;
          }
        }
      }
      //if( (SLSum.min() < MIN_NODE1 || SRSum.min() < MIN_NODE1) && ( nel < MIN_NODE2 || nelr-nel < MIN_NODE2) )
      if(  (SLSum(0) < MIN_NODE1 || SRSum(0) < MIN_NODE1) ) {
        dICONTemp = 0;
      } else {
        dICONTemp = arma::sum(abs(fLSum%SRSum - fRSum%SLSum))/2;
        for(size_t i = 0; i < fmatTerm.n_cols; i++) {
          if( i != nd ) {
            for(size_t k = 0; k < K; k++) {
              //  DOUBLE CHECK !!!
              if( (fLSum(k) + fRSum(k)) * SmatTerm(k,i) <= fmatTerm(k,i) * (SLSum(k) + SRSum(k)) ) {
                if( fmatTerm(k,i)*SLSum(k) <  fLSum(k)*SmatTerm(k,i) )
                {// lambdaR < lambdand <= lambdai < lambdaL
                  dICONTemp = dICONTemp +  fLSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SLSum(k) ;
                }
                if( fmatTerm(k,i)*SRSum(k) <  fRSum(k)*SmatTerm(k,i) )
                {// lambdaL < lambdand <= lambdai < lambdaR
                  dICONTemp = dICONTemp + fRSum(k)*SmatTerm(k,i) - fmatTerm(k,i)*SRSum(k) ;
                }
              } else {
                if( fmatTerm(k,i)*SLSum(k) >  fLSum(k)*SmatTerm(k,i) )
                {// lambdaL < lambdai < lambdand < lambdaR
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SLSum(k) - fLSum(k)*SmatTerm(k,i) ;
                }
                if( fmatTerm(k,i)*SRSum(k) >  fRSum(k)*SmatTerm(k,i) )
                {// lambdaR < lambdai < lambdand < lambdaL
                  dICONTemp = dICONTemp + fmatTerm(k,i)*SRSum(k) - fRSum(k)*SmatTerm(k,i) ;
                }
              }
            }
          }
        }
      }
      if(dICONTemp>dICONmax) {
        dICONmax = dICONTemp;
        varsp = p;
        cutsp = cu;
        fmat.col(ndcount + 1) = fLSum;
        fmat.col(ndcount + 2) = fRSum;
        Smat.col(ndcount + 1) = SLSum;
        Smat.col(ndcount + 2) = SRSum;
      }
    }
  }
  arma::ivec vecsp(3);
  if(varsp == -1) {
    vecsp(0) = 0;
    vecsp(1) = -1;
    vecsp(2) = 0;
  } else {
    vecsp(0) = 1;
    vecsp(1) = varsp;
    vecsp(2) = cutsp;
  }
  return vecsp;
}
#include "TreePrediction.h"
#include "globals.h"

void TreePrediction::transformZ(const arma::mat& z,
                                arma::umat& z2,
                                const arma::mat& dat,
                                const arma::uvec& e,
                                const arma::vec& breaks,
                                const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = dat.col(p);
      // BE CAREFUL
      int j = 0;
      double z2ecdf;
      for(int i = 0; i < n; i++) {
        if(e(i) == 1) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,i) )- zpref.begin())/(n-i+0.0)  ;
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1;
          j++;
        }
        ind = ind + 1;
        if(ind.n_elem > 0){ ind.shed_row(0);}
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}

void TreePrediction::transformZH(const arma::mat& z,  // z on tg
                                 const arma::vec& tg, //grid of time point
                                 arma::umat& z2,
                                 const arma::mat& dat, // dat is predictor matrix
                                 const arma::vec& y,
                                 const arma::uvec& e,
                                 const arma::vec& breaks,
                                 const arma::uvec& disc)
{
  int P = z.n_rows;
  int n = e.n_elem;
  int NG = tg.n_elem;
  arma::vec::const_iterator bb = breaks.begin();
  arma::vec::const_iterator be = breaks.end();
  for(int p = 0; p < P; p++) {
    if(disc(p) == 0) {
      arma::uvec ind = arma::cumsum(arma::regspace<arma::uvec>(0,n-1));
      arma::vec zp = dat.col(p);
      // BE CAREFUL
      int j = 0;
      double z2ecdf;
      for(int i = 0; i < n-1; i++) {
        if( y(i) <= tg(j) && tg(j) < y(i+1) ) {
          arma::vec zpref = arma::sort(zp.elem( ind ));
          z2ecdf = (std::lower_bound(zpref.begin(), zpref.end(), z(p,j) )- zpref.begin())/(n-i+0.0)  ;
          z2(p,j) = std::distance(bb, std::upper_bound(bb, be, z2ecdf ) ) + 1;
          j++;
          if(j == NG) break;
        }
        ind = ind + 1;
        if(ind.n_elem > 0){ ind.shed_row(0);}
      }
    } else {
      z2.row(p) = arma::conv_to<arma::urowvec>::from(z.row(p));
    }
  }
}


TreePrediction::TreePrediction(const arma::umat& zy,
                               const arma::field<arma::umat>& zt,
                               const arma::uvec& vars,
			       const arma::uvec& values,
			       const arma::uvec& lcs,
			       const arma::uvec& rcs,
			       const arma::uvec& il)
{
  arma::uvec ndy(zy.n_cols);
  for(size_t it = 0; it != zy.n_cols; it++) {
    arma::uvec zyi = (zy.col(it));
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zyi(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    ndy(it) = k;
  }
  size_t nT = zt.size();
  //arma::field<arma::uvec> ndst(nT);
  int nNd = arma::accu(il);
  //Rcpp::Rcout << arma::find(il == 1);
  arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
  //Rcpp::Rcout << MAX_NODE;
  arma::uvec tnd2 = tnd(arma::find(il == 1));
  arma::umat ndsz = arma::zeros<arma::umat>(nNd, nT);
  arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
  tnd3( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
  for(arma::uword c=0; c < nT; c++) {
    arma::umat m = zt(c);
    int it_end = m.n_cols;
    arma::uvec ndt(it_end);
    for(int it = 0; it != it_end; it++) {
      arma::uvec zti = (m.col(it));
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      int k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zti(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndt(it) = k;
      //Rcpp::Rcout << k;
      ndsz(tnd3(k),c)++;
    }
    //ndst(c) = ndt;
  }
  this->nodeLabel = ndy;
  this->nodeSize = ndsz;
  this->tnd3 = tnd3;
}


arma::vec TreePrediction::getSurvival(const arma::umat& zt2,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::uvec& vars,
				      const arma::uvec& values,
				      const arma::uvec& lcs,
				      const arma::uvec& rcs,
				      const arma::uvec& il)
{
  arma::vec y2 = y( arma::find(e == 1));
  int NE = y2.n_elem;
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(int i = 0; i < NE; i++) {
    arma::uvec zY = zt2.col(i);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zY(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    size_t nd = k;
    if(nodeLabel(i) == nd) {
      w2(i) = w2(i) + 1.0/nodeSize(tnd3(nd),i);
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( arma::find(e == 1)) = w2;
  return exp(-cumsum(w));
}

arma::vec TreePrediction::getSurvival(const arma::umat& zt2,
                                      const arma::vec& y,
                                      const arma::uvec& e,
                                      const arma::umat& nodeSize0,
                                      const arma::uvec& nodeLabel0,
                                      const arma::uvec& tnd30,
                                      const arma::umat& treeMat)
{
  arma::vec y2 = y( arma::find(e == 1));
  arma::uvec vars = treeMat.col(0);
  arma::uvec values = treeMat.col(1);
  arma::uvec lcs = treeMat.col(2);
  arma::uvec rcs = treeMat.col(3);
  arma::uvec il = treeMat.col(4);
  int NE = y2.n_elem;
  arma::vec w2 = arma::zeros<arma::vec>(NE);
  for(int i = 0; i < NE; i++) {
    arma::uvec zY = zt2.col(i);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    size_t k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zY(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    // size_t nd = k;
    if(nodeLabel0(i) == k) {
      w2(i) = w2(i) + 1.0/nodeSize0(tnd30(k),i);
    }
  }
  arma::vec w = arma::zeros<arma::vec>(y.n_elem);
  w( arma::find(e == 1)) = w2;
  return exp(-cumsum(w));
}


arma::vec TreePrediction::getHazard(const arma::umat& ztvec,
                                    const arma::vec& tg,
                                    const arma::vec& y,
                                    const arma::uvec& e,
                                    const arma::mat& fy2,
                                    const double h,
                                    const arma::umat& nodeSize0,
                                    const arma::uvec& nodeLabel0,
                                    const arma::uvec& tnd30,
                                    const arma::umat& treeMat)
{
  arma::vec y2 = y( arma::find(e == 1));
  arma::uvec vars = treeMat.col(0);
  arma::uvec values = treeMat.col(1);
  arma::uvec lcs = treeMat.col(2);
  arma::uvec rcs = treeMat.col(3);
  arma::uvec il = treeMat.col(4);
  int NG = tg.n_elem;
  int NE = y2.n_elem;
  arma::vec hz = arma::zeros<arma::vec>(NG);
  for(int j = 0; j < NG; j++) {
    arma::uvec zj = ztvec.col(j);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    size_t k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zj(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    arma::vec v0 = arma::zeros<arma::vec>(NE);
    arma::vec v1 = arma::zeros<arma::vec>(NE);
    for(int i = 0; i < NE; i++) {
      if( y2(i) >= tg(j) - h && y2(i) <= tg(j) + h && nodeLabel0(i) == k) {
        //hz(j) = hz(j) + fy2(i,j)/nodeSize0(tnd30(k),i);
        v0(i)++;
      }
      v1(i) += nodeSize0(tnd30(k),i);
    }
    v1(arma::find(v1 == 0)).ones();
    hz(j) = arma::sum( fy2.col(j) %  v0 / v1 );
  }
  return hz;
}

TreePrediction::TreePrediction(const Data2* dat2,
                               const arma::uvec& vars,
			       const arma::uvec& values,
			       const arma::uvec& lcs,
			       const arma::uvec& rcs,
			       const arma::uvec& il)
{
  arma::umat zy = dat2->get_zy();
  arma::field<arma::umat> zt = dat2->get_zt();
  arma::uvec ndy(zy.n_cols);
  for(size_t it = 0; it != zy.n_cols; it++) {
    arma::uvec zyi = zy.col(it);
    int isl = 0;
    int varsp = 0;
    size_t cutsp = 0;
    int k = 0;
    while(isl == 0) {
      varsp = vars(k);
      cutsp = values(k);
      if(zyi(varsp) > cutsp) {
        k = rcs(k);
      } else {
        k = lcs(k);
      }
      isl = il(k);
    }
    ndy(it) = k;
  }
  size_t nT = zt.size();
  //arma::field<arma::uvec> ndst(nT);
  int nNd = arma::accu(il);
  arma::uvec tnd = arma::regspace<arma::uvec>(0, il.n_elem-1);
  arma::uvec tnd2 = tnd(arma::find(il == 1));
  arma::umat ndsz = arma::zeros<arma::umat>(nNd, nT);
  arma::uvec tnd3 = arma::zeros<arma::uvec>(il.n_elem);
  tnd3.elem( tnd2 ) = arma::regspace<arma::uvec>(0, nNd-1);
  for(arma::uword c=0; c < nT; c++) {
    arma::umat m = zt(c);
    int it_end = m.n_cols;
    arma::uvec ndt(it_end);
    for(int it = 0; it != it_end; it++) {
      arma::uvec zti = (m.col(it));
      int isl = 0;
      int varsp = 0;
      size_t cutsp = 0;
      int k = 0;
      while(isl == 0) {
        varsp = vars(k);
        cutsp = values(k);
        if(zti(varsp) > cutsp) {
          k = rcs(k);
        } else {
          k = lcs(k);
        }
        isl = il(k);
      }
      ndt(it) = k;
      ndsz(tnd3(k),c)++;
    }
  }
  this->nodeLabel = ndy;
  this->nodeSize = ndsz;
  this->tnd3 = tnd3;
}
