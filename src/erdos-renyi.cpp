#include <RcppArmadillo.h>
using namespace Rcpp;

//' Erdős-Renyi
//' @param n Size.
//' @param m Number of edges in the `G(n,m)` model.
//' @param loop Logical. When `TRUE`, no loops.
//' @name Erdos-Renyi
//' @export
// [[Rcpp::export(name = "rgraph_er_gnm")]]
arma::imat rgraph_erdos_renyi_gnm(double n, double m, bool loop = false) {

  arma::imat edgelist(m, 2);
  int ego, alter;
  for (int i = 0; i < m; i++) {

    ego   = floor(unif_rand()*n) + 1;
    alter = floor(unif_rand()*n) + 1;

    if (!loop)
      while (alter == ego)
        alter = floor(unif_rand()*n) + 1;

    edgelist.at(i, 0) = ego;
    edgelist.at(i, 1) = alter;

  }

  return edgelist;

}

//' @rdname Erdos-Renyi
//' @export
//' @param p Probability as in the `G(n, p)` model.
// [[Rcpp::export(name = "rgraph_er_gnp")]]
arma::imat rgraph_erdos_renyi_gnp(double n, double p, bool loop = false) {

  int m = R::rbinom(n*(n - (double) loop), p);

  return rgraph_erdos_renyi_gnm(n, m, loop);

}
