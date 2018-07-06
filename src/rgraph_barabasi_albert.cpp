#include <RcppArmadillo.h>
using namespace Rcpp;

//' Barabasi-Albert
//' @param t Time steps.
//' @param m Number of ties added per node.
//' @param alpha Degree
//' @param a Baseline probability
//' @param loop Logical. When `TRUE` it allows for auto-links (loops).
//' @export
// [[Rcpp::export(name = "rgraph_ba")]]
arma::imat rgraph_barabasi_albert(
    int t,
    int m,
    double alpha = 1.0,
    double a     = 1.0,
    bool loop    = false
  ) {

  arma::imat edgelist(m*(t - 1), 2);

  // Vector of indegree
  std::vector< double > k(t + 1);
  k.at(0) = 0;

  double sum_of_k = a, u;

  int j, idx = -1;
  for (int i = 1; i < t; i++) {

    // Initializing appealingness, and updating the denominator of the prob
    k.at(i) = 0;

    // If loop, then add it before
    if (loop)
      sum_of_k += a;

    // Drawing number
    for (int ii = 0; ii < m; ii++) {

      u = unif_rand()*sum_of_k;

      j = -1;
      double cumprob = 0.0;

      while (++j < i) {
        cumprob += powf(k.at(j), alpha) + a;
        if (u < cumprob)
          break;
      }

      if (!loop && i == j)
        j--;

      edgelist.at(++idx, 0) = i;
      edgelist.at(idx, 1)   = j;

      // Updating sum_of_k
      sum_of_k -= powf(k.at(j), 1.0/alpha);
      k.at(j) += 1.0;
      sum_of_k += powf(k.at(j), alpha);

    }

    // If no loop, then add it after
    if (!loop)
      sum_of_k += a;

  }

  return edgelist;

}
