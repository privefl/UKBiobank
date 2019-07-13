#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dist2(const NumericMatrix& U,
                    const NumericVector& x) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector dist(m);

  for (int j = 0; j < m; j++) {
    double d = 0;
    for (int i = 0; i < n; i++) {
      double u = U(i, j) - x[i];
      d += u*u;
    }
    dist[j] = d;
  }

  return dist;
}

// [[Rcpp::export]]
NumericVector dist_sort(const NumericMatrix& U,
                        const NumericVector& x,
                        int nfirst) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector dist(m);

  for (int j = 0; j < m; j++) {
    double d = 0;
    for (int i = 0; i < n; i++) {
      double u = U(i, j) - x[i];
      d += u*u;
    }
    dist[j] = d;
  }

  std::partial_sort(dist.begin(), dist.begin() + nfirst, dist.end());

  return dist;
}

// [[Rcpp::export]]
double dist_kNN(const NumericMatrix& U,
                const NumericVector& x,
                int nfirst) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector dist(m);

  for (int j = 0; j < m; j++) {
    double d = 0;
    for (int i = 0; i < n; i++) {
      double u = U(i, j) - x[i];
      d += u*u;
    }
    dist[j] = d;
  }

  std::partial_sort(dist.begin(), dist.begin() + nfirst, dist.end());

  return std::accumulate(dist.begin(), dist.begin() + nfirst, 0.0);
}

// [[Rcpp::export]]
NumericVector dist_kNN_all(const NumericMatrix& U,
                           int nfirst) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector dist(m), res(m);

  for (int k = 0; k < m; k++) {

    for (int j = 0; j < m; j++) {
      double d = 0;
      for (int i = 0; i < n; i++) {
        double u = U(i, j) - U(i, k);
        d += u*u;
      }
      dist[j] = d;
    }

    std::partial_sort(dist.begin(), dist.begin() + nfirst, dist.end());

    res[k] = std::accumulate(dist.begin(), dist.begin() + nfirst, 0.0);
  }

  return res;
}


// [[Rcpp::export]]
NumericVector dist_kNN_all_clever(const NumericMatrix& U,
                                  int nfirst) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector min_dist(nfirst), res(m);

  for (int k = 0; k < m; k++) {

    for (int l = 0; l < nfirst; l++) min_dist[l] = R_PosInf;

    int pos_max = which_max(min_dist);
    double d_max = min_dist[pos_max];

    bool go_up = true, go_down = true;
    int j, l = 1;

    while (go_up || go_down) {

      j = k + l;
      go_up = (j < m);
      if (go_up) {
        double u, d = 0;
        int i = 0;
        while (i < n && d < d_max) {
          u = U(i, j) - U(i, k);
          d += u*u;
          i++;
        }
        if (d < d_max) { // new max amongst min distances
          min_dist[pos_max] = d;
          pos_max = which_max(min_dist);
          d_max = min_dist[pos_max];
        }
      }

      j = k - l;
      go_down = (j >= 0);
      if (go_down) {
        double u, d = 0;
        int i = 0;
        while (i < n && d < d_max) {
          u = U(i, j) - U(i, k);
          d += u*u;
          i++;
        }
        if (d < d_max) { // new max amongst min distances
          min_dist[pos_max] = d;
          pos_max = which_max(min_dist);
          d_max = min_dist[pos_max];
        }
      }

      l++;

    }

    res[k] = Rcpp::sum(min_dist);
  }

  return res;
}

// [[Rcpp::export]]
NumericMatrix dist_kNN_all_clever2(const NumericMatrix& U,
                                   int nfirst) {

  int n = U.nrow();
  int m = U.ncol();

  NumericVector min_dist(nfirst);
  NumericMatrix res(nfirst, m);

  for (int k = 0; k < m; k++) {

    for (int l = 0; l < nfirst; l++) min_dist[l] = R_PosInf;

    int pos_max = which_max(min_dist);
    double d_max = min_dist[pos_max];

    bool go_up = true, go_down = true;
    int j, l = 1;

    while (go_up || go_down) {

      j = k + l;
      go_up = (j < m);
      if (go_up) {
        double u, d = 0;
        int i = 0;
        while (i < n && d < d_max) {
          u = U(i, j) - U(i, k);
          d += u*u;
          i++;
        }
        if (d < d_max) { // new max amongst min distances
          min_dist[pos_max] = d;
          pos_max = which_max(min_dist);
          d_max = min_dist[pos_max];
        }
      }

      j = k - l;
      go_down = (j >= 0);
      if (go_down) {
        double u, d = 0;
        int i = 0;
        while (i < n && d < d_max) {
          u = U(i, j) - U(i, k);
          d += u*u;
          i++;
        }
        if (d < d_max) { // new max amongst min distances
          min_dist[pos_max] = d;
          pos_max = which_max(min_dist);
          d_max = min_dist[pos_max];
        }
      }

      l++;

    }

    std::sort(min_dist.begin(), min_dist.end());
    for (int l = 0; l < nfirst; l++) res(l, k) = min_dist[l];
  }

  return res;
}
