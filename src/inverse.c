#include "determinant.h"

float* adjoint(float* M, int n) {
  float* adj = malloc(n*n*sizeof(float));
  if (n==1) {
    adj[0] = 1;
  } else {
    for (int i=0; i<n; ++i) {
      for (int j=0; j<n; ++j) {
        adj[i+n*j] = cofactor(M, n, i, j);
      }
    }
  }
  return adj;
}

float* inverse(float* M, int n) {
  float det = determinant(M, n);
  if (det==0) {
    printf("Matrix not invertible.");
    assert(0);
  }
  float* adj = adjoint(M, n);
  float* inv = malloc(n*n*sizeof(float));
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      inv[j+n*i] = adj[j+n*i] / det;
    }
  }
  free(adj);
  return inv;
}

void test_adjoint() {
  int n = 3;
  float M[9] = {
    3, 1,  1,
    1, 3, -1,
    2, 4,  1
  };
  float* actual_adj = adjoint(M, n);
  float expected_adj[9] = {
     7,   3, -4,
    -3,   1,  4,
    -2, -10,  8
  };
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      assert(actual_adj[j+n*i] == expected_adj[j+n*i]);
    }
  }
  free(actual_adj);
}

void test_inverse() {
  int n = 3;
  float M[9] = {
    1, 2, 3,
    4, 5, 6,
    7, 2, 9
  };
  float* actual_inv = inverse(M, n);
  float expected_inv[9] = {
    -11/12., 1/3., 1/12.,
    -1/6., 1/3., -1/6.,
    3/4., -1/3., 1/12.
  };
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      assert(actual_inv[j+n*i] == expected_inv[j+n*i]);
    }
  }
  free(actual_inv);
}

//int main() { test_adjoint(); test_inverse(); }