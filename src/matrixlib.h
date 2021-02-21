#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

float determinant(float* M, int n);
float cofactor(float* M, int n, int i, int j);
float minor(float* M, int n, int i, int j);
float* submatrix(float* M, int n, int i, int j);
float* adjoint(float* M, int n);
float* inverse(float* M, int n);
float* multiply(float* A, float* v, int m, int n);

void test_adjoint();
void test_inverse();
void test_determinant();
void test_submatrix();

float determinant(float* M, int n) {
  if (n==1) return M[0];
  int sum = 0;
  // cofactor expansion along the first row
  for (int j=0; j<n; ++j) {
    sum += M[j] * cofactor(M, n, 0, j);
  }
  return sum;
}

float cofactor(float* M, int n, int i, int j) {
  int sign = (i+j)%2 == 0 ? 1 : -1;
  return sign * minor(M, n, i, j);
}

float minor(float* M, int n, int i, int j) {
  float* sub = submatrix(M, n, i, j);
  float minor = determinant(sub, n-1);
  free(sub);
  return minor;
}

float* submatrix(float* M, int n, int i, int j) {
  float* C = malloc((n-1)*(n-1)*sizeof(float));
  for (int k=0; k<i; ++k) {
    memcpy(C + (n-1)*k, M + n*k, j*sizeof(float));
    memcpy(C + j + (n-1)*k, M + j + 1 + n*k, (n-j-1)*sizeof(float));
  }
  for (int k=i+1; k<n; ++k) {
    memcpy(C + (n-1)*(k-1), M + n*k, j*sizeof(float));
    memcpy(C + j + (n-1)*(k-1), M + j + 1 + n*k, (n-j-1)*sizeof(float));
  }
  return C;
}

void test_submatrix() {
  int n = 3;
  float M[9] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  };
  float* S;
  for (int r=0; r<n-1; ++r) {
    for (int c=0; c<n-1; ++c) {
      S = submatrix(M, n, r, c);
      for (int i=0; i<n-1; ++i) {
        for (int j=0; j<n-1; ++j) {
          assert( S[j + (n-1)*i] == M[(j<c ? j : j+1) + n*(i<r ? i : i+1)] );
        }
      }
      free(S);
    }
  }
}

void test_determinant() {
  int n = 4;
  float M[16] = {
    1, 3, 5, 9,
    1, 3, 1, 7,
    4, 3, 9, 7,
    5, 2, 0, 9
  };
  float det = determinant(M, n);
  assert(det == -376);
}

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
  if (det==0) return NULL;
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

float* multiply(float* M, float* v, int nrow, int ncol) {
  float* u = malloc(nrow*sizeof(float));
  for (int i=0; i<nrow; ++i) {
    u[i] = 0;
    for (int j=0; j<ncol; ++j) {
      u[j] += M[i*ncol+j] * v[j];
    }
  }
  return u;
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

//int main() { test_adjoint(); test_inverse(); test_submatrix(); test_determinant(); }