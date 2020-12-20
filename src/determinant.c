#include "determinant.h"

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

//int main() { test_submatrix(); test_determinant(); }