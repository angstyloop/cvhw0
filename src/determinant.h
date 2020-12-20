#ifndef ASSERT_H
#define ASSERT_H
#include <assert.h>
#endif

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif

#ifndef DETERMINANT_H
#define DETERMINANT_H

float determinant(float* M, int n);
float cofactor(float* M, int n, int i, int j);
float minor(float* M, int n, int i, int j);
float* submatrix(float* M, int n, int i, int j);
void test_determinant();
void test_submatrix();

#endif //DETERMINANT_H
