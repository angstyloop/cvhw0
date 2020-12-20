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

#include "determinant.h"

#ifndef INVERSE_H
#define INVERSE_H

float* adjoint(float* M, int n);
float* inverse(float* M, int n);
void test_adjoint();
void test_inverse();

#endif //VISION_HW0_INVERSE_H
