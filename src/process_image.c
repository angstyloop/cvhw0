#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "image.h"

// half-open intervals by default [l, r)

int clamp_to_range(int x, int l, int r) {
  int x_;
  if (x < l)
    x_ = l;
  else if (x >= r)
    x_ = r-1;
  else
    x_ = x;
  return x_;
}

int clamp(int x, int r) {
  return clamp_to_range(x, 0, r);
}

// half open range [l, r)
bool in_range(int x, int l, int r) {
  return l<=x && x<r;
}

// When getting pixels, use "clamp" strategy - for out-of-bounds coords, return the nearest boundary value.
float get_pixel(image im, int x, int y, int c)
{
  return im.data[clamp(x, im.w) + im.w*clamp(y, im.h) + im.w*im.h*clamp(c, 3)];
}

// When setting pixels, just don't set anything if out of bounds.
void set_pixel(image im, int x, int y, int c, float v)
{
  if (in_range(x, 0, im.w) && in_range(y, 0, im.h) && in_range(c, 0, im.c))
    im.data[x + im.w*y + im.w*im.h*c] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w*im.h*im.c*sizeof(float));
    return copy;
}

float luma(float r, float g, float b) {
  return 0.299*r + 0.587*g + 114*b;
}

image rgb_to_grayscale(image im)
{
    float sum;
    float luma[] = {1,2,3}; // coefficients in a sum over the gamma-encoded channels
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int x=0; x<im.w; x++) {
      for (int y=0; y<im.h; y++) {
        sum = 0;
        for (int c=0; c<im.c; c++) {
          sum += luma[c] * get_pixel(im, x, y, c);
        }
        set_pixel(gray, x, y, 0, sum); // set the single channel to the luma value
      }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
}

void clamp_image(image im)
{
    // TODO Fill this in
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
}
