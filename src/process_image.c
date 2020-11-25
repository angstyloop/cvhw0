#include <string.h>
#include <assert.h>
#include <math.h>
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

typedef enum {false, true} boolean;

// half open range [l, r)
boolean in_range(int x, int l, int r) {
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

// convention in video image processing
float luma(float r, float g, float b) {
  return 0.299*r + 0.587*g + 114*b;
}

image rgb_to_grayscale(image im)
{
    float sum;
    float luma[] = {0.299, 0.587, 0.114}; // coefficients in a sum over the gamma-encoded channels
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

// scale the intensity of channel c in each pixel by a constant amount
void scale_image(image im, int c, float v) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      set_pixel(im, x, y, c, v * get_pixel(im, x, y, c));
    }
  }
}

// shift the intensity of channel c in each pixel by a constant amount
void shift_image(image im, int c, float v)
{
    for (int x=0; x<im.w; x++) {
      for (int y=0; y<im.h; y++) {
        set_pixel(im, x, y, c, get_pixel(im, x, y, c) + v);
      }
    }
}

// clamp pixel values to [0,1]
void clamp_image(image im)
{
  float v;
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      for (int z=0; z<im.c; z++) {
        v = get_pixel(im, x, y, z);
        if (v<0)
          set_pixel(im, x, y, z, 0);
        else if (v>1)
          set_pixel(im, x, y, z, 1);
      }
    }
  }
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

// compute hue from value, chroma, and (R,G,B) values
// convention of our implementation is "C==0 => H==0"
float hue(float V, float C, float R, float G, float B)
{
  float H_, H;
  H_ = H = 0;
  if (C > 0) {
    if (V == R)
      H_ = (G - B)/C;
    else if (V == G)
      H_ = (B - R)/C + 2;
    else if (V == B)
      H_ = (R - G)/C + 4;
  }
  if (H_ < 0)
    H = H_/6 + 1;
  else
    H = H_/6;
  return H;
}

void rgb_to_hsv(image im)
{
  float R, G, B, H, S, V, C, m;
  R = G = B = H = V = C = m = 0;
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      // read rgb values
      R = get_pixel(im, x, y, 0);
      G = get_pixel(im, x, y, 1);
      B = get_pixel(im, x, y, 2);

      // computed HSV values
      V = three_way_max(R, G, B);
      m = three_way_min(R, G, B);
      C = V - m;
      S = V>0 ? C/V : 0;
      H = hue(V, C, R, G, B);

      // write hsv values instead of rgb
      set_pixel(im, x, y, 0, H);
      set_pixel(im, x, y, 1, S);
      set_pixel(im, x, y, 2, V);
    }
  }
}

void hsv_to_rgb(image im)
{
  float R_, G_, B_, R, G, B, H, S, V, C, m, X;
  R_ = G_ = B_ = R = G = B = H = S = V = C = m = 0;
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      H = get_pixel(im, x, y, 0);
      S = get_pixel(im, x, y, 1);
      V = get_pixel(im, x, y, 2);
      C = V * S;
      m = V - C;
      X = C * (float)( 1 - fabs(fmod(6*H, 2) - 1) );
      if (0/6.<=H && H<1/6.) {
        R_ = C;
        G_ = X;
        B_ = 0;
      } else if (1/6.<=H && H<2/6.) {
        R_ = X;
        G_ = C;
        B_ = 0;
      } else if (2/6.<=H && H<3/6.) {
        R_ = 0;
        G_ = C;
        B_ = X;
      } else if (3/6.<=H && H<4/6.) {
        R_ = 0;
        G_ = X;
        B_ = C;
      } else if (4/6.<=H && H<5/6.) {
        R_ = X;
        G_ = 0;
        B_ = C;
      } else if (5/6.<=H && H<6/6.) {
        R_ = C;
        G_ = 0;
        B_ = X;
      }

      R = (R_ + m);
      G = (G_ + m);
      B = (B_ + m);

      set_pixel(im, x, y, 0, R);
      set_pixel(im, x, y, 1, G);
      set_pixel(im, x, y, 2, B);
    }
  }
}

float gamma_compress_value(float u) {
  return (u<=0.0031308) ? (float)(12.92*u) : (float)(1.055*pow(u, 1/2.4)-0.055);
}

float gamma_decompress_value(float u) {
  return (u<=0.04045) ? (float)(u/12.92) : (float)pow((u+0.055)/1.055, 2.4);
}

float get_pixel_compressed(image im, int x, int y, int z) {
  return gamma_compress_value(get_pixel(im, x, y, z));
}

float get_pixel_decompressed(image im, int x, int y, int z) {
  gamma_decompress_value(get_pixel(im, x, y, z));
}

void gamma_compress_image(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      for (int z=0; z<im.c; z++) {
        set_pixel(im, x, y, z, get_pixel_compressed(im, x, y, z));
      }
    }
  }
}

void gamma_decompress_image(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      for (int z=0; z<im.c; z++) {
        set_pixel(im, x, y, z, get_pixel_decompressed(im, x, y, z));
      }
    }
  }
}

float* transform_ciexyz_to_srgb(float* v) {

}

void ciexyz_to_srgb(image im) {

}

float* transform_srgb_to_ciexyz(float* v) {
  float w[3];
  static float M[3][3] = {0.41239080, 0.35758434, 0.18048079, 0.21263901, 0.71516868, 0.07219232, 0.01933082, 0.11919478, 0.95053215};
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++)
      w[i] = M[i][j] * v[j];
  }
  memcpy(v, w, sizeof(float)*3);
  return v;
}

void srgb_to_ciexyz(image im) {
  float v[3];
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      for (int z=0; z<im.c; z++)
        v[z] = get_pixel(im, x, y, z);
      transform_srgb_to_ciexyz(v);
      for (int z=0; z<im.c; z++)
        set_pixel(im, x, y, z, v[z]);
    }
  }
}

void cieluv_to_ciexyz(image im) {}

void ciexyz_to_cieluv(image im) {}

void cieluv_to_hcl(image im) {}

void hcl_to_cieluv(image im) {}