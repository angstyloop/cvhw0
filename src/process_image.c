#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"
#include "matrixlib.h"

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

float* get_tsv(image im, int x, int y) {
  float* v = calloc(im.c, sizeof(float));
  for (int z=0; z<im.c; z++) {
    v[z] = get_pixel(im, x, y, z);
  }
  return v;
}

void set_tsv(image im, int x, int y, float* rgb) {
  for (int z=0; z<im.c; z++) {
    set_pixel(im, x, y, z, rgb[z]);
  }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w*im.h*im.c*sizeof(float));
    return copy;
}

// convention in video image processing
// https://en.wikipedia.org/wiki/Luma_(video)
float luma(float r, float g, float b) {
  return 0.299*r + 0.587*g + .114*b;
}

image convert_image_rgb_to_grayscale(image im)
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
// based on https://en.wikipedia.org/wiki/HSL_and_HSV
// (a bit different since all our rgb/hsv values are in [0,1])
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

void convert_image_rgb_to_hsv(image im)
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

void convert_image_hsv_to_rgb(image im)
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

float* get_tsv_compressed(image im, int x, int y) {
  float* tsv_compressed = calloc(im.c, sizeof(float));
  for (int z=0; z<im.c; ++z) {
    tsv_compressed[z] = get_pixel_compressed(im, x, y, z);
  }
  return tsv_compressed;
}

float get_pixel_decompressed(image im, int x, int y, int z) {
  return gamma_decompress_value(get_pixel(im, x, y, z));
}

float* get_tsv_decompressed(image im, int x, int y) {
  float* tsv_decompressed = calloc(im.c, sizeof(float));
  for (int z=0; z<im.c; ++z) {
    tsv_decompressed[z] = get_pixel_compressed(im, x, y, z);
  }
  return tsv_decompressed;
}

void gamma_compress_image(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* tsv = get_tsv_compressed(im, x, y);
      set_tsv(im, x, y, tsv);
      free(tsv);
    }
  }
}

void gamma_decompress_image(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* tsv = get_tsv_decompressed(im, x, y);
      set_tsv(im, x, y, tsv);
      free(tsv);
    }
  }
}

float* rgb_to_ciexyz_transform_matrix(float* r, float* g, float* b, float* w) {
  float A[9];
  A[0] = r[0];
  A[3] = r[1];
  A[6] = 1 - r[0] - r[1];
  A[1] = g[0];
  A[4] = g[1];
  A[7] = 1- g[0] - g[1];
  A[2] = b[0];
  A[5] = b[1];
  A[8] = 1 - b[0] - b[1];
  float* A_inv = inverse(A, 3);
  float* S = multiply(A_inv, w, 3, 3);
  free(A_inv);
  float* M = malloc(9*sizeof(float));
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      M[i*3+j] = S[(i*3+j)%3] * A[i*3+j];
    }
  }
  free(S);
  return M;
}

// v = rgb tristimulus vector
// r, g, b = the (2-d) xy chromacity coordinates of the nominal primaries: red, green, blue.
// w = the (3-d) 1931 CIE color space coordinates of the nominal white point.

// v = xyz tristimulus vector (CIE 1931 XYZ color space)
float* ciexyz_to_rgb(float* v, float* r, float* g, float* b, float* w) {
  float* M = rgb_to_ciexyz_transform_matrix(r, g, b, w);
  float* M_inv = inverse(M, 3);
  free(M);
  float* u = multiply(M_inv, v, 3, 3);
  free(M_inv);
  return u;
}

float* rgb_to_ciexyz(float* v, float* r, float* g, float* b, float* w) {
  float* M = rgb_to_ciexyz_transform_matrix(r, g, b, w);
  float* u = multiply(M, v, 3, 3);
  free(M);
  return u;
}

void convert_image_ciexyz_to_rgb(image im, float* r, float* g, float* b, float* w) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* v = get_tsv(im, x, y);
      float* u = ciexyz_to_rgb(v, r, g, b, w);
      free(v);
      set_tsv(im, x, y, u);
      free(u);
    }
  }
}

void convert_image_rgb_to_ciexyz(image im, float* r, float* g, float* b, float* w) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* v = get_tsv(im, x, y) ;
      float* u = rgb_to_ciexyz(v, r, g, b, w);
      free(v);
      set_tsv(im, x, y, u);
      free(u);
    }
  }
}

float* ciexyz_to_uv_chromacity(float xyz[3]) {
  float* uv = malloc(2 * sizeof(float));
  uv[0] = 4*xyz[0]/(xyz[0]+15*xyz[1]+3*xyz[2]);
  uv[1] = 9*xyz[1]/(xyz[0]+15*xyz[1]+3*xyz[2]);
  return uv;
}

// p = cieluv coordinate vector (Luv), w = nominal white point, in XYZ.
// in place, returns pointer
float* cieluv_to_ciexyz(float luv[3], float w[3]) {
  float* xyz = calloc(3, sizeof(float));
  float* uv = ciexyz_to_uv_chromacity(w);
  float u_ = luv[1]/(13*luv[0])+uv[0];
  float v_ = luv[2]/(13*luv[0])+uv[1];
  free(uv);
  xyz[1] = luv[0] > 8 ? w[1]*pow((luv[0]+16)/116.,3) : w[1]*luv[0]*pow(13/29.,3);
  xyz[0] = xyz[1]*9*u_/(4*v_);
  xyz[2] = xyz[1]*(12-3*u_-20*v_)/(4*v_);
  return xyz;
}

// transform an image from CIELUV color space to CIEXYZ color space, with respect to a nominal white point w in CIEXYZ color space
void convert_image_cieluv_to_ciexyz(image im, float w[3]) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++){
      float* tsv = get_tsv(im, x, y);
      tsv = cieluv_to_ciexyz(tsv, w);
      set_tsv(im, x, y, tsv);
      free(tsv);
    }
  }
}

// p = ciexyz coordinate vector, w = nominal white point (in ciexyz)
float* ciexyz_to_cieluv(float xyz[3], float w[3]) {
  float* luv = calloc(3, sizeof(float));
  // set uv to the (u_,v_) chromacity of the given CIEXYZ point p
  float* uv = ciexyz_to_uv_chromacity(xyz);
  // save these and reuse uv;
  float u_ = uv[0];
  float v_ = uv[1];
  // set uv to the (u_,v_) chromacity of the CIEXYZ whitepoint w
  ciexyz_to_uv_chromacity(w);
  luv[0] = xyz[1]/w[1] > pow(6/29.,3) ? 116*pow(xyz[1]/w[1],1/3.)-16 : pow(29/3.,3)*xyz[1]/w[1];
  luv[1] = 13*luv[0]*(u_-uv[0]);
  luv[2] = 13*luv[0]*(v_-uv[1]);
  free(uv);
  return luv;
}

void convert_image_ciexyz_to_cieluv(image im, float w[3]) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* tsv = get_tsv(im, x, y);
      tsv = ciexyz_to_cieluv(tsv, w);
      set_tsv(im, x, y, tsv);
      free(tsv);
    }
  }
}

float* cieluv_to_hcl(float* luv) {
  float* hcl = calloc(3, sizeof(float));
  hcl[0] = atan2(luv[2], luv[1]);
  hcl[1] = sqrt(pow(luv[1], 2) + pow(luv[2], 2));
  hcl[2] = luv[0];
  return hcl;
}

//keep going
void convert_image_cieluv_to_hcl(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* tsv_in = get_tsv(im, x, y);
      float* tsv_out = cieluv_to_hcl(tsv_in);
      free(tsv_in);
      set_tsv(im, x, y, tsv_out);
      free(tsv_out);
    }
  }
}

float* hcl_to_cieluv(float* hcl) {
  float* luv = calloc(3, sizeof(float));
  luv[0] = hcl[2];
  luv[1] = hcl[1] * cos(hcl[0]);
  luv[2] = hcl[1] * sin(hcl[0]);
  return luv;
}

void convert_image_hcl_to_cieluv(image im) {
  for (int x=0; x<im.w; x++) {
    for (int y=0; y<im.h; y++) {
      float* tsv_in = get_tsv(im, x, y);
      float* tsv_out = hcl_to_cieluv(tsv_in);
      free(tsv_in);
      set_tsv(im, x, y, tsv_out);
      free(tsv_out);
    }
  }
}