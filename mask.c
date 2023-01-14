#include "mask.h"

// -------------------- Parameters --------------------
void maskMatriceProd(mat *res, mat *L, mat *Ni, mat *N, int n, int t){
  // LN = Ni * L.augment(matrix(F, t, n-t)).stack(matrix(F, n-t, t).augment(matrix(F, n-t, n-t, 1))) * N
  // Ai.augment(matrix(F, d+1, d)).stack(matrix(F, d, d+1).augment(matrix(F, d, d, 1)))*M
  mat tmp1, tmp2, tmp3, tmp4;
  setZero(&tmp1, t, n-t);
  concatCols(&tmp2, L, &tmp1);
  setZero(&tmp1, n-t, t);
  identityMat(&tmp3, n-t);
  concatCols(&tmp4, &tmp1, &tmp3);
  concatRows(&tmp1, &tmp2, &tmp4);
  multMat(&tmp2, Ni, &tmp1);
  multMat(res, &tmp2, N);
}

void initParam(param *p, int d, int t){
  int i, j, k;
  mat tmp1, tmp2, Uj, U1;
  mat Ai;
  int n = 2*d+1;
  if((t != 4) && (t != 8) && (t != 16)){
    printf("Le paramètre t est incorrecte : Les valeurs possibles sont {4, 8, 16}\n");
    exit(EXIT_FAILURE);
  }
  p->d = d;
  p->t = t;
  p->n = n;
  p->m = 16/t;
  p->alpha = 0x13; // Primitf element
  p->omega = 0x0D; // alpha^(255/n)
  if(n == 51)
    p->omega = 0x4D; // alpha^(255/n)
  p->M  = malloc(sizeof(mat));  // n x n
  p->Mi = malloc(sizeof(mat));  // n x n
  p->A  = malloc(sizeof(mat));  // (d+1) x (d+1)
  p->A1 = malloc(sizeof(mat));  // t x t
  p->A2 = malloc(sizeof(mat));  // (d+1-t) x t
  p->A3 = malloc(sizeof(mat));  // (d+1-t) x t
  p->A1i= malloc(sizeof(mat));  // t x t 
  p->U  = malloc(sizeof(mat));  // 1 x n
  p->G  = malloc(sizeof(mat));  // (t-1) x n 
  p->N  = malloc(sizeof(mat));  // n x n
  p->Ni = malloc(sizeof(mat));  // n x n
  p->LN = malloc(sizeof(mat));  // n x n
  initMat(p->M, n, n);
  initMat(p->Mi, n, n);  
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      (p->M)->m[i][j] = expo(p->omega, i*j);
      (p->Mi)->m[i][j] = expo(p->omega, (i*j*255-(i*j))%255);
    }
  }
  // ------------------------- calcul A ------------------
  initMat(p->A, d+1, d+1);
  for(i = 0; i < d+1; i++){
    for(j = 0; j < d+1; j++){
      p->A->m[i][j] = expo(p->alpha, i*j);
    }
  }

  inverseMat(&Ai, p->A);
  
  subMat(p->A1, p->A, 0, t, 0, t);
  subMat(p->A2, p->A, t, d+1-t, 0, t);

  inverseMat(p->A1i, p->A1);

  multMat(p->A3, p->A2, p->A1i);

  // ------------------------- calcul U ------------------
  setZero(&tmp2, t, t);
  for(i = 0; i < t; i++){
    tmp2.m[i][i] = expo(p->alpha, i*d);
  }
  multMat(&Uj, &tmp2, p->A1i);
  setZero(&U1, 1, n);
  for(j = 0; j < t; j++){
    for(i = 0; i < t; i++){
      U1.m[0][j] ^= Uj.m[i][j];
    }
  }

  multMat(p->U, &U1, p->M);
  
  // ------------------------- calcul G ------------------

	// # G_i(x) = sum_{j=1}^t U_j(x) . alpha^(jd-jt+j+ij)
	// G1 = matrix(F, t-1, t)
  setZero(&tmp1, t-1, n);

  for(i = 1; i < t; i++){
    for(j = 0; j < t; j++){
      for(k = 0; k < t; k++){
        tmp1.m[i-1][k] ^= multByte(p->A->m[j][d-t+1+i], Uj.m[j][k]);
      }
    }
  }

  multMat(p->G, &tmp1, p->M);

  // ------------------------- calcul N ------------------

  // N = Ai.augment(matrix(F, d+1, d)).stack(matrix(F, d, d+1).augment(matrix(F, d, d, 1)))*M
  // mask(x, r, 0) = (x, r, 0) . N
  identityMat(&tmp1, n);
  maskMatriceProd(p->N, &Ai, &tmp1, p->M, n, d+1);
  inverseMat(p->Ni, p->N);

  // ------------------------- calcul LN ------------------

  initMat(&tmp1, p->t, p->t);
  for(i = 0; i < ((p->t)/4); i++){
    tmp1.m[i*4+0][i*4] = 0x02; tmp1.m[i*4+0][i*4+1] = 0x01; tmp1.m[i*4+0][i*4+2] = 0x01; tmp1.m[i*4+0][i*4+3] = 0x03;
    tmp1.m[i*4+1][i*4] = 0x03; tmp1.m[i*4+1][i*4+1] = 0x02; tmp1.m[i*4+1][i*4+2] = 0x01; tmp1.m[i*4+1][i*4+3] = 0x01;
    tmp1.m[i*4+2][i*4] = 0x01; tmp1.m[i*4+2][i*4+1] = 0x03; tmp1.m[i*4+2][i*4+2] = 0x02; tmp1.m[i*4+2][i*4+3] = 0x01;
    tmp1.m[i*4+3][i*4] = 0x01; tmp1.m[i*4+3][i*4+1] = 0x01; tmp1.m[i*4+3][i*4+2] = 0x03; tmp1.m[i*4+3][i*4+3] = 0x02;
  }
	// LN = Ni * L.augment(matrix(F, t, n-t)).stack(matrix(F, n-t, t).augment(matrix(F, n-t, n-t, 1))) * N
  maskMatriceProd(p->LN, &tmp1, p->Ni, p->N, n, t);

  
  // ------------------------- calcul P ------------------
  p->P = (mat **) malloc(3 * sizeof(mat));

  int power[3] = {2, 4, 16};
  for(i = 0; i < 3; i++){
    p->P[i] = (mat*) malloc(sizeof(mat));
    mapPowerMat(&tmp1, p->Ni, power[i]);
    multMat(p->P[i], &tmp1, p->N);
  }
  // ------------------------- calcul pour shiftRows ------------------

  if(t == 4){
    p->I = (mat **) malloc(16 * sizeof(mat));
    setZero(&tmp1, t, t);
    for(i = 0; i < t; i++){
      p->I[i] = (mat*) malloc(sizeof(mat));
      tmp1.m[i][i] = 0x01;
      maskMatriceProd(p->I[i], &tmp1, p->Ni, p->N, p->n, t);
      tmp1.m[i][i] = 0x00;
    }
  }
  if(t == 16){
    p->Shift = malloc(sizeof(mat));
    setZero(&tmp1, t, t);
    tmp1.m[ 0][ 0] = 1; tmp1.m[ 5][ 1] = 1; tmp1.m[10][ 2] = 1; tmp1.m[15][ 3] = 1;
    tmp1.m[ 4][ 4] = 1; tmp1.m[ 9][ 5] = 1; tmp1.m[14][ 6] = 1; tmp1.m[ 3][ 7] = 1;
    tmp1.m[ 8][ 8] = 1; tmp1.m[13][ 9] = 1; tmp1.m[ 2][10] = 1; tmp1.m[ 7][11] = 1;
    tmp1.m[12][12] = 1; tmp1.m[ 1][13] = 1; tmp1.m[ 6][14] = 1; tmp1.m[11][15] = 1;
    maskMatriceProd(p->Shift, &tmp1, p->Ni, p->N, n, t);  
  }
  if(t == 8){
    p->Shift = malloc(sizeof(mat));
    p->Permute1 = malloc(sizeof(mat));
    p->Permute2 = malloc(sizeof(mat));
    setZero(&tmp1, t, t);
    tmp1.m[ 0][ 0] = 1; tmp1.m[ 5][ 1] = 1; tmp1.m[ 2][ 2] = 1; tmp1.m[ 7][ 3] = 1;
    tmp1.m[ 4][ 4] = 1; tmp1.m[ 1][ 5] = 1; tmp1.m[ 6][ 6] = 1; tmp1.m[ 3][ 7] = 1;

    setZero(&tmp2, t, t);
    tmp2.m[ 2][ 2] = 1; tmp2.m[ 3][ 3] = 1; tmp2.m[ 5][ 5] = 1; tmp2.m[ 6][ 6] = 1;
    maskMatriceProd(p->Permute1, &tmp2, p->Ni, p->N, n, t);

    setZero(&tmp2, t, t);
    tmp2.m[ 0][ 0] = 1; tmp2.m[ 1][ 1] = 1; tmp2.m[ 4][ 4] = 1; tmp2.m[ 7][ 7] = 1;
    maskMatriceProd(p->Permute2, &tmp2, p->Ni, p->N, n, t);
    maskMatriceProd(p->Shift, &tmp1, p->Ni, p->N, n, t);  
  }
}

void printParam(param *p){
  int i;
  int power[3] = {2, 4, 16};
  printf("d = %d\n", p->d);
  printf("t = %d\n", p->t);
  printf("n = %d\n", p->n);
  printMat("M" , p->M);
  printMat("Mi", p->Mi);
  // printMat("A" , p->A);
  // printMat("A1", p->A1);
  // printMat("A2", p->A2);
  // printMat("A3", p->A3);
  // printMat("A1i",p->A1i);
  // printMat("N" , p->N);
  // printMat("Ni", p->Ni);
  // printMat("LN", p->LN);
  printMat("G" , p->G);
  printMat("U" , p->U);
  // if(p->t == 4){
  //   for(i = 0; i < p->t; i++){
  //     printf("I %d ", i); printMat("", p->I[i]);
  //   }
  // }
  // for(i = 0; i < 3; i++){
  //   printf("P %d ", power[i]); printMat("", p->P[i]);
  // }
}

void freeParam(param *p){
  int i;
  freeMat(p->M);
  freeMat(p->Mi);
  freeMat(p->A);
  freeMat(p->A1);
  freeMat(p->A2);
  freeMat(p->A3);
  freeMat(p->A1i);
  freeMat(p->N);
  freeMat(p->Ni);
  freeMat(p->LN);
  freeMat(p->G);
  freeMat(p->U);
  if(p->t == 4) for(i = 0; i < p->t; i++) freeMat(p->I[i]);
  if(p->t == 8) {
    freeMat(p->Permute1);
    freeMat(p->Permute2);
    freeMat(p->Shift);
  }
  if(p->t == 16) freeMat(p->Shift);
  for(i = 0; i < 3; i++) freeMat(p->P[i]);
}

// --------------------------FFT -----------------------------

byte mul(byte a, byte b)
{
  if (a && b)
    return (Alogtable[(Logtable[a] + Logtable[b]) % 255]);
  else
    return (0);
}

/* Primitive = α = 1+x */
/* ω = α^{17} = 1 + x^5 + x^6 + x^7 = 225 */

/* degre 14 pol modulo  x 8 + x 4 + x 2 + x + 1 =
(a 7 + a 11 + a 13 + a 14 )x 7 + (a 6 + a 10 + a 12 + a 13 )x 6 + (a 5 + a 9 + a 11 + a 12 )x 5 +
(a 4 + a 8 + a 10 + a 11 )x 4 + (a 3 + a 9 + a 10 + a 11 + a 13 + a 14 )x 3 +
(a 2 + a 8 + a 9 + a 10 + a 12 + a 13 )x 2 + (a 1 + a 8 + a 9 + a 12 + a 13 + a 14 )x+
(a 0 + a 8 + a 12 + a 14 )

   degre 15 pol modulo   x^8 + x^4 + x^2 + x =
(a 7 + a 11 + a 13 + a 14 )x 7 + (a 6 + a 10 + a 12 + a 13 + a 14 )x 6 +
(a 5 + a 9 + a 11 + a 12 + a 13 )x 5 + (a 4 + a 8 + a 10 + a 11 + a 12 )x 4 +
(a 3 + a 9 + a 10 + a 13 + a 14 )x 3 + (a 2 + a 8 + a 9 + a 12 + a 13 + a 14 )x 2 +
(a 1 + a 8 + a 12 + a 14 )x + a 0
*/

void modeigth(byte *out, byte *in)
{

  out[0] = in[0];
  out[1] = in[1] ^ in[8] ^ in[12] ^ in[14];
  out[2] = in[2] ^ in[8] ^ in[9] ^ in[12] ^ in[13] ^ in[14];
  out[3] = in[3] ^ in[9] ^ in[10] ^ in[13] ^ in[14];
  out[4] = in[4] ^ in[8] ^ in[10] ^ in[11] ^ in[12];
  out[5] = in[5] ^ in[9] ^ in[11] ^ in[12] ^ in[13];
  out[6] = in[6] ^ in[10] ^ in[12] ^ in[13] ^ in[14];
  out[7] = in[7] ^ in[11] ^ in[13] ^ in[14];

  out[8] = in[0] ^ in[8] ^ in[12] ^ in[14];
  out[9] = in[1] ^ in[8] ^ in[9] ^ in[12] ^ in[13] ^ in[14];
  out[10] = in[2] ^ in[8] ^ in[9] ^ in[10] ^ in[12] ^ in[13];
  out[11] = in[3] ^ in[9] ^ in[10] ^ in[11] ^ in[13] ^ in[14];
  out[12] = in[4] ^ in[8] ^ in[10] ^ in[11];
  out[13] = in[5] ^ in[9] ^ in[11] ^ in[12];
  out[14] = in[6] ^ in[10] ^ in[12] ^ in[13];
  out[15] = in[7] ^ in[11] ^ in[13] ^ in[14];

  return;
}

/* degre 7 pol modulo x^4+x+λ =
(a 3 + a 6 + λa 7 )x 3 + (a 2 + a 5 + λa 6 )x 2 + (a 1 + a 4 + a 7 + λa 5 )x + a 0 + λ(a 4 + a 7 ).
 */

void modfour(byte *out, byte *in, byte l)
{

  out[0] = in[0] ^ mul(l, in[4] ^ in[7]);
  out[1] = in[1] ^ in[4] ^ in[7] ^ mul(l, in[5]);
  out[2] = in[2] ^ in[5] ^ mul(l, in[6]);
  out[3] = in[3] ^ in[6] ^ mul(l, in[7]);

  return;
}

/* degre 3 pol modulo x 2 + x + λ
(a 1 + a 2 + a 3 + λa 3 )x + a 0 + λ(a 2 + a 3 )
*/

void modtwo(byte *out, byte *in, byte l)
{

  out[0] = in[0] ^ mul(l, in[2] ^ in[3]);
  out[1] = in[1] ^ in[2] ^ in[3] ^ mul(l, in[3]);

  return;
}

/* degre 1 pol modulo x + λ */

void modone(byte *out, byte *in, byte l)
{

  *out = in[0] ^ mul(l, in[1]);

  return;
}

/*
ω = 225;ω^2 = 92;ω^3 = 12;ω^4 = 224;ω^5 = 189;ω^6 = 80;ω^7 = 236;ω^8 = 93;ω^9 = 237;ω^10 = 188;ω^11 = 177;ω^12 = 176;ω^13 = 81;ω^14 = 13;ω^15 = 1;

0 , 0 , 5 , 10 , 1 , 4 , 2 , 8 , 7 , 9 , 6  , 13 , 3  , 14 , 11 , 12
0 , 1 , 2 , 3  , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15
*/

void IFFT_CT(param *p, mat *out, mat *in)
{
  byte in1[16];
  byte fout[16];
  int i;
  initMat(out, p->m, p->n);
  for(i = 0; i < p->m; i++){
    modeigth(fout, in->m[i]);

    modfour(in1     , fout    , 0);
    modfour(&in1[4] , fout    , 1);
    modfour(&in1[8] , &fout[8], 189);
    modfour(&in1[12], &fout[8], 188);

    modtwo(fout     , in1     , 0  );
    modtwo(&fout[2] , in1     , 1  );
    modtwo(&fout[4] , &in1[4] , 189);
    modtwo(&fout[6] , &in1[4] , 188);
    modtwo(&fout[8] , &in1[8] , 225);
    modtwo(&fout[10], &in1[8] , 224);
    modtwo(&fout[12], &in1[12], 92 );
    modtwo(&fout[14], &in1[12], 93 );

    modone(&(in1[0]),  fout     , 0);
    modone(&(in1[1]),  fout     , 1);
    modone(&(in1[2]),  &fout[2] , 189);
    modone(&(in1[3]),  &fout[2] , 188);
    modone(&(in1[4]),  &fout[4] , 225);
    modone(&(in1[5]),  &fout[4] , 224);
    modone(&(in1[6]),  &fout[6] , 92);
    modone(&(in1[7]),  &fout[6] , 93);
    modone(&(in1[8]),  &fout[8] , 236);
    modone(&(in1[9]),  &fout[8] , 237);
    modone(&(in1[10]), &fout[10], 80);
    modone(&(in1[11]), &fout[10], 81);
    modone(&(in1[12]), &fout[12], 12);
    modone(&(in1[13]), &fout[12], 13);
    modone(&(in1[14]), &fout[14], 177);
    modone(&(in1[15]), &fout[14], 176);

    out->m[i][0]  = in1[1];
    out->m[i][1]  = in1[4];
    out->m[i][2]  = in1[6];
    out->m[i][3]  = in1[12];
    out->m[i][4]  = in1[5];
    out->m[i][5]  = in1[2];
    out->m[i][6]  = in1[10];
    out->m[i][7]  = in1[8];
    out->m[i][8]  = in1[7];
    out->m[i][9]  = in1[9];
    out->m[i][10] = in1[3];
    out->m[i][11] = in1[14];
    out->m[i][12] = in1[15];
    out->m[i][13] = in1[11];
    out->m[i][14] = in1[13];
  }

}

// 0,225,92,12,224,189,80,236,93,237,188,177,176,81,13,1,
// 0,1,225,92,12,224,189,80,236,93,237,188,177,176,81,13,

/*   0  1   2   3   4   5   6  7  8   9  10  11  12  13  14
w^   0  1   2   3   4   5   6  7  8   9  10  11  12  13  14
     0 -14 -13 -12 -11 -10 -9 -8 -7  -6  -5  -4  -3  -2  -1
     1 -1   -2
*/

void FFT_CT(param *p, mat *out, mat *in)
{

  mat in1;
  int i, j;
  initMat(&in1, p->m, p->n);
  for(j = 0; j < p->m; j++){
    for(i = 0; i < p->n; i++){
      in1.m[j][i]  = in->m[j][((p->n)-i) %(p->n)];
    }
  }
  IFFT(p, out, &in1);
}

void FFT(param *p, mat *out, mat *in)
{
  multMat(out, in, p->M);
}

void IFFT(param *p, mat *out, mat *in)
{
  multMat(out, in, p->Mi);
}

void mask(param *p, mat *z, mat *x)
{
  mat a, a1, a2, b;

  fillMat(&a2, p->m, (p->d) + 1 - (p->t));
  multMat(&b, &a2, p->A2); // b = a2.A2
  addInMat(&b, x); // b = x + a2.A2

  // a1 = b.A1^-1  = (x + a2.A2). A1^-1
  multMat(&a1, &b, p->A1i);
  concatCols(&a, &a1, &a2);
  setZero(&a1, p->m, p->d);
  concatCols(&b, &a, &a1);  // b = (a1, a2, zero)
  FFT(p, z, &b);
}

void unmask(param *p, mat *x, mat *z)
{
  mat a, b, c;
  IFFT(p, &a, z);
  subMat(&b, &a, 0, p->m, 0, p->t);
  multMat(&c, &b, p->A1);
  subMat(&b, &a, 0, p->m, p->t, (p->d) + 1 - (p->t));
  multMat(&a, &b, p->A2);
  addMat(x, &c, &a);
}

void testCalcul(param *p, char* name, mat * z)
{
  mat res;
  unmask(p, &res, z);
  printMat(name, &res);
}

void extractLastCoefficients(param *p, mat *r, mat *y)
{
  int i, j;
  mat mu, la, nu, ifft_mu, ifft_la;

  initMat(&mu, p->m, p->n);
  initMat(&la, p->m, p->n);
  initMat(&nu, p->m, p->d);
  for(j = 0; j < p->m; j++){
    for (i = 0; i < p->d; i++)
    {
      mu.m[j][i] = mul(y->m[j][i], p->Mi->m[i][(p->d) + 1]);
      la.m[j][i] = mul(y->m[j][i + (p->d)], p->Mi->m[i + (p->d)][(p->d) + 1]);
      nu.m[j][i] = mul(y->m[j][p->d], p->Mi->m[p->d][(p->d) + 1 + i]);
    }
    mu.m[j][i] = mul(y->m[j][p->d], p->Mi->m[p->d][p->d + 1]);
    la.m[j][i] = mul(y->m[j][2 * p->d], p->Mi->m[2 * p->d][p->d + 1]);
  }

  IFFT(p, &ifft_mu, &mu);
  IFFT(p, &ifft_la, &la);
  initMat(r, p->m, p->d);

  for(j = 0; j < p->m; j++){
    for (i = 0; i < p->d; i++)
    {
      r->m[j][i] = ifft_mu.m[j][i] ^ nu.m[j][i] ^ mul(p->Mi->m[p->d][i], ifft_la.m[j][i]);
    }
  }

}

void sMult(param *p, mat *z, mat *z1, mat *z2)
{
  int i, j, k;
  mat y;
  mat c, c1, c2;
  mat fft_c1, fft_c2;
  mat jj;

  mapMultMat(&y, z1, z2);
  extractLastCoefficients(p, &c, &y);
  setZero(&c1, p->m, p->n);
  setZero(&c2, p->m, p->n);
  copySubMat(&c1, 0, 1, &c, 0, 0, p->m, (p->d) + 1 - (p->t));
  copySubMat(&c2, 0, (p->d)+1, &c, 0, 0, p->m, (p->d));
  setZero(&jj, p->m, p->n);
  for(k = 0; k < p->m; k++){
    for (i = 0; i < p->t - 1; i++)
    {
      for (j = 0; j < p->n; j++)
      {
        jj.m[k][j] ^= mul(p->G->m[i][j], c.m[k][(p->d) - (p->t) + 1 + i]);
      }
    }
  }
  FFT(p, &fft_c1, &c1);
  FFT(p, &fft_c2, &c2);

  initMat(z, p->m, p->n);
  for(k = 0; k < p->m; k++){
    for (i = 0; i < p->n; i++)
    {
      z->m[k][i] = y.m[k][i] ^ fft_c2.m[k][i] ^ mul(fft_c1.m[k][i], p->U->m[0][i]) ^ jj.m[k][i];
    }
  }
}

// z = z + maskBloc(0)
void refresh(param *p, mat *z)
{
  mat a, a1, a2, b, z0;

  fillMat(&a2, p->m, (p->d) + 1 - (p->t));
  multMat(&b, &a2, p->A2); // b = 0 + a2.A2

  // a1 = b.A1^-1  = (0 + a2.A2). A1^-1
  multMat(&a1, &b, p->A1i);
  concatCols(&a, &a1, &a2);
  setZero(&a1, p->m, p->d);
  concatCols(&b, &a, &a1);  // b = (a1, a2, zero)
  FFT(p, &z0, &b);
  addInMat(z, &z0);
}
