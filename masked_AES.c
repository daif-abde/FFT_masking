#include "masked_AES.h"


// =============================tests =============================

void initAESParam(param *p){
  mat s0, x;
  int i, j, l;

  initMat(&s0, p->m, p->t);
  for(j = 0; j < p->m; j++){
    for(i = 0; i < p->t; i++){
      s0.m[j][i] = 0x63;
    }
  }
  mask(p, &S0, &s0);
  initMat(&x, p->m, p->t);
  for(i = 0; i < 11; i++){
    for(j = 0; j < p->m; j++){
      for(l = 0; l < p->t; l++){
        x.m[j][l] = key[i][j*(p->t) + l];
      }
    }
    mask(p, &maskedKey[i], &x);
    
  }
}

// calculate mask(x ^ q) where q in {2, 4, 16}
void power(param *p, mat *y, mat *z, int q)
{
  mat u;
  int j;
  copyMat(&u, z);
  switch (q)
  {
  case 16:
    mapMultMatIn(&u, &u);    // u <- z^2
    mapMultMatIn(&u, &u);    // u <- z^4
    mapMultMatIn(&u, &u);    // u <- z^8
    mapMultMatIn(&u, &u);    // u <- z^16
    j = 2;
    break;
  case 4 :
    mapMultMatIn(&u, &u);    // u <- z^2
    mapMultMatIn(&u, &u);    // u <- z^4
    j = 1;
    break;
  case 2 :
    mapMultMatIn(&u, &u);    // u <- z^2
    j = 0;
    break;
  default:
    printf("Error q must be equal to : 2 or 4 or 16\n");
    exit(-1);
  }
  multMat(y, &u, p->P[j]);
  
}

void power254(param *p, mat *y, mat *z)
{
    mat w;
    mat v;
    mat u;
    mat h;
    power(p, &v,  z,  2);
    sMult(p, &u,  z, &v);
    power(p, &w, &u,  4);
    sMult(p,  y, &u, &w);
    power(p, &h,  y, 16);
    sMult(p, &u, &h, &w);
    sMult(p,  y, &u, &v);
}


// z = mask(SybBytes(x)) where mask(x) = z
void subBytes(param *p, mat *z)
{
    mat u;
    mat v;
    mat y;

    setZero(&y, p->m, p->n);
    addInMat(&y, &S0); // y = mask(63, ..., 63)
    power254(p, &u, z);   // u = z^-1
    mapMultMatByte(&v, &u, 0x05); // v = 05 . z^-1
    addInMat(&y, &v);

    power(p, &v, &u, 2);   // v = z^-2
    mapMultMatByte(&u, &v, 0x09); // u = 09 . z^-2
    addInMat(&y, &u);

    power(p, &u, &v, 2);   // u = z^-4
    mapMultMatByte(&v, &u, 0xF9); // v = F9. z^-4
    addInMat(&y, &v);

    power(p, &v, &u, 2);   // v = z^-8
    mapMultMatByte(&u, &v, 0x25); // u = 25 . z^-8
    addInMat(&y, &u);

    power(p, &u, &v, 2);   // u = z^-16
    mapMultMatByte(&v, &u, 0xF4); // v = F4. z^-16
    addInMat(&y, &v);

    power(p, &v, &u, 2); // v = z^-32
    addInMat(&y, &v);

    power(p, &u, &v, 2);   // u = z^-64
    mapMultMatByte(&v, &u, 0xB5); // v = B5. z^-64
    addInMat(&y, &v);

    power(p, &v, &u, 2);   // v = z^-128
    mapMultMatByte(&u, &v, 0x8F); // u = 8F . z^-128
    addInMat(&y, &u);
    copyMat (z, &y);
}

void mixColumn(param *p, mat *z){
    mat y;
    multMat(&y, z, p->LN);
    copyMat(z, &y);
}

// return mask(0, 0, x_pos, 0, 0)
void extract_position(param *p, mat *z, int pos){
    mat y;
    multMat(&y, z, p->I[pos]);
    copyMat(z, &y);
}

void switchPosition(param *p, mat *res, mat *z, int ij, int dj, int shift){
  mat I, Ni, N;
  setZero(&I, p->n, p->n);
  I.m[ij][dj] = 0x01;
  multMat(&Ni, p->Ni, &I);
  multMat(&N, &Ni, p->N);
  multMat(res, z, &N);

  shiftRowsMat(res, shift);


}

void shiftRows(param *p, mat *z)
{
  mat zz, tmp1, tmp2;
  int j;
  setZero(&zz, p->m, p->n);
  initMat(&tmp1, p->m, p->n);
  initMat(&tmp2, p->m, p->n);
  switch (p->t)
  {
  case 4:
    multMat(&tmp1, z, p->I[0]);                         addInMat(&zz, &tmp1);
    multMat(&tmp1, z, p->I[1]); shiftRowsMat(&tmp1, 1); addInMat(&zz, &tmp1);
    multMat(&tmp1, z, p->I[2]); shiftRowsMat(&tmp1, 2); addInMat(&zz, &tmp1);
    multMat(&tmp1, z, p->I[3]); shiftRowsMat(&tmp1, 3); addInMat(&zz, &tmp1);
    copyMat(z, &zz);
    break;
  case 8:
    multMat(&zz, z, p->Shift);
    multMat(&tmp1, &zz, p->Permute1);
    multMat(&tmp2, &zz, p->Permute2);
    swapRows(&tmp1, 0, 1);
    addMat(z, &tmp1, &tmp2);
    break;
  case 16:
    multMat(&zz, z, p->Shift);
    copyMat(z, &zz);
    break;
  }
}

void aes_encrypt(param *p, mat *z){
    addInMat(z, &maskedKey[0]);
    for(int i = 1; i < 10; i++){
        subBytes(p, z);
        shiftRows(p, z);
        mixColumn(p, z);
        addInMat(z, &maskedKey[i]);
    }
    subBytes(p, z);
    shiftRows(p, z);
    addInMat(z, &maskedKey[10]);
}

void test_mixColumn(param *p){
  mat x, z;
  fillMat(&x, p->m, p->t);
  printMat("x", &x);
  mask(p, &z, &x);
  // printMat("z", &z);

  mixColumn(p, &z);
  testCalcul(p, "MixColumn(x)", &z);

}

void test_mask(param *p)
{
  mat in, res, out;
  fillMat(&in , p->m, p->t);
  printMat("in", &in);

  mask(p, &out, &in);
  testCalcul(p, "in", &out);

}


void test_mult(param *p){
  mat z1, z2, z, out;
  mat x1, x2, x;

  setZero(&x1 , p->m, p->t);
  setZero(&x2 , p->m, p->t);
  // printMat("x1", &x1);
  // printMat("x2", &x2);

  mask(p, &z1, &x1);
  mask(p, &z2, &x2);
  testCalcul(p, "x1", &z1);
  testCalcul(p, "x2", &z2);

  sMult(p, &z, &z1, &z2);
  testCalcul(p, "x  ", &z);
  mapMultMat(&x, &x1, &x2);
  printMat("x  ", &x);
  printf("----------------------------\n");
}


void test_power(param *p){
  mat x, z, y, u;
  fillMat(&x, p->m, p->t);
  // printMat("x", &x);
  mask(p, &z, &x);

  power254(p, &y, &z);
  // power(p, &y,  &z,  16);
  testCalcul(p, "y", &y);
  mapPowerMat(&u, &x, 254);
  printMat("u", &u);
}

void test_subBytes(param *p){
  mat x, z;
  fillMat(&x, p->m, p->t);
  printMat("x", &x);
  mask(p, &z, &x);
  // printMat("z", &z);

  subBytes(p, &z);
  testCalcul(p, "S(x)", &z);
}

void test_shiftRows(param *p){
  mat x, z;
  fillMat(&x, p->m, p->t);
  printMat("x", &x);
  mask(p, &z, &x);

  shiftRows(p, &z);
  testCalcul(p, "Shift(x)", &z);

}

void test_switchPosition(param *p){
  mat x, z, zz;
  fillMat(&x, p->m, p->t);
  printMat("x", &x);
  mask(p, &z, &x);
  initMat(&zz, p->m, p->n);
  switchPosition(p, &zz, &z, 1, 5, 0);
  
  testCalcul(p, "switch position(x)", &zz);
}

void test_enc(param *p){
  int i;
  mat x, z;
  // fillMat(&x, p->m, p->t);
  setZero(&x, p->m, p->t);
  // printMat("x", &x);
  mask(p, &z, &x);
  
  aes_encrypt(p, &z);
  testCalcul(p, "Enc(x)", &z);

  // x
  // 00 01 02 03 
  // 00 01 02 03 
  // 00 01 02 03 
  // 00 01 02 03 
  // Enc(x)
  // BB E6 A2 34 
  // 04 B6 95 50 
  // 8A F9 FE AF 
  // 03 E8 F4 6C 


  // 000000 --> 66e94bd4ef8a2c3b884cfa59ca342b2e
}

void test(int d, int t){
  param p;
  // printf("------------ t = %d -----------\n", t);
  initParam(&p, d, t);
  // printParam(&p);
  initAESParam(&p);
  test_mask(&p);
  // test_mult(&p);
  // testCalcul(&p, "key", &maskedKey[1]);
  // test_power(&p);
  // test_subBytes(&p);
  // test_mixColumn(&p);
  // test_shiftRows(&p);
  // test_enc(&p);
  freeParam(&p);
}


void test_perfs(){
  param p;
  int d = 25;
  int t = 4;
  int i, j, k;
  double elapsed;
  clock_t start, end;
  mat x, z;
  int tt[3] = {4, 8, 16};
  int dd[2] = {7, 25};
  for(k = 0; k < 2; k++){
    d = dd[k];
    // if(d != 7) continue;
    for(j = 0; j < 3; j++){
      t = tt[j];
      if(t > d) break;
      
      initParam(&p, d, t);
      initAESParam(&p);

      fillMat(&x, p.m, t);
      mask(&p, &z, &x);
      
      start = clock();
      for(i = 0; i < 50; i++)
        aes_encrypt(&p, &z);
      end = clock();
      elapsed = ((double)end - start) / CLOCKS_PER_SEC;
      printf("50xAES       \t :\t t = %d\t d = %d\t n = %d\t\t Elapsed time is: %.2f s.\n", t, d, p.n, elapsed);


      start = clock();
      for(i = 0; i < 500; i++)
        subBytes(&p, &z);
      end = clock();
      elapsed = ((double)end - start) / CLOCKS_PER_SEC;
      printf("500xSubBytes \t :\t t = %d\t d = %d\t n = %d\t\t Elapsed time is: %.2f s.\n", t, d, p.n, elapsed);

      start = clock();
      for(i = 0; i < 500; i++)
        shiftRows(&p, &z);
      end = clock();
      elapsed = ((double)end - start) / CLOCKS_PER_SEC;
      printf("500xShiftRows \t :\t t = %d\t d = %d\t n = %d\t\t Elapsed time is: %.2f s.\n", t, d, p.n, elapsed);

      start = clock();
      for(i = 0; i < 500; i++)
        mixColumn(&p, &z);
      end = clock();
      elapsed = ((double)end - start) / CLOCKS_PER_SEC;
      printf("500xMixColumn \t :\t t = %d\t d = %d\t n = %d\t\t Elapsed time is: %.2f s.\n", t, d, p.n, elapsed);

      start = clock();
      for(i = 0; i < 500; i++)
        addInMat(&z, &maskedKey[i%11]);
      end = clock();
      elapsed = ((double)end - start) / CLOCKS_PER_SEC;
      printf("500xAdd \t :\t t = %d\t d = %d\t n = %d\t\t Elapsed time is: %.2f s.\n", t, d, p.n, elapsed);
      printf("-------------------------\n");
      freeParam(&p);
    }
  }

}


int main()
{
  // srand(time(NULL));
  // test_perfs();

  test(7, 4);
  test(25, 4);
  test(25, 8);
  test(25, 16);

  return 0;
}