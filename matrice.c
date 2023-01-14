#include "matrice.h"

// ============================= F_256 ===========================
byte xtime(byte a)
{
    if ((a >> 7) & 1)
        return (a << 1) ^ 0x1B; // 0x1B = x^4+x^3+x+1
    return (byte)(a << 1);
}

byte xpower(byte a, int q)
{
    int i;
    for (i = 0; i < q; i++)
        a = xtime(a);
    return a;
}

byte multByte(byte a, byte b)
{
    int i;
    byte sum = 0;
    for (i = 0; i < 8; i++)
        if ((b >> i) & 1)
            sum ^= xpower(a, i);
    return sum;
}

byte invByte(byte a)
{
    int i;
    byte y = a;
    if ((a == 0) || (a == 1))
        return a;
    for (i = 0; i < 6; i++)
        y = multByte(multByte(y, y), a);
    return multByte(y, y);
}

byte expo(byte a, int q)
{
    byte p = (byte)0x01;
    for (int i = 0; i < q; i++)
        p = multByte(p, a);
    return p;
}

// =====================================================================

void initMat(mat *m, int l, int c){
  m->l = l;
  m->c = c;
}

void setZero(mat*m, int l, int c){
  int i;
  m->l = l;
  m->c = c;
  for(i = 0; i < l; i++){
    memcpy(m->m[i], zero, c);
  }
}

void freeMat(mat *m){
  m->l = 0;
  m->c = 0;
}

void fillMat(mat *m, int l, int c){
  int i, j;
  initMat(m, l, c);
  int ind = 0;
  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      // m->m[i][j] = (byte) (rand()%0x100);
      m->m[i][j] = alea[ind%64]; ind ++;
    }
  }
}

void printMat(char* name, mat *m){
  int i, j;
  printf("%s\t\n", name);
  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      printf("%02X ", m->m[i][j]);
    }
    printf("\n");
  }
}

void copyExistingMat(mat *cp, mat *original){
  int i;
  cp->c = original->c;
  cp->l = original->l;
  for(i = 0; i < original->l; i++){
    memcpy(cp->m[i], original->m[i], original->c);
  }
}

void copyMat(mat *cp, mat *original){
  int i;
  cp->l = original->l;
  cp->c = original->c;
  for(i = 0; i < original->l; i++){
    memcpy(cp->m[i], original->m[i], original->c);
  }
}

void copySubMat(mat *cp, int idebut, int jdebut, mat *original, int iodebut, int jodebut, int isize, int jsize){
  int i;
  for(i = 0; i < isize; i++){
    memcpy(cp->m[idebut+i]+jdebut, original->m[iodebut+i]+jodebut, jsize);
  }
}

void transposeMat(mat *m){
  int i, j;
  byte tmp;
  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      tmp = m->m[i][j];
      m->m[i][j] = m->m[j][i];
      m->m[j][i] = tmp;
    }
  }
}

void identityMat(mat *m, int l){
  int i;
  m->l = l;
  m->c = l;
  for(i = 0; i < l; i++){
    memcpy(m->m[i], zero, l);
    m->m[i][i] = 0x01;
  }
}

void zeroMat(mat *m, int l, int c){
  int i;
  m->l = l;
  m->c = c;
  for(i = 0; i < l; i++){
    memcpy(m->m[i], zero, c);
  }
}

void concatCols(mat *res, mat *m1, mat *m2){
  int i;
  if((m1->l) != (m2->l)){
    printf("Erreur concatCols\n");
    exit(-1);
  }
  res->l = m1->l;
  res->c = (m1->c) + (m2->c);
  for(i = 0; i < m1->l; i++){
    memcpy(res->m[i], m1->m[i], m1->c);
    memcpy(res->m[i]+(m1->c), m2->m[i], m2->c);
  }
}

void concatRows(mat *res, mat *m1, mat *m2){
  int i;
  if((m1->c) != (m2->c)){
    printf("Erreur concatRows\n");
    exit(-1);
  }
  res->l = (m1->l) + (m2->l);
  res->c = m1->c;
  for(i = 0; i < m1->l; i++){
    memcpy(res->m[i], m1->m[i], m1->c);
  }
  for(i = 0; i < m2->l; i++){
    memcpy(res->m[i+(m1->l)], m2->m[i], m2->c);
  }
}

void swapRows(mat *m, int i1, int i2){
  byte tmp[SIZE_L];
  memcpy(tmp, m->m[i1], m->c);
  memcpy(m->m[i1], m->m[i2], m->c);
  memcpy(m->m[i2], tmp, m->c);
}

void shiftRowsMat(mat *m, int shift){
  mat cp;
  int i;
  if(!shift) return;
  copyMat(&cp, m);
  for(i = 0; i < cp.l; i++){
    memcpy(m->m[i], cp.m[(cp.l + i + shift)%(cp.l)], cp.c);
  }
}

void subMat(mat *res, mat *m, int idebut, int isize, int jdebut, int jsize){
  int i;
  res->l = isize;
  res->c = jsize;
  for(i = 0; i < isize; i++){
    memcpy(res->m[i], m->m[i+idebut]+jdebut, jsize);
  }
}

void addMat(mat *res, mat *m1, mat *m2){
  int i, j;
  if(((m1->c) != (m2->c)) || ((m1->l) != (m2->l))){
    printf("Erreur addMat\n");
    exit(-1);
  }
  res->l = m1->l;
  res->c = m1->c;
  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m1->c; j++){
      res->m[i][j] = (m1->m[i][j]) ^ (m2->m[i][j]);
    }
  }
}

// m1 = m1 + m2
void addInMat(mat *m1, mat *m2){
  int i, j;
  if(((m1->c) != (m2->c)) || ((m1->l) != (m2->l))){
    printf("Erreur addInMat\n");
    exit(-1);
  }
  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m1->c; j++){
      m1->m[i][j] = (m1->m[i][j]) ^ (m2->m[i][j]);
    }
  }
}

void multMat(mat *res, mat *m1, mat *m2){
  int i, j, k;
  if(((m1->c) != (m2->l))){
    printf("Erreur multMat\n");
    exit(-1);
  }
  setZero(res, m1->l, m2->c);
  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m2->c; j++){
      for(k = 0; k < m1->c; k++){
        (res->m[i][j]) ^= multByte(m1->m[i][k],  m2->m[k][j]);
      }
    }
  }
}

void mapPowerMat(mat *res, mat *m, int q){
  int i, j;
  res->l = m->l;
  res->c = m->c;
  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      res->m[i][j] = expo(m->m[i][j], q);
    }
  }
}

void mapMultMatByte(mat *res, mat *m, byte x){
  int i, j;
  res->l = m->l;
  res->c = m->c;

  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      res->m[i][j] = multByte(m->m[i][j], x);
    }
  }
}

void multTransposeMat(mat *res, mat *m1, mat *m2){
  int i, j, k;
  if(((m1->c) != (m2->c))){
    printf("Erreur multTransposeMat\n");
    exit(-1);
  }
  setZero(res, m1->l, m1->c);

  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m2->l; j++){
      for(k = 0; k < m1->c; k++){
        (res->m[i][j]) ^= multByte(m1->m[i][k],  m2->m[j][k]);
      }
    }
  }
}

// multiplication terme à terme dans res (dimesions de m1 = dimensions de m2)
void mapMultMat(mat *res, mat *m1, mat* m2){
  int i, j;
  if(((m1->c) != (m2->c)) || ((m1->l) != (m2->l))){
    printf("Erreur mapMultMat\n");
    exit(-1);
  }
  res->l = m1->l;
  res->c = m1->c;
  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m1->c; j++){
      res->m[i][j] = multByte(m1->m[i][j], m2->m[i][j]);
    }
  }
}

// m1 = mapMultMat(m1, m2)
void mapMultMatIn(mat *m1, mat *m2){
  int i, j;
  if(((m1->c) != (m2->c)) || ((m1->l) != (m2->l))){
    printf("Erreur mapMultMatIn\n");
    exit(-1);
  }
  for(i = 0; i < m1->l; i++){
    for(j = 0; j < m1->c; j++){
      m1->m[i][j] = multByte(m1->m[i][j], m2->m[i][j]);
    }
  }
}

void inverseMat(mat *inv, mat *m){
  int i, j, jj, k, r = 0;
  byte x;
  mat n;
  if(m->c != m->l){
    printf("Erreur inverseMat\n");
    exit(-1);
  }
  identityMat(inv, m->c);
  copyMat(&n, m);

  for(j = 0; j < n.c; j++){
    k = r;
    for(i = r; i < n.l; i++){
      if(n.m[i][j] > n.m[k][j])
        k = i;
    }
    if(n.m[k][j]){
      x = invByte(n.m[k][j]);
      for(i = 0; i < m->c; i++){
        n.m[k][i] = multByte(n.m[k][i], x);
        inv->m[k][i]  = multByte(inv->m[k][i], x);
      }
      if(k != r){
        swapRows(&n, k, r);
        swapRows(inv, k, r);
      }
      for(i = 0; i < n.l; i++){
        if( i!= r){
          x = n.m[i][j];
          for(jj = 0; jj < n.c; jj++){
            n.m[i][jj] = n.m[i][jj] ^ multByte(n.m[r][jj], x);
            inv->m[i][jj]  = inv->m[i][jj] ^  multByte(inv->m[r][jj], x);
          }
        }
      }
      r++;
    }
  }
}

// m x inv = id : l <= c
void pseudoInverseRightMat(mat *inv, mat *m){
  mat tmp1;
  mat tmp2;
  if(m->l > m->c){
    printf("Erreur pseudoInverseRightMat\n");
    exit(-1);
  }

  //return m^T . (m.m^T)^-1
  multTransposeMat(&tmp1, m, m);
  inverseMat(&tmp2, &tmp1);
  copyMat(&tmp1, m);
  transposeMat(&tmp1);
  multMat(inv, &tmp1, &tmp2);

}

// inv x m = id : c <= l
void pseudoInverseLeftMat(mat *inv, mat *m){
  mat tmp1;
  mat tmp2;

  if(m->c > m->l){
    printf("Erreur pseudoInverseLeftMat\n");
    exit(-1);
  }
  //return (M^T.M)^-1 . M^T
  copyMat(&tmp1, m);
  transposeMat(&tmp1);
  multMat(&tmp2, &tmp1, m);
  inverseMat(&tmp1, &tmp2);
  multTransposeMat(inv, &tmp1, m);

}

// m[i,j] = m[i,j]^q
void powerMat(mat *m, int q){
  int i, j;
  for(i = 0; i < m->l; i++){
    for(j = 0; j < m->c; j++){
      m->m[i][j] = expo(m->m[i][j], q);
    }
  }
}