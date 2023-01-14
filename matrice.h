#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

typedef uint8_t byte;
#define SIZE_L 800
#define SIZE_C 800

typedef struct _matrice{
  byte m[SIZE_L][SIZE_C]; // matrice
  int l;    // Lignes
  int c;    // colonnes
} mat;

byte zero[SIZE_L] = {0x00};

// gcc -c -Wall -Werror -fpic matrice.c
// gcc -shared -o matrice.so matrice.o

// ============================= F_256 ===========================
byte xtime(byte a);
byte xpower(byte a, int q);
byte multByte(byte a, byte b);
byte invByte(byte a);
byte expo(byte a, int q);

// =====================================================================

byte alea[64] = {
  0x68, 0x1F, 0x32, 0x4B, 0x69, 0x8D, 0x5B, 0x08, 0x17, 0x27, 0x4E, 0x42, 0xE1, 0x30, 0x22, 0x25,
  0x35, 0x76, 0x26, 0xB9, 0x3E, 0xC1, 0x50, 0x35, 0xB6, 0x75, 0xB3, 0x7A, 0xEA, 0x5D, 0x26, 0x7D,
  0x7C, 0xB9, 0xE9, 0x9C, 0x53, 0xAF, 0x00, 0x2C, 0xD0, 0x0B, 0x1A, 0xDD, 0xA2, 0xE8, 0xDC, 0x9E,
  0x83, 0x2C, 0xCD, 0xD3, 0x43, 0xFA, 0xCE, 0x3D, 0x5C, 0x51, 0x8F, 0x74, 0x49, 0x12, 0x13, 0x56};

void initMat(mat *m, int l, int c);

void freeMat(mat *m);

void fillMat(mat *m, int l, int c);

void printMat(char* name, mat *m);

void copyMat(mat *cp, mat *original);

void transposeMat(mat *m);

void identityMat(mat *m, int l);

void zeroMat(mat *m, int l, int c);

void concatCols(mat *res, mat *m1, mat *m2);

void concatRows(mat *res, mat *m1, mat *m2);

void swapRows(mat *m, int i1, int i2);

void subMat(mat *res, mat *m, int idebut, int isize, int jdebut, int jsize);

void addMat(mat *res, mat *m1, mat *m2);

// m1 = m1 + m2
void addInMat(mat *m1, mat *m2);

void multMat(mat *res, mat *m1, mat *m2);

void multTransposeMat(mat *res, mat *m1, mat *m2);

void inverseMat(mat *inv, mat *m);

// m x inv = id : l <= c
void pseudoInverseRightMat(mat *inv, mat *m);

// inv x m = id : c <= l
void pseudoInverseLeftMat(mat *inv, mat *m);