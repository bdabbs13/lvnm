/*****  helper.h
 *****  Header File for helper.cpp defining multiple useful functions
 *****/

#ifndef LVNM_HELPER_R
#define LVNM_HELPER_R
#define MIN_LOG 0.0000000001
#define MAX_LOG 1.0 - MIN_LOG//0.999999



void shiftFlatTable(int shift_size, int flatLength, int total, 
		    double *flatTable);

void RprintDoubleMat(int rows, int cols, double *mat);

void RprintIntMat(int rows, int cols, int *mat);

void rdirichlet(int k, double *alpha, double *x);

void getMeanVar(double *vec, int lower, int upper,
		double *Mean_t, double *Var_t, double * Len_t);

int convergenceCheck(double *logLik, int total, double qq);

void colSums(double *mat, int rows, int cols, double *totals);

void rowSums(double *mat, int rows, int cols, double *totals);

void logZeroFix(double *vec, int ll);

  

#endif

