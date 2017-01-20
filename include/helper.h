/*****  helper.h
 *****  Header File for helper.cpp defining multiple useful functions
 *****/

#ifndef LVNM_HELPER_R
#define LVNM_HELPER_R
#include<vector>
//#define MIN_LOG 0.0000000001
#define MIN_LOG DBL_EPSILON
#define MAX_LOG 1.0 - MIN_LOG


//  Function for determining missingness.  Currently < 0 indicates missingness
bool isMissing(int ii);

void shiftFlatTable(int shift_size, int flatLength, int total,
		    double *flatTable);

void shiftFlatTable(int shift_size, int flatLength, int total,
		    int *flatTable);

void RprintDoubleMat(int rows, int cols, double *mat);

void RprintIntMat(int rows, int cols, int *mat);

void rdirichlet(int k, double *alpha, double *x);

extern "C" {
   void RGammaPrior(int *total, int *thin, int *burnin,
		    double *alpha_init, double *beta_init,
		    double *p, double *q, double *r, double *s,
		    double *prop_sd, double *alpha, double *beta);
}



void rGammaPrior(double alpha_init, double beta_init,
		 double p, double q, double r, double s,
		 double prop_sd, double *alpha, double *beta);

void rGammaPriorStep(double *alphabeta,
		     double p, double q, double r, double s,
		     double prop_sd);

void rGammaPrior(int total, int thin, int burnin,
		 double alpha_init, double beta_init,
		 double p, double q, double r, double s,
		 double prop_sd, double *alpha, double *beta);


double ldGammaPrior(double alpha, double beta,
		    double p, double q, double r, double s);

void getMeanVar(double *vec, int lower, int upper,
		double *Mean_t, double *Var_t, double * Len_t);

int convergenceCheck(double *logLik, int total, double qq);



void colSums(std::vector<std::vector<double> > const &mat,
	     std::vector<double> &totals);
void colSums(std::vector<std::vector<double> > const &mat, double *totals);
void colSums(double *mat, int rows, int cols, std::vector<double> totals);
void colSums(double *mat, int rows, int cols, double *totals);

void rowSums(double *mat, int rows, int cols, double *totals);
void rowSums(std::vector<std::vector<double> > const &mat,
	     std::vector<double> &totals);


void normalizeVec(std::vector<double> &vec);
void normalizeVec(double *vec, int ll);

void logZeroFix(std::vector<double> &vec);
void logZeroFix(double *vec, int ll);

double logCheck(double val);



#endif

