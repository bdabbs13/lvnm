//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <R.h>
#include <Rmath.h>
#include "helper.h"
//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;



//bool isMissing(int ii){
//   return (ii < 0);
//}

/******************************************/
void shiftFlatTable(int shift_size, int flatLength, int total,
		    double *flatTable){
  int flatStart = shift_size * flatLength;
  int flatEnd = total * flatLength;
  std::copy(flatTable + flatStart, flatTable + flatEnd, flatTable);
}

void shiftFlatTable(int shift_size, int flatLength, int total,
		    int *flatTable){
  int flatStart = shift_size * flatLength;
  int flatEnd = total * flatLength;
  std::copy(flatTable + flatStart, flatTable + flatEnd, flatTable);
}





/********************************************
 **********  Output Functions  **************
 *******************************************/

/*  These functions allow easy printing of
 *  matrix-like objects, given the number of
 *  rows and columns.
 */

void RprintDoubleMat(int rows, int cols, double *mat){
  int ii, jj;
  for(ii = 0 ; ii < rows ; ii++){
    for(jj = 0 ; jj < cols ; jj++){
      Rprintf("%.3f ",mat[ii*cols + jj]);
    }
    Rprintf("\n");
  }
}

void RprintIntMat(int rows, int cols, int *mat){
  int ii, jj;
  for(ii = 0 ; ii < rows ; ii++){
    for(jj = 0 ; jj < cols ; jj++){
      Rprintf("%d ",mat[ii*cols + jj]);
    }
    Rprintf("\n");
  }
}

//  Saves the current MMSBM state to the flatTable
void rdirichlet(int k, double *alpha, double *x){
  int ii;
  double draw[k];
  double total = 0.0;

  for(ii=0 ; ii < k ; ii++){
    draw[ii] = rgamma(alpha[ii],1);
    //draw[ii] = gsl_ran_gamma(rng,alpha[ii],1);
    total = total + draw[ii];
  }

  for(ii = 0; ii < k; ii++){
    x[ii] = draw[ii] / total;
  }
}


void getMeanVar(double *vec, int lower, int upper,
		double *Mean_t, double *Var_t, double * Len_t){


  double Var;
  double Mean = 0.0;
  double Len = upper - lower;
  double SumSq = 0.0;
  int ii;
  for(ii = lower ; ii < upper ; ii++){
    Mean = Mean + vec[ii];
    SumSq = SumSq + vec[ii] * vec[ii];
  }
  Mean = Mean / Len;
  Var = ((SumSq/Len) - (Mean * Mean))*(Len / (Len-1));
  *Mean_t = Mean;
  *Var_t = Var;
  *Len_t = Len;

}

int convergenceCheck(double *logLik, int total, double qq){

  double llMean1, llVar1, llLen1;
  double llMean2, llVar2, llLen2;
  double sdTot;

  //Rprintf("lower1 = %d, upper1 = %d\n",0,total/10);
  getMeanVar(logLik,0,total/10,
	     &llMean1, &llVar1, & llLen1);
  //Rprintf("lower2 = %d, upper2 = %d\n",total/2,total);
  getMeanVar(logLik,total/2,total,
	     &llMean2, &llVar2, & llLen2);
  sdTot = sqrt((llVar1/llLen1) + (llVar2/llLen2));

  //Rprintf("Mean 1 = %f \t Mean 2 = %f \t sd.est = %f\n",
  //llMean1,llMean2,sdTot);


  if((llMean2 - llMean1) <= (qq * sdTot)){
    if((llMean1 - llMean2) <= (qq * sdTot)){
      return(1);
    }
  }
  return(0);

}


void colSums(std::vector<std::vector<double> > const &mat, double *totals){
   int ii, kk;

   int nrow = mat.size();
   if(nrow < 1){
      std::cout << "Matrix has no rows" << std::endl;
      return;
   }
   int ncol = mat[0].size();
   if(ncol < 1){
      std::cout << "Matrix has no columns" << std::endl;
      return;
   }

   for(kk = 0 ; kk < ncol ; kk++){
      totals[kk] = 0.0;
      for(ii = 0 ; ii < nrow ; ii++){
	 totals[kk] += mat[ii][kk];
      }
   }
   return;
}


void colSums(std::vector<std::vector<double> > const &mat,
	     std::vector<double> &totals){
   int ii, kk;

   int nrow = mat.size();
   if(nrow < 1){
      std::cout << "Matrix has no rows" << std::endl;
      return;
   }
   int ncol = mat[0].size();
   if(ncol < 1){
      std::cout << "Matrix has no columns" << std::endl;
      return;
   }

   for(kk = 0 ; kk < ncol ; kk++){
      totals[kk] = 0.0;
      for(ii = 0 ; ii < nrow ; ii++){
	 totals[kk] += mat[ii][kk];
      }
   }
   return;
}





void colSums(double *mat, int rows, int cols, double *totals){
  int ii, jj, rr;
  for(ii = 0; ii < cols ; ii++){
    totals[ii] = 0.0;
  }

  for(ii = 0 ; ii < rows ; ii++){
    rr = ii * cols;
    for(jj = 0 ; jj < cols; jj++){
      totals[jj] = totals[jj] + mat[rr + jj];
    }
  }
}

void colSums(double *mat, int rows, int cols, std::vector<double> totals){
  int ii, jj, rr;
  for(ii = 0; ii < cols ; ii++){
    totals[ii] = 0.0;
  }

  for(ii = 0 ; ii < rows ; ii++){
    rr = ii * cols;
    for(jj = 0 ; jj < cols; jj++){
      totals[jj] = totals[jj] + mat[rr + jj];
    }
  }
}

void rowSums(double *mat, int rows, int cols, double *totals){
  int ii, jj, rr;
  for(jj = 0 ; jj < rows ; jj++){
    totals[jj] = 0.0;
    rr = jj * cols;
    for(ii = 0 ; ii < cols ; ii++){
      totals[jj] = totals[jj] + mat[rr + ii];
    }
  }
}

void rowSums(std::vector<std::vector<double> > const &mat,
	     std::vector<double> &totals){
   int ii, kk;

   int nrow = mat.size();
   if(nrow < 1){
      std::cout << "Matrix has no rows" << std::endl;
      return;
   }
   int ncol = mat[0].size();
   if(ncol < 1){
      std::cout << "Matrix has no columns" << std::endl;
      return;
   }

   for(kk = 0 ; kk < ncol ; kk++){
      totals[kk] = 0.0;
      for(ii = 0 ; ii < nrow ; ii++){
	 totals[kk] += mat[kk][ii];
      }
   }
   return;
}




void normalizeVec(double *vec, int ll){
  int ii;
  double total = 0.0;
  for(ii = 0 ; ii < ll ; ii++){
    total = total + vec[ii];
  }
  for(ii = 0 ; ii < ll ; ii++){
    vec[ii] = vec[ii] / total;
  }
}

void normalizeVec(std::vector<double> &vec){
   std::vector<double>::iterator it;
   double vec_total = 0.0;
   for(it = vec.begin(); it != vec.end() ; it++){
      vec_total += *it;
   }
   for(it = vec.begin(); it != vec.end() ; it++){
      *it /= vec_total;
   }

}

void logZeroFix(double *vec, int ll){
  int ii;
  double total = 0.0;

  for(ii = 0 ; ii < ll ; ii++){
    if(vec[ii] <= MIN_LOG){
      vec[ii] = MIN_LOG;
    }
    total = total + vec[ii];
  }

  if(total > 1.0){
    // Renormalize
    for(ii = 0 ; ii < ll ; ii++){
      vec[ii] = vec[ii] / total;
    }
  }
}

void logZeroFix(std::vector<double> &vec){
   std::vector<double>::iterator it;
   double vec_total = 0.0;

   for(it = vec.begin() ; it != vec.end() ; it++){
      if(*it <= MIN_LOG){
	 *it = MIN_LOG;
      }
      vec_total += *it;
   }

   if(vec_total > 1.0){
      // Renormalize
      for(it = vec.begin() ; it != vec.end() ; it++){
	 *it /= vec_total;
      }
   }

}

double logCheck(double val){
  if(val > MAX_LOG){
    return(MAX_LOG);
  }else if(val < MIN_LOG){
    return(MIN_LOG);
  }else{
    return(val);
  }
}

/*  OLD VERSION
void logZeroFix(double *vec, int ll){
  int ii;
  double fixVal = 0.0;
  double zeroCount = 0.0;
  for(ii = 0 ; ii < ll ; ii++){
    if(vec[ii] <= MIN_LOG){
      zeroCount = zeroCount + 1.0;
      fixVal  = fixVal + (MIN_LOG - vec[ii]);
    }
  }
  if(fixVal > 0.0){
    fixVal = fixVal / (ll - zeroCount);
    for(ii = 0 ; ii < ll ; ii++){
      if(vec[ii] <= MIN_LOG){
	vec[ii] = MIN_LOG;
      }else{
	vec[ii] = vec[ii] - fixVal;
      }
    }
  }
}

*/

/*
void updateFlatTable(int iter, int nn, int dd, double *BB, double *PP,
		     double *flatTable){
  int ii, offset;
  offset = iter * (dd * (nn + dd));
  for(ii = 0 ; ii < dd*dd ; ii++){
    flatTable[offset + ii] = BB[ii];
  }
  offset = offset + dd*dd;
  for(ii = 0 ; ii < nn*dd ; ii++){
    flatTable[offset + ii] = PP[ii];
  }

}

*/


/*
//  Prints all of the matrices involved in an MMSBM step
void printMMSBM(int nn, int dd, double *BB, double *PP,
int *sendMat, int *recMat){

printf("\nBlock Matrix\n");
printDoubleMat(dd,dd,BB);

printf("\nPP Matrix\n");
printDoubleMat(nn,dd,PP);

printf("\nSender Matrix:\n");
printIntMat(nn,nn,sendMat);

printf("Receiver Matrix:\n");
printIntMat(nn,nn,recMat);

}
*/
/*
  void printTableMMSBM(int nn, int dd, double *BB, double *PP,
  int *sendMat, int *recMat){
  int ii;
  for(ii = 0 ; ii < dd*dd ; ii++){
  //    printf("%f\t",BB[ii]);
  myfile << BB[ii] << "\t";
  }
  for(ii = 0 ; ii < nn*dd ; ii++){
  //    printf("%f\t",PP[ii]);
  myfile << PP[ii] << "\t";
  }
  myfile << "\n";
  }
*/




