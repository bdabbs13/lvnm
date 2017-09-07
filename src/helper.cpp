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

extern "C" {
    void RGammaPrior(int *total, int *thin, int *burnin,
		     double *alpha_init, double *beta_init,
		     double *p, double *q, double *r, double *s,
		     double *prop_sd, double *alpha, double *beta){

	GetRNGstate();

	rGammaPrior(*total,*thin,*burnin,*alpha_init,*beta_init,
		    *p,*q,*r,*s,*prop_sd,alpha,beta);

	PutRNGstate();

    }


    void RGammaPriorNorm(int *total, int *thin, int *burnin,
			 double *alpha_init, double *beta_init,
			 double *p, double *q, double *r, double *s,
			 double *prop_sd, double *alpha, double *beta){

	GetRNGstate();

	rGammaPriorNorm(*total,*thin,*burnin,*alpha_init,*beta_init,
			*p,*q,*r,*s,*prop_sd,alpha,beta);

	PutRNGstate();

    }


    void RGammaPriorMV(int *total, int *thin, int *burnin,
		       double *alpha_init, double *beta_init,
		       double *p, double *q, double *r, double *s,
		       double *prop_sd, double *alpha, double *beta){

	GetRNGstate();

	rGammaPriorMV(*total,*thin,*burnin,*alpha_init,*beta_init,
		      *p,*q,*r,*s,*prop_sd,alpha,beta);

	PutRNGstate();

    }


    void ldGammaPrior(double *alpha, double *beta,
		      double *p, double *q, double *r, double *s,
		      double *result){
	double log_num = (alpha[0] - 1)*log(p[0]) - beta[0]*q[0] + alpha[0]*s[0]*log(beta[0]);
	double log_den = r[0] * lgammafn(alpha[0]);
	// Rprintf("alpha = %.2f, lgamma(alpha) = %.2f\n",alpha[0],lgammafn(alpha[0]));
	result[0] = log_num - log_den;
    }


}

//  Function for drawing a single sample from the GammaPrior distribution
void rGammaPrior(double alpha_init, double beta_init,
		 double p, double q, double r, double s,
		 double prop_sd, double *alpha, double *beta){

    return rGammaPrior(1,1,0,alpha_init,beta_init,p,q,r,s,prop_sd,alpha,beta);

}

void rGammaPriorStep(double *alphabeta,
		     double p, double q, double r, double s,
		     double prop_sd){

    return rGammaPrior(alphabeta[0],alphabeta[1],
		       p,q,r,s,prop_sd,
		       alphabeta,&(alphabeta[1]));

}

//  Function for drawing multiple samples from the GammaPrior distribution
void rGammaPrior(int total, int thin, int burnin,
		 double alpha_init, double beta_init,
		 double p, double q, double r, double s,
		 double prop_sd, double *alpha, double *beta){

    double alpha_cur = alpha_init;
    double beta_cur = beta_init;
    double alpha_prop, beta_prop, la_prob, lunif;
    //   double alpha_prop_mean, beta_prop_mean,alpha_cur_mean,beta_cur_mean;
    double ll_cur, ll_prop, lp_cur, lp_prop;

    int ii, save_iter;
    int long_total = total * thin + burnin;

    for(ii  = 0 ; ii < long_total ; ii++){
	ll_cur = ldGammaPrior(alpha_cur,beta_cur,p,q,r,s);

	//  Drawing From Proposal Distribution
	alpha_prop = rgamma(alpha_cur/prop_sd,prop_sd);
	while(alpha_prop <= 0) alpha_prop = rgamma(alpha_cur/prop_sd,prop_sd);
	beta_prop = rgamma(beta_cur/prop_sd,prop_sd);
	while(beta_prop <= 0) beta_prop = rgamma(beta_cur/prop_sd,prop_sd);

	ll_prop = ldGammaPrior(alpha_prop,beta_prop,p,q,r,s);

	//  Calculating Probability of Proposal from Current
	lp_prop = dgamma(alpha_prop,alpha_cur/prop_sd,prop_sd,1);
	lp_prop += dgamma(beta_prop,beta_cur/prop_sd,prop_sd,1);

	//  Calculating Probability of Current from Proposal
	lp_cur = dgamma(alpha_cur,alpha_prop/prop_sd,prop_sd,1);
	lp_cur += dgamma(beta_cur,beta_prop/prop_sd,prop_sd,1);

	//  Calculating Acceptance Probability
	la_prob = lp_cur - lp_prop + ll_prop - ll_cur;

	//  Accepting or Rejecting
	lunif = log(unif_rand());
	if(lunif < la_prob){
	    alpha_cur = alpha_prop;
	    beta_cur = beta_prop;
	}

	//  Saving Thinned/Burned In Draws
	if(ii >= burnin && (((ii - burnin) % thin) == 0)){
	    save_iter = (ii - burnin)/thin;
	    alpha[save_iter] = alpha_cur;
	    beta[save_iter] = beta_cur;
	}

    }
}


void GammaPriorMHDraw(double &alpha_cur, double &beta_cur,
		      double p, double q, double r, double s,
		      double &alpha_sd, double &beta_sd, double rho){


    //  Drawing From Truncated Normal
    double alpha_prop = -1;
    double beta_prop = -1;
    int draw_iter = 0;
    while(alpha_prop <= 0 || beta_prop <= 0){
	alpha_prop = rnorm(alpha_cur,alpha_sd);
	if(alpha_prop > 0){
	    double beta_prop_mean = beta_cur + rho*(alpha_prop - alpha_cur)*beta_sd/alpha_sd;
	    double beta_prop_sd = sqrt(1 - rho*rho)*beta_sd;
	    beta_prop = rnorm(beta_prop_mean,beta_prop_sd);
	}
	draw_iter++;
	if(draw_iter > 1000){
	    Rprintf("ERROR\n");
	    break;
	}
    }

    double ll_cur = ldGammaPrior(alpha_cur,beta_cur,p,q,r,s);
    double ll_prop = ldGammaPrior(alpha_prop,beta_prop,p,q,r,s);

    double la_prob = ll_prop - ll_cur;

    //  Accepting or Rejecting
    double lunif = log(unif_rand());
    if(lunif < la_prob){
	alpha_cur = alpha_prop;
	beta_cur = beta_prop;
    }

}


double GetDrawSD(double sse, int nn, int dd, double epsilon){
    return 2.38 * sqrt(((sse / nn) + epsilon)/ dd);
}


void GammaPriorMHDraw(double alphabeta[2], int &nn,
		      double p, double q, double r, double s,
		      double covInfo[5],double epsilon){

    double alpha_sd = GetDrawSD(covInfo[0],nn,2,epsilon);
    double beta_sd = GetDrawSD(covInfo[1],nn,2,epsilon);
    double rho = covInfo[2]/sqrt(covInfo[0] * covInfo[1]);

    GammaPriorMHDraw(alphabeta[0],alphabeta[1],
		     p,q,r,s,alpha_sd,beta_sd,rho);
}






//  Function for drawing multiple samples from the GammaPrior distribution
void rGammaPriorNorm(int total, int thin, int burnin,
		     double alpha_init, double beta_init,
		     double p, double q, double r, double s,
		     double prop_sd, double *alpha, double *beta){

    double alpha_cur = alpha_init;
    double beta_cur = beta_init;
    //   double alpha_prop, beta_prop, la_prob, lunif;
    //   double beta_prop_mean, beta_prop_sd;
    //   double ll_cur, ll_prop;


    double alpha_sd = prop_sd;
    double beta_sd = prop_sd;
    double rho = 0.0;
    int update_iter = 1;
    int update_start = 100;
    int update_count = 0;

    double alpha_mean = 0 ,beta_mean = 0;
    double alpha_sum = 0, beta_sum = 0;
    double alpha_sse = 0 ,beta_sse = 0, ab_csse = 0;

    int ii, save_iter;
    int long_total = total * thin + burnin;

    for(ii  = 0 ; ii < long_total ; ii++){

	// Drawing from MH Sampler
	GammaPriorMHDraw(alpha_cur,beta_cur,p,q,r,s,alpha_sd,beta_sd,rho);

	update_count++;
	updateSingleSSE(alpha_cur,beta_cur,update_count,
			alpha_sse,beta_sse,ab_csse,
			alpha_mean,beta_mean);

	if(update_count >= update_start){

	    alpha_sd = 2.38*sqrt(alpha_sse / update_count / 2.0);
	    beta_sd = 2.38*sqrt(beta_sse / update_count / 2.0);
	    rho = ab_csse / sqrt(alpha_sse * beta_sse);

	}

	//  Saving Thinned/Burned In Draws
	if(ii >= burnin){
	    if(((ii - burnin) % thin) == 0){

		save_iter = (ii - burnin)/thin;
		alpha[save_iter] = alpha_cur;
		beta[save_iter] = beta_cur;
	    }
	}
    }

}





void getCovariance(double *v1, double *v2, int start, int end,
		   double &var1, double &var2, double &cov,
		   double &m1, double &m2){

    m1 = 0.0;
    m2 = 0.0;
    cov = 0.0;
    double bar1_2 = 0.0;
    double bar2_2 = 0.0;

    for(int ii = start ; ii < end ; ii++){
	m1 += v1[ii];
	m2 += v2[ii];
	bar1_2 += v1[ii] * v1[ii];
	bar2_2 += v2[ii] * v2[ii];
	cov += v1[ii] * v2[ii];
    }

    var1 = (bar1_2 - m1*m1/end) / (end - 1);
    var2 = (bar2_2 - m2*m2/end) / (end - 1);
    cov = (cov - (m1 * m2 / end)) /(end - 1);
}

void updateCovariance(double *v1, double *v2, int start, int end,
		      double &var1, double &var2, double &cov,
		      double &m1, double &m2){
    double delta1, delta2;
    if((end - start) < 1){
	return;
    }else if((end - start) < 2){

	delta1 = v1[start] - m1;
	delta2 = v2[start] - m2;
	cov = cov*(start-1)/start + delta1*delta2/end;

	m1 += delta1/end;
	delta1 *= (v1[start] - m1);
	var1 = ((var1*(start - 1)) + delta1) / start;

	m2 += delta2/end;
	delta2 *= (v2[start] - m2);
	var2 = ((var2*(start - 1)) + delta2) / start;

	return;
    }else{
	return;
    }
}





void getSSE(double *v1, double *v2, int start, int end,
	    double &sse1, double &sse2, double &csse,
	    double &m1, double &m2){

    int nn = end - start;
    m1 = 0.0;
    m2 = 0.0;
    csse = 0.0;
    double bar1_2 = 0.0;
    double bar2_2 = 0.0;

    for(int ii = start ; ii < end ; ii++){
	m1 += v1[ii];
	m2 += v2[ii];
	bar1_2 += v1[ii] * v1[ii];
	bar2_2 += v2[ii] * v2[ii];
	csse += v1[ii] * v2[ii];
    }

    sse1 = (bar1_2 - m1*m1/nn);
    sse2 = (bar2_2 - m2*m2/nn);
    csse = (csse - (m1 * m2 / nn));
    m1 = m1/nn;
    m2 = m2/nn;

}

void updateSSE(double *v1, double *v2, int start, int end,
	       double &sse1, double &sse2, double &csse,
	       double &m1, double &m2){
    double delta1, delta2;
    if((end - start) < 1){
	return;
    }else if((end - start) < 2){
	updateSingleSSE(v1[start],v2[start],end,sse1,sse2,csse,m1,m2);
	return;
    }else{
	double sse1_new, sse2_new, csse_new, m1_new, m2_new;
	getSSE(v1,v2,start,end,sse1_new,sse2_new,csse_new,m1_new,m2_new);

	delta1 = m1_new - m1;
	m1 += delta1 * (end - start) / end;

	delta2 = m2_new - m2;
	m2 += delta2 * (end - start) / end;

	sse1 += sse1_new + delta1*delta1 * (end - start)*start/ end;
	sse2 += sse2_new + delta2*delta2 * (end - start)*start/ end;
	csse += csse_new + delta1*delta2 * (end-start)*start/end;

	return;
    }
}

void updateSingleSSE(const double &v1, const double &v2, int &nn,
		     double &sse1, double &sse2, double &csse,
		     double &m1, double &m2){


    double delta1 = v1 - m1;
    double delta2 = v2 - m2;

    m1 += delta1/nn;
    m2 += delta2/nn;

    csse = csse + delta1*(v2 - m2);

    delta1 *= (v1 - m1);
    delta2 *= (v2 - m2);

    sse1 = sse1 + delta1;
    sse2 = sse2 + delta2;

    return;

}



void updateSingleSSE(const double v[2], int &nn,
		     double covInfo[5]){

    updateSingleSSE(v[0],v[1],nn,
		    covInfo[0],covInfo[1],covInfo[2],
		    covInfo[3],covInfo[4]);
    return;

}


// void updateSingleSSE(double v1, double v2, int nn,
// 		     double &sse1, double &sse2, double &csse,
// 		     double &sum1, double &sum2){

//    double delta1, delta2;
//    if(nn > 1){
//       delta1 = v1 - sum1/(nn-1);
//       delta2 = v2 - sum2/(nn-1);
//    }else{
//       delta1 = v1;
//       delta2 = v2;
//    }
//    sum1 += delta1;
//    sum2 += delta2;

//    csse = csse + delta1*(v2 - sum2/nn);

//    delta1 *= (v1 - sum1/nn);
//    delta2 *= (v2 - sum2/nn);

//    sse1 = sse1 + delta1;
//    sse2 = sse2 + delta2;

//    return;

// }


void rGammaPriorMV(int total, int thin, int burnin,
		   double alpha_init, double beta_init,
		   double p, double q, double r, double s,
		   double prop_sd, double *alpha, double *beta){

    double mu_cur = alpha_init/beta_init;
    double sigma_cur = mu_cur/beta_init;
    double mu_prop, sigma_prop, la_prob, lunif;
    double mu_prop_mean, sigma_prop_mean,mu_cur_mean,sigma_cur_mean;
    double ll_cur, ll_prop, lp_cur, lp_prop;

    double mu_sd = prop_sd;
    double sigma_sd = prop_sd;

    int ii, save_iter;
    int long_total = total * thin + burnin;

    for(ii  = 0 ; ii < long_total ; ii++){
	ll_cur = ldGammaPrior(mu_cur*mu_cur/sigma_cur,mu_cur/sigma_cur,p,q,r,s);

	//  Drawing From Proposal Distribution
	mu_prop = rgamma(mu_cur/prop_sd,prop_sd);
	sigma_prop = rgamma(sigma_cur/prop_sd,prop_sd);
	ll_prop = ldGammaPrior(mu_prop*mu_prop/sigma_prop,mu_prop/sigma_prop,p,q,r,s);

	//  Calculating Probability of Proposal from Current
	lp_prop = dgamma(mu_prop,mu_cur/prop_sd,prop_sd,1);
	lp_prop += dgamma(sigma_prop,sigma_cur/prop_sd,prop_sd,1);

	//  Calculating Probability of Current from Proposal
	lp_cur = dgamma(mu_cur,mu_prop/prop_sd,prop_sd,1);
	lp_cur += dgamma(sigma_cur,sigma_prop/prop_sd,prop_sd,1);

	//  Calculating Acceptance Probability
	la_prob = lp_cur - lp_prop + ll_prop - ll_cur;

	//  Accepting or Rejecting
	lunif = log(unif_rand());
	if(lunif < la_prob){
	    mu_cur = mu_prop;
	    sigma_cur = sigma_prop;
	}

	//  Saving Thinned/Burned In Draws
	if(ii >= burnin && (((ii - burnin) % thin) == 0)){
	    save_iter = (ii - burnin)/thin;
	    alpha[save_iter] = mu_cur*mu_cur/sigma_cur;
	    beta[save_iter] = mu_cur/sigma_cur;
	    // if((save_iter % 100) == 0){
	    //    mu_sd = 2.38*2.38*var(alpha,save_iter)/2.0;
	    //    sigma_sd = 2.38*2.38*var(beta,save_iter)/2.0;
	    //    Rprintf("mu_sd = %.4f, sigma_sd = %.4f\n",
	    // 	    mu_sd,sigma_sd);
	    // }
	}

    }
}


double var(double *vec, int end){
    double bar = 0.0;
    double bar_2 = 0.0;
    for(int ii = 0 ; ii < end ; ii++){
	bar += vec[ii];
	bar_2 += (vec[ii] * vec[ii]);
    }
    bar = bar / end;
    bar_2 = bar_2 / end;
    return (bar_2 - bar*bar)*end/(end - 1);
}


double ldGammaPrior(double alpha, double beta,
		    double p, double q, double r, double s){
    double log_num = (alpha - 1)*log(p) - beta*q + alpha*s*log(beta);
    double log_den = r * lgammafn(alpha);

    return log_num - log_den;
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

    getMeanVar(logLik,0,total/10,
	       &llMean1, &llVar1, & llLen1);
    getMeanVar(logLik,total/2,total,
	       &llMean2, &llVar2, & llLen2);
    sdTot = sqrt((llVar1/llLen1) + (llVar2/llLen2));

    if((llMean2 - llMean1) <= (qq * sdTot)){
	if((llMean1 - llMean2) <= (qq * sdTot)){
	    return(1);
	}
    }
    return(0);

}

bool convergenceCheck(std::vector<double> logLik, double qq){
    double m1, v1, l1;
    double m2, v2, l2;
    int total = logLik.size();

    getMeanVar(logLik,0,total/10,    m1,v1,l1);
    getMeanVar(logLik,total/2,total, m2,v2,l2);

    double sdTot = sqrt((v1/l1) + (v2/l2));
    if(((m2 - m1) <= (qq * sdTot)) &&
       ((m1 - m2) <= (qq * sdTot))){
	return true;
    }else{
	return false;
    }
}

void getMeanVar(std::vector<double> vec, int lower, int upper,
		double &mm, double &vv, double &ll){

    int ii;
    double SS;

    ll = upper - lower;
    mm = 0.0;
    SS = 0.0;
    for(ii = lower ; ii < upper ; ii++){
	mm = mm + vec[ii];
	SS = SS + vec[ii] * vec[ii];
    }
    mm = mm / ll;
    vv = ((SS/ll) - (mm*mm))*(ll / (ll-1));

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




/////  MCMC TEMPLATE FUNCTION  /////



