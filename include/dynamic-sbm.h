/*****  wsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef DYNSBM_R
#define DYNSBM_R
#include <vector>




//  Definition for WSBM class
class CDynSBM {
 public:
    CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses,
	     int total, int mImpute);
    ~CDynSBM();


    // Loading Functions
    void LoadDataR(int *AdjMat, int *rTimeMap, double *rHours,
		   double *rHyperSender, double *rHyperReceiver,
		   double *rHyperBlockMat);

    void LoadStateR(double *rSenderEffects, double *rReceiverEffects,
		    double *rBlockEffects, int *rBlockMemb,
		    double *rPosteriorMemb, int update_mmb,
		    double *rPriorSender, double *rPriorReceiver,
		    double *rPriorBlockMat, double *rPriorBlockMemb,
		    double *rLogLik);

    /*****  MCMC Interface Functions  *****/
    //  MCMC Functions
    void step();
    void adapt();

    // Writing Function
    void write(int iter);

    // Log-Likleihood Function
    double LogLike();

    //  Retrieval Functions
    int GetNodes() const { return aNodes;}
    int GetBlocks() const { return aBlocks;}
    int GetTimes() const {return aTimes;}
    int GetTimeClasses() const {return aClasses;}

    // Printing Functions (Primarily for Debugging)
    void printPriors();
    void printAllWSBM(bool printNetworks);
    void printBlockMemb();
    void PrintCovariance(int tt, int ii);


 private:

    /**********  Parameters and Containers  **********/

    // Dimensional Parameters
    int aNodes;
    int aBlocks;
    int aTimes;
    int aClasses;
    int aTotal;

    // Flags
    int missingVal;
    bool update_mmb;
    bool aImputeFlag;
    bool is_BlockMat_logged;

    // Adaptive MCMC Parameters
    double mhEpsilon;
    double mhSD;
    double mhStart;


    /***** WSBM Objects  *****/
    std::vector<CWSBM *> aWsbmList;

    /***** Time Hash Functions *****/
    std::vector<double> aHours;
    std::vector<int> aTimeMap;

    /*****  Block Membership Params  *****/
    std::vector<int> aBlockMemb; // aNodes
    std::vector<std::vector<double> > aPosteriorMemb; // aNodes x aBlocks

    /*****  Prior Containers  *****/
    std::vector<double> aPriorBlockMemb;
    std::vector<std::vector<double *> *> aPriorSender;
    std::vector<std::vector<double *> *> aPriorReceiver;
    std::vector< std::vector<std::vector<double *> > *> aPriorBlockMat;

    /*****  Adaptive Sampling Containers  *****/
    int covCount;
    std::vector<std::vector<double *> *> aPriorSenderCov;
    std::vector<std::vector<double *> *> aPriorReceiverCov;
    std::vector< std::vector<std::vector<double *> > *> aPriorBlockMatCov;

    //  Hyperpriors
    double aHyperSender[4];
    double aHyperReceiver[4];
    double aHyperBlockMat[4];

    /*****  Writer Information  *****/
    enum WriterType { R_WRITER, TEXT_WRITER};
    WriterType aWriter;

    /*****  Pointers for R_WRITER  *****/
    double *aRPriorSender;
    double *aRPriorReceiver;
    double *aRPriorBlockMat;
    int *aRBlockMemb;
    double *aRPosteriorMemb;
    double *aRLogLike;


    /**********  Internal Functions  **********/

    // LoadData Functions
    void LoadTime(int *rTimeMap, double *rHours);
    void LoadAdjacencyMatrices(int *AdjMat);
    void LoadHyperPriors(double *rHyperSender, double *rHyperReceiver,
			 double *rHyperBlockMat);

    // LoadStateR Functions
    void LoadParameters(double *rSenderEffects, double *rReceiverEffects,
			double *rBlockEffects, int *rBlockMemb,
			double *rPosteriorMemb, int update_mmb);

    void LoadPriors(double *rPriorSender, double *rPriorReceiver,
		    double *rPriorBlockMat, double *rPriorBlockMemb);
    // Internal Loading Functions
    void LoadPriorSender();
    void LoadPriorReceiver();
    void LoadPriorBlockMat();

    void LoadLogLike(double *rLogLik);
    void PassReferences();

    // MCMC Functions
    void DrawPriors();
    // Subfunctions of Draw Priors
    void DrawPriorSender();
    void DrawPriorReceiver();
    void DrawPriorBlockMat();

    void DrawParameters();
    // Subfunction of DrawParameters
    void DrawBlockMemb();
    //void imputeMissingValues();  //Need to be updated

    // Log-Likelihood Functions
    double nodeLogLike(int ii);


    // Writing Functions for R
    void writeR(int iter);

    void writeRPriorSender(int iter);
    void writeRPriorReceiver(int iter);
    void writeRPriorBlockMat(int iter);
    void writeRBlockMemb(int iter);
    void writeRPosteriorMemb(int iter);
    void writeRLogLike(int iter);

    // Printing Functions
    void printPriorSender();
    void printPriorReceiver();
    void printPriorBlockMat();

};

//  Function for performing the MCMC algorithm
void dynSBMMCMC(CDynSBM *myDynSBM, int start, int total, int burnIn, int thin,
		int shift_size, int extend_max, double qq, int verbose);


#endif

