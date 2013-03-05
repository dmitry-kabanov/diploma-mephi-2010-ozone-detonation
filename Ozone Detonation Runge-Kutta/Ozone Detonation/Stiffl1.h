// Stiffl1.h: interface for the Stiffl class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STIFFL1_H__F45EFB30_3AAA_4BAD_8D88_AB2AE2BB2757__INCLUDED_)
#define AFX_STIFFL1_H__F45EFB30_3AAA_4BAD_8D88_AB2AE2BB2757__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef min
	#define min(a, b)  (((a) < (b)) ? (a) : (b)) 
	#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#endif

class Stiffl  
{
public:
	Stiffl(int NYDIM_PAR, double *  pDInitialValues,double t_begin, double t_end, double t_step_begin);
	virtual ~Stiffl();
	int STIFFL();
	int Call_Calculation_of_relaxation_Time();
	void Settfin(double dPar);
	double Gettfin();
	static void SetMAX_ITERATION( double); 
	static void SetANOISE( double); 
	static void SetDELTZR( double );
	static void SetMAXITE( int );
	static void SetMAXFAI( int );
	static void SetRMXINI( double );
	static void SetRMXNOR( double );
	static void SetRMXFAI( double );
	static void SetIDELAY( int );
	static void SetRHCORR( double );
	static void SetRHERR3( double );
	static void SetTHRSHL( double );
	static void SetRCTEST( double );
	static void SetBIAS1( double);
	static void SetBIAS2( double );
	static void SetBIAS3( double );
	static void SetCMIN( double )  ;

	void static SetMAXDER(int);
	void static SetNYDIM(int);

	static void SetHMIN( double )  ;
	static void Seteps( double )  ;
	static void SetMFF( char )  ;

	static double GetMAX_ITERATION();
	static double GetANOISE(); 
	static double GetDELTZR( );
	static int GetMAXITE(  );
	static int GetMAXFAI(  );
	static double GetRMXINI(  );
	static double GetRMXNOR(  );
	static double GetRMXFAI(  );
	static int GetIDELAY(  );
	static double GetRHCORR(  );
	static double GetRHERR3(  );
	static double GetTHRSHL(  );
	static double GetRCTEST(  );
	static double GetBIAS1( );
	static double GetBIAS2(  );
	static double GetBIAS3(  );
	static double GetCMIN(  )  ;

	static int GetMAXDER();
	static int GetNYDIM();
	static double Geteps(  );
	static double GetHMIN(  );
	static char GetMFF(  )  ;
	void SetConc(double *);
	void GetConc(double *);
	void SetInitialTimeStep(double h);
	void SetInitialTime(double par);
	double GetInitialTime();
	void SetFinalTime(double par);
protected:
	double static ANOISE; 
	double static DELTZR ;
	int static MAXITE ;
	int static MAXFAI ;
	double static RMXINI ;
	double static RMXNOR ;
	double static RMXFAI ;
	int static IDELAY ;
	double static RHCORR ;
	double static RHERR3 ;
	double static THRSHL ;
	double static RCTEST ;
	double static BIAS1  ;
	double static BIAS2  ;
	double static BIAS3  ;
	double static CMIN   ;
	double static MAX_ITERATION;

	int static MAXDER;
	int static NYDIM;

	char static	MFF;
	double static	HMIN;
	double static	eps;
	long double STPF(long double X,long double Y);

	
    double* FSAVE;
    double* Z;
    double* VERROR;
    double* YMAX;        
	double** Y;
	double ***A;

	int  *IPIV,n,i,j,k,l,NQ,IDOUB;//IPIV[NYDIM]
	long double h,t,tfin,EL[6],
		   RC,RH,RMAX,SUM;
	char ped,EVALJA,CONVER,ISIGN,JSTART;
	
	char TREUG();
	void PREDIC();
	void SOL();
	void RESCAL();
	void COSET();


	virtual int DIFFUN( double **YY, double *F) = 0;
	virtual void PEDERV() =0;
	virtual int IFNSH()  =0;
    void prepareStiffl();

struct
	{
		long double EDN,E,EUP,BND;
	} ST;

	
public:
	int nNumberOfIterations;
};







#endif // !defined(AFX_STIFFL1_H__F45EFB30_3AAA_4BAD_8D88_AB2AE2BB2757__INCLUDED_)
