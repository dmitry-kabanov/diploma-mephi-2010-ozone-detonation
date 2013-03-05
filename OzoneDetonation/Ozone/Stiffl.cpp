// Stiffl1.cpp: implementation of the Stiffl class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"

#include "Stiffl.h"
#include <cmath>
using namespace std;

//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif 

double  Stiffl:: MAX_ITERATION= 100000;

double Stiffl:: ANOISE= 2.0E-16;
double Stiffl:: DELTZR =1.0E-70;
int    Stiffl:: MAXITE =3;
int    Stiffl:: MAXFAI =3;
double Stiffl:: RMXINI =1.0E+4;
double Stiffl:: RMXNOR =10.0E0;
double Stiffl:: RMXFAI =2.0E0;
int    Stiffl:: IDELAY =9;
double Stiffl:: RHCORR =0.25;
double Stiffl:: RHERR3 =0.1E0;
double Stiffl:: THRSHL =0.0;
double Stiffl:: RCTEST =0.3;
double Stiffl:: BIAS1  =1.3;
double Stiffl:: BIAS2  =1.2;
double Stiffl:: BIAS3  =1.4;
double Stiffl:: CMIN   =1.0E3;

int    Stiffl:: MAXDER =5;
int Stiffl::NYDIM=256;

char		 Stiffl::MFF=0;
double		Stiffl::HMIN=1e-15;
double		Stiffl::eps=1e-5;

	void Stiffl:: SetMAX_ITERATION( double dPar){MAX_ITERATION = dPar; }
	void Stiffl:: SetANOISE( double dPar){ANOISE = dPar; }
	void Stiffl:: SetDELTZR( double dPar ){DELTZR = dPar; }
	void Stiffl:: SetMAXITE( int dPar ){MAXITE = dPar; }
	void Stiffl:: SetMAXFAI( int dPar ){MAXFAI = dPar; }
	void Stiffl:: SetRMXINI( double dPar ){RMXINI = dPar; }
	void Stiffl:: SetRMXNOR( double dPar ){RMXNOR = dPar; }
	void Stiffl:: SetRMXFAI( double dPar ){RMXFAI = dPar; }
	void Stiffl:: SetIDELAY( int dPar ){ IDELAY= dPar; }
	void Stiffl:: SetRHCORR( double dPar ){RHCORR = dPar; }
	void Stiffl:: SetRHERR3( double dPar ){RHERR3 = dPar; }
	void Stiffl:: SetTHRSHL( double dPar ){THRSHL = dPar; }
	void Stiffl:: SetRCTEST( double dPar ){RCTEST = dPar; }
	void Stiffl:: SetBIAS1( double dPar){ BIAS1= dPar; }
	void Stiffl:: SetBIAS2( double dPar ){BIAS2 = dPar; }
	void Stiffl:: SetBIAS3( double dPar ){BIAS3 = dPar; }
	void Stiffl:: SetCMIN( double dPar )  { CMIN= dPar; }

	void Stiffl:: SetMAXDER(int dPar){ MAXDER= dPar; }
	void Stiffl:: SetNYDIM(int dPar){ NYDIM= dPar; }

	void Stiffl:: SetHMIN( double dPar )  { HMIN= dPar; }
	void Stiffl:: Seteps( double dPar )  { eps= dPar; }
	void Stiffl:: SetMFF( char dPar )  { MFF= dPar; }

	double Stiffl:: GetMAX_ITERATION(){ return MAX_ITERATION   ;} 
	double Stiffl:: GetANOISE(){ return ANOISE   ;} 
  	double Stiffl:: GetDELTZR( ){ return DELTZR   ;}
	int Stiffl:: GetMAXITE(  ){ return MAXITE   ;}
	int Stiffl:: GetMAXFAI(  ){ return MAXFAI   ;}
	double Stiffl:: GetRMXINI(  ){ return RMXINI   ;}
	double Stiffl:: GetRMXNOR(  ){ return RMXNOR   ;}
	double Stiffl:: GetRMXFAI(  ){ return RMXFAI   ;}
	int Stiffl:: GetIDELAY(  ){ return  IDELAY  ;}
	double Stiffl:: GetRHCORR(  ){ return RHCORR   ;}
	double Stiffl:: GetRHERR3(  ){ return RHERR3   ;}
	double Stiffl:: GetTHRSHL(  ){ return THRSHL   ;}
	double Stiffl:: GetRCTEST(  ){ return  RCTEST  ;}
	double Stiffl:: GetBIAS1( ){ return BIAS1   ;}
 	double Stiffl:: GetBIAS2(  ){ return BIAS2   ;}
	double Stiffl:: GetBIAS3(  ){ return BIAS3   ;}
	double Stiffl:: GetCMIN(  )  { return CMIN   ;}

	int Stiffl:: GetMAXDER(){ return MAXDER   ;}
	int Stiffl:: GetNYDIM(){ return NYDIM   ;}

	double Stiffl:: GetHMIN(  )  { return HMIN   ;}
	double Stiffl:: Geteps(  )  { return eps   ;}
	char Stiffl:: GetMFF(  )  { return MFF   ;}
	void Stiffl:: Settfin(double dPar){ tfin= dPar; }
	double Stiffl:: Gettfin(){ return tfin; }
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Stiffl::Stiffl(int NYDIM_PAR,double*  pdInitialValues, double t_begin, double t_end, double t_step_begin)
{
	
		NYDIM=NYDIM_PAR+1;
		t = t_begin;
		tfin = t_end;
		h = t_step_begin;
		n=NYDIM_PAR-1;
		nNumberOfIterations=0;

        FSAVE= new double [NYDIM];
		for(int i=0;i<NYDIM;i++)
			FSAVE[i]=0;
        Z= new double [NYDIM];
		for( i=0;i<NYDIM;i++)
			Z[i]=0;
        VERROR= new double [NYDIM];
			for( i=0;i<NYDIM;i++)
				VERROR[i]=0;
        YMAX= new double [NYDIM];
		for( i=0;i<NYDIM;i++)
				YMAX[i]=0;
	
        Y = new double* [MAXDER+1];

		for( i=0;i<MAXDER+1;i++)
			Y[i] = new double [NYDIM];

		for(int j=0;j<MAXDER+1;j++)
			for( i=0;i<NYDIM;i++)
				(Y[j])[i]=0.0;
		A    = new double ** [NYDIM];
		for(   i=0;i<NYDIM;i++)
			A[i] =  new double *  [NYDIM];
		for(   i=0;i<NYDIM;i++)
			*A[i] =  new double   [NYDIM];

		for( j=0;j<NYDIM;j++)
			for( i=0;i<NYDIM;i++)
				(*A[i])[j]=0.0;
		
		
		IPIV =  new int   [NYDIM];
			for( i=0;i<NYDIM;i++)
				IPIV[i] =0;

        for(i=0;i<NYDIM;i++)
		(Y[0])[i]=pdInitialValues[i];
		
		
	
		
}

Stiffl::~Stiffl()
{
	delete [] FSAVE;
	delete [] Z;
	delete [] VERROR;
	delete [] YMAX;
	delete [] IPIV;
	for(int i=0;i<MAXDER+1;i++)
			delete [] Y[i] ;
	delete [] Y;

	for(   i=0;i<NYDIM;i++)
		delete [] *A[i];
	for(   i=0;i<NYDIM;i++)
			delete [] A[i] ;
	delete [] A;


}

long double Stiffl::STPF(long double X,long double Y)
{
return(pow(X,Y));
}

char Stiffl::TREUG() /*A[]^[],n*/
/*****************************************************************/
// DEC IS CALLED BY STIFFL AND IMPLEMENTS LU-FACTORIZATION FOR
//LINEAR SYSTEM WITCH IS SOLVING BY THE METHOD OF GAUSS ELEMINATION.
/*****************************************************************/
{
int i,j,k,l,m;
double TT;
static double *RRZ, *ZRZ, *XRZ;
static double   **AAZ, **AMZ;

IPIV[n]=1;
for(k=0;k<n;k++)
    {
    AMZ=A[k];
    l=k+1;
    m=k;
    for(i=l;i<=n;i++) if(fabs((*AMZ)[i])>fabs((*AMZ)[m]) ) m=i;
    IPIV[k]=m;
    TT=(*AMZ)[m];
    if(m!=k)
       {
       IPIV[n]=-IPIV[n];
       (*AMZ)[m]=(*AMZ)[k];
       (*AMZ)[k]=TT;
       };
//    if(TT==0.0) goto A80;
	if(TT==0.0) return(1);
    for(i=l;i<=n;i++)
       {
       RRZ=&((*AMZ)[i]);
       (*RRZ)=-(*RRZ)/TT;
       }
    for(j=l;j<=n;j++)
       {
        AAZ=A[j];
        RRZ=&((*AAZ)[m]);
        ZRZ=&((*AAZ)[k]);
		TT =(*RRZ);
		(*RRZ)=(*ZRZ);
		(*ZRZ)=TT;
		if(TT!=0.0)
	   for(i=l;i<=n;i++)
			{
			  XRZ=&((*AAZ)[i]);
              (*XRZ)=(*XRZ)+(*AMZ)[i]*TT;
			}
	};
    };
	if((*(A[n]))[n]==0.0)
		return(1); //A80:return(1);
	return(0);
};


void Stiffl::SOL()  /*A[]^[],n,FSAVE[]*/
/*****************************************************************/
/*  SOL IS CALLED BY STIFFL AND DO THE REVERSE GAUSS ELIMENATION */
/*****************************************************************/
{
int  k,i;
long double TT;
static double  *RRZ, *ZRZ, *XRZ;
static double   **AMZ; //!!!!!!

for(k=0;k<n;k++)
    {
    XRZ=&(FSAVE[k]);
    ZRZ=&(FSAVE[IPIV[k]]);
    TT =(*ZRZ);
    (*ZRZ)=(*XRZ);
    (*XRZ)=TT;
    for(i=k+1;i<=n;i++)
       {
       RRZ=&(FSAVE[i]);
       (*RRZ)=(*RRZ)+(*(A[k]))[i]*TT;
       };
    };
for(k=n;k>=1;k--)
    {
    AMZ=A[k];
    RRZ=&(FSAVE[k]);
    (*RRZ)=(*RRZ)/(*AMZ)[k];
    TT=-(*RRZ);
    for(i=0;i<k;i++)
       {
       RRZ=&(FSAVE[i]);
       (*RRZ)=(*RRZ)+(*AMZ)[i]*TT;
       };
    };
FSAVE[0]=FSAVE[0]/(*(A[0]))[0];
};


void Stiffl:: RESCAL()
/*****************************************************************/
/* RESCAL IS CALLED BY STIFFL AND RESCALES ARRAY OF Y IF THE TIME*/
/* STEP WILL CHANGED				                 */
/*****************************************************************/
{
long double R1;
int  i,j;
static double  *RRZ;

RH=min(min(max(RH,HMIN/fabs(h)),fabs((tfin-t)/h)),RMAX);
R1=1.0;
for(j=1;j<=l;j++)
    {
    R1=R1*RH;
    for(i=0;i<=n;i++)
       {
       RRZ=&((Y[j])[i]);
       (*RRZ)=(*RRZ)*R1;
       }
    };
h=h*RH;
RC=RC*RH;
IDOUB=l+1;
};

void Stiffl::PREDIC()
/*******************************************************************/
/* PREDIC IS CALLED BY STIFFL AND CALCULATES VALUES OF PREDICTOR   */
/* UNDER EFFECTIVE MULTIPLICATION OF Y-ARRAY AND TRIANGULAR        */
/* PASCAL'S MATRIX.					           */
/*******************************************************************/
{
int i,j,j1;
static double  *RRZ;

for(j1=0;j1<=NQ;j1++)
 for(j=NQ;j>=j1;j--)
     if(ISIGN)
	 for(i=0;i<=n;i++)
	   {
	   RRZ=&((Y[j])[i]);
	   (*RRZ)=(*RRZ)+(Y[j+1])[i];
	   }
     else
	for(i=0;i<=n;i++)
	   {
	   RRZ=&((Y[j])[i]);
	   (*RRZ)=(*RRZ)-(Y[j+1])[i];
	   };
};


/*STIFFL needes in following parameters: n,Y[0][0..n],h,t,tfin,EPS,MFF*/
/*and memory for array A[0..n]^[0..n]*/
int Stiffl::STIFFL()
{
 static char KFLAG,ITER;
static int  NEWQ;
static long double PR1,PR2,D1,D,FACTOR,YJ,R,
		   PR3,HOLD,OLDL0,TREND,EPSJ,TOLD,YS,EPSCOM;
static double  *RRZ, *ZRZ;

//My!! Buff_F=Buff_Y;
nNumberOfIterations=0;
if(fabs(tfin-t)<fabs(h)) h=tfin-t;
EPSCOM=ANOISE/eps;
JSTART=0;
HMIN=DELTZR;
ped=0;

REPEAT:
	  nNumberOfIterations++;
	  if(nNumberOfIterations>MAX_ITERATION) return nNumberOfIterations;
      if(h==0.0) return 0;
      if ((t+h)>tfin) h=tfin-t;
//	 	printf("t %g h %g T %f conc %f\n",t,h,(Y[0])[9], (Y[0])[7]); 
/***********************************************************;*/
      KFLAG=0;
      TOLD=t;
      if(JSTART==0)
	 {
	 if(DIFFUN(&Y[0],FSAVE)) return -1;
	 for(i=0;i<=n;i++) 
		 (Y[1])[i]=FSAVE[i]*h;
	 NQ=0;
	 l=1;
	 IDOUB=l+1;
	 RMAX=RMXINI;
	 EPSJ=sqrtl(ANOISE);
	 TREND=1.0;
	 OLDL0=1.0;
	 RC=0.0;
	 HOLD=h;
	 EVALJA=1;
	 CONVER=0;
	 COSET();
	 };
      if((h!=HOLD)&(JSTART<0))
	 {
	 RH=h/HOLD;
	 h=HOLD;
	 RESCAL();
	 };
A200: if(fabs(RC-1.0)>RCTEST) EVALJA=1;
      t=t+h;
      ISIGN=1;
/************************************************/
      YS=CMIN;
      for(i=0;i<=n;i++)
	  {
	  if(i==n) YS=DELTZR;
	  YMAX[i]=fabs((Y[0])[i])+fabs((Y[1])[i]*h)+YS;
	  };
/************************************************/
      PREDIC();
A220: for(i=0;i<=n;i++) VERROR[i]=0;
      ITER=0;
      ped=EVALJA & MFF;
      if(DIFFUN(&Y[0],Z)) return -1;
      if(!EVALJA)  goto A460;
      if(MFF)
	 {
//My!!!!	 PEDERV();
	 ped=0;
	 R=-EL[0]*h;
	 for(i=0;i<=n;i++)
	    for(j=0;j<=n;j++)
	       (*RRZ)=(*(RRZ=&((*(A[j]))[i])))*R;
	 }
      else
	 {
	 for(j=0;j<=n;j++)
	     {
	     R=max(fabs(h*Z[j])*1.0E3*ANOISE,EPSJ*YMAX[j]);
	     YJ=(*(RRZ=&((Y[0])[j])));
	     (*RRZ)=YJ+R;
	     FACTOR=-EL[0]*h/R;
	     if(DIFFUN(&Y[0],FSAVE)) return -1;
	     for(i=0;i<=n;i++) (*(A[j]))[i]=(FSAVE[i]-Z[i])*FACTOR;
	     (*RRZ)=YJ;
	     };
	 };

      for(i=0;i<=n;i++)
	 (*RRZ)=(*(RRZ=&((*(A[i]))[i])))+1.0;
	  //!!!!!!
      EVALJA=0;
      CONVER=0;
      RC=1.0;
      if(TREUG()) 
		  throw 201;;
A460: for(i=0;i<=n;i++) FSAVE[i]=Z[i]*h-(Y[1])[i]-VERROR[i];
      SOL();
      D=0.0;
      for(i=0;i<=n;i++)
	 {
	  (*RRZ)=(*(RRZ=&(VERROR[i])))+(*(ZRZ=&(FSAVE[i])));
/*!!!*/	  D=D+((*ZRZ)/YMAX[i])*((*ZRZ)/YMAX[i]);
	  (*ZRZ)=(Y[0])[i]+EL[0]*(*RRZ);
	 };
      if(ITER!=0) TREND=max(0.9*TREND,D/D1);
      if(D*min(1.0,2.0*TREND)>ST.BND)
	 {
	 D1=D;
	 ITER++;
	 if(ITER!=MAXITE)
	    {
	    if(DIFFUN(&FSAVE,Z))
			return -1;
	    goto A460;
	    };
	 if(CONVER)
	    {
	    EVALJA=1;
	    goto A220;
	    };
	 t=TOLD;
	 RMAX=RMXFAI;
	 ISIGN=0;
	 PREDIC();
	 if(fabs(h)<=HMIN*1.00001) goto A820;
	 RH=RHCORR;
	 RESCAL();
	 goto A200;
	 };
      D=0.0;
/*!!!*/
      for(i=0;i<=n;i++) D=D+(VERROR[i]/YMAX[i])*(VERROR[i]/YMAX[i]);
      CONVER=1;
      if(D>ST.E)
	 {
	 KFLAG++;
	 t=TOLD;
	 ISIGN=0;
	 PREDIC();
	 RMAX=RMXFAI;
	 if(fabs(h)<=HMIN*1.00001E0)
	    {
	    KFLAG=1;
	    goto A840;
	    };
	 if(KFLAG<MAXFAI)
	    {
	    PR2=1.0/(STPF(D/ST.E,0.5/(l+1))*BIAS2+DELTZR);
	    if( NQ!=0 )
	       {
	       SUM=0.0;
/*!!!*/	       for(i=0;i<=n;i++) SUM=SUM+((Y[l])[i]/YMAX[i])*((Y[l])[i]/YMAX[i]);
	       PR1=1.0/(STPF(SUM/ST.EDN,0.5/(NQ+1))*BIAS1+DELTZR);
	       if(PR1>PR2) goto A640;
	       };
	    RH=PR2;
	    RESCAL();
	    goto A200;
A640:       l=NQ;
	    NQ--;
	    RH=PR1;
	    COSET();
	    RC=RC*EL[0]/OLDL0;
	    OLDL0=EL[0];
	    RESCAL();
	    goto A200;
	    };
	 RH=RHERR3;
	 RH=max(HMIN/fabs(h),RH);
	 h=h*RH;
	 if(DIFFUN(&Y[0],FSAVE)) return -1;
	 for(i=0;i<=n;i++) (Y[1])[i]=h*FSAVE[i];
	 EVALJA=1;
	 IDOUB=IDELAY;
	 if(NQ==0) goto A200;
	 NQ=0;
	 l=1;
	 COSET();
	 OLDL0=EL[0];
	 goto A200;
	 };
      KFLAG=0;
      for(j=0;j<=l;j++) for(i=0;i<=n;i++)
	  (*RRZ)=(*(RRZ=&((Y[j])[i])))+EL[j]*VERROR[i];
      if(IDOUB!=0)
	 {
	 IDOUB--;
	 if((IDOUB<=0)&(NQ!=MAXDER-1))
	    for(i=0;i<=n;i++) (Y[MAXDER])[i]=VERROR[i];
	 goto A840;
	 };
      PR3=DELTZR;
      if(NQ!=MAXDER-1)
	 {
	 SUM=0.0;
/*!!!*/	 for(i=0;i<=n;i++) SUM=SUM+((VERROR[i]-(Y[MAXDER])[i])/YMAX[i])*((VERROR[i]-(Y[MAXDER])[i])/YMAX[i]);
	 PR3=1.0/(STPF(SUM/ST.EUP,0.5/(l+2))*BIAS3+DELTZR);
	 };
      PR2=1.0/(STPF(D/ST.E,0.5/(l+1))*BIAS2+DELTZR);
      PR1=DELTZR;
      if(NQ!=0)
	 {
	 SUM=0.0;
/*!!!*/	 for(i=0;i<=n;i++) SUM=SUM+((Y[l])[i]/YMAX[i])*((Y[l])[i]/YMAX[i]);
	 PR1=1.0/(STPF(SUM/ST.EDN,0.5/(NQ+1))*BIAS1+DELTZR);
	 };
      if((PR3<=PR1)|(PR3<=PR2))
	 {
	 if(PR1<=PR2)
	    {
	    RH=PR2;
	    if(RH<THRSHL) goto A800; else goto A790;
	    };
	 NEWQ=NQ-1;
	 RH=PR1;
	 if(RH<THRSHL) goto A800; else goto A780;
	 };
      NEWQ=l;
      RH=PR3;
      if(RH<THRSHL) goto A800;
      for(i=0;i<=n;i++) (Y[NEWQ+1])[i]=VERROR[i]*EL[l]/(l+1);
A780: NQ=NEWQ;
      l=NQ+1;
      COSET();
      RC=RC*EL[0]/OLDL0;
      OLDL0=EL[0];
A790: RESCAL();
     goto A830;
	  RMAX=RMXNOR;
	  goto A840;	  
A800: IDOUB=IDELAY;
      goto A840;
A820: KFLAG=2;
      goto A840;

A830: RMAX=RMXNOR;
A840: HOLD=h;
      if(KFLAG)
	 {
	 JSTART=-1;
	 h=h/10.0;
	 }
      else
	 {
	 JSTART=NQ+1;
     if(0!=IFNSH()) return 11;
	 };

if(t<tfin) goto REPEAT;          
      return 0;   
}

void Stiffl::COSET()
/*********************************************************************/
// COSET IS CALLED BY STIFFL.
// VECTOR EL DETERMINES THE METHOD COEFFICIENTS.
// VECTOR TQ IS CONNECTES TIME STEP WITH ROUNDING ERROR.
/*********************************************************************/
{
long double *STM;
static long double PERTST[3][5]={
      {1,1,0.5,1./6.,1./24.},
      {2,4.5,22./3.,125./12.,13.7},
      {3,6,55./6.,12.5,1}};
static long double ELCNST[5][6]={
      {1,1,0,0,0,0},{2./3.,1,1./3.,0,0,0},
      {6./11.,1,6./11.,1./11.,0,0},{0.48,1,0.7,0.2,0.02,0},
      {120./274.,1,225./274.,85./274.,15./274.,1./274.}};

int  k;

for(k=0;k<=5;k++) EL[k]=ELCNST[NQ][k];
STM=(long double *)&ST;
for(k=0;k<=2;k++)
   {
   (*STM)=(eps*PERTST[k][NQ])*(eps*PERTST[k][NQ]);   /*!!!*/
   STM++;
   }
ST.BND=(0.5/(NQ+3)*eps)*(0.5/(NQ+3)*eps);
};
void Stiffl::SetConc(double *yy)
{
	for(int i=0;i<NYDIM;i++)
		(Y[0])[i]=yy[i];
}
void Stiffl::GetConc(double *yy)
{
	for(int i=0;i<NYDIM;i++)
		yy[i]=(Y[0])[i];
}
void Stiffl::SetInitialTimeStep(double par)
{
	h=par;
}
void Stiffl::SetInitialTime(double par)
{
	t=par;
}
void Stiffl::SetFinalTime(double par)
{
	tfin=par;
}
double Stiffl::GetInitialTime()
{
	return t;
}

void Stiffl::prepareStiffl()
{
    nNumberOfIterations=0;

    for(int i=0;i<NYDIM;i++)
        FSAVE[i]=0;

    for( i=0;i<NYDIM;i++)
        Z[i]=0;

    for( i=0;i<NYDIM;i++)
        VERROR[i]=0;

    for( i=0;i<NYDIM;i++)
        YMAX[i]=0;

    for(int j=0;j<MAXDER+1;j++)
        for( i=0;i<NYDIM;i++)
            (Y[j])[i]=0.0;

    for( j=0;j<NYDIM;j++)
        for( i=0;i<NYDIM;i++)
            (*A[i])[j]=0.0;

    for( i=0;i<NYDIM;i++)
        IPIV[i] =0;
}