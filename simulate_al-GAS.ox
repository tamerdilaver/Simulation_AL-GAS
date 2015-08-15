#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iB;	 					//Repeats
static decl iSIZE;					//Size of time series
static decl iSTEPS;					//#Steps to divide the size
static decl iSIMS;					//# of Zt ~ N(0,1)
static decl iPLOT;
static decl iDIST;
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl iPARS;					//number of parameters
static decl vSTD_ALD;				// Zt ~ ALD(0, b, p)
static decl s_vY; 					//Simulated returns
//static decl	bSTD_ERROR;				//0 or 1 boalean



findicator(const vX){
	decl vReturn, iN;
	iN = sizerc(vX);

	if(rows(vX)!=1){
		vReturn = zeros(iN,1);
	}else{
		vReturn = zeros(1,iN);
	}
	
//	print(vReturn);
	
	for(decl i=0;i<iN;i++){
		if(vX[i]<=0){
			vReturn[i]=1;
		}
	}
	return vReturn; 
}

fSimALGAS(const dAlpha, const dBeta, const dOmega, const dGamma, const dP, const avReturns){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_ALD[i];
		vH[i+1] = dOmega + dBeta*vH[i] + dAlpha*((dP-findicator(vTemp[i]))*(2*sqrt(1-2*dP+2*sqr(dP))/((1-dP)*dP))*vTemp[i]*sqrt(vH[i])-2*vH[i]) ;
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vH;
	return 1;
}

main()
{

	decl vRanExp, vRanExp2, vRanBern, vRanALD, vRanN;
	decl dB;  //dB > 0
	decl dP;  //0 < dP < 1
	
	dB = 1;
	dP = 0.3;
	iDIST = 10000000;
	iSIZE = 1000;
	iPLOT = 1000;
	vRanExp = ranexp(iDIST, 1, 1);
	vRanExp2 = ranexp(iDIST, 1, 1);
//	vRanBern = ranbinomial(iDIST, 1, 1, dP);
	vRanALD = vRanExp/dP-vRanExp2/(1-dP);	 //~ALD(0,1,p)
	vRanALD = dB*vRanALD;					 //~ALD(0,b,p)

	vSTD_ALD = (1-dP)*dP/(dB*sqrt(1-2*dP+2*sqr(dP)))*vRanALD;
	print(variance(vSTD_ALD));
	vRanN = rann(iDIST,1);

	decl vReturns;
	fSimALGAS(0.1, 0.99, 0.01, 0.1, 0.45, &vReturns);
	//print(vReturns);
																		

	SetDrawWindow("0_SIM_asymmetric_laplace-GAS");
	DrawDensity(0, (vRanN)', {"Density zt~ N(0,1)"});
	DrawDensity(0, (vSTD_ALD)', {"Density zt~ st_ALD(0,b,p)"});
	Draw(1, (vSTD_ALD[0:iPLOT])');
	Draw(2, (vReturns)');
	ShowDrawWindow();
	//print(vRanALD);
}
