/*
**	Master Thesis Econometrics
**
**  Purpose:
**  	For some fixed parameter values simulate and estimate GARCH model parameters
**		with Maximum Likelikhood many times. s.t. Elog(alpha_0 z_t^2 + beta_0) = 0.
**
**  Date:
**    	30/05/2015
**
**  Author:
**	  	Tamer Dilaver
**
**	Supervisor:
**		Fransisco Blasques
**
*/


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
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl dPROB;
static decl iPARS;					//number of parameters
static decl vSTD_ALD;				// Zt ~ N(0,1)
static decl s_vY; 					//Simulated returns
static decl	bSTD_ERROR;				//0 or 1 boalean

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


/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/
fTransform2(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);

	return 1;
}

/*
**  Function:	Transform parameters back
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/
fTransformBack2(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adBeta, vTheta
**
**  Output: 	1 
*/
fGetPars2(const adBeta, const vTheta){

	adBeta[0] = exp(vTheta[0]);
	return 1;
}

/*
**  Function:	Calculate targetvalue -[ E log(alpha_0 z_t^2 + beta_0) ]^2 for given parameters
**
**  Input:		vTheta [parametervalues], adFunc [adres functionvalue], avScore [the score],  amHessian [hessianmatrix]
**
**  Output:		1
**
*/
fExpectation(const vTheta, const adFunc, const avScore, const amHessian){
	decl dBeta;
	fGetPars2( &dBeta, vTheta);

	//NOTICE: Since maxBFGS() is a function that maximimises I have squared the target function
	//and then I have put a negative sign in front of it. This ensures that maximising will
	//effectively search for the value for when the target function is zero.
	adFunc[0] = - (meanc(log(fabs(dALPHA*( 2*sqrt(1-2*dPROB+2*sqr(dPROB))/((1-dPROB)*dPROB)*(dPROB-findicator(vSTD_ALD)).*vSTD_ALD - 2)+ dBeta))))^2; 						
	return 1;
}

/*
**  Function:	Get the value for beta subject to Elog(alpha_0 z_t^2 + beta_0) = 0 for given parameters	and simulated Zt's
**
**  Input:		iSims, adBeta  
**
**  Output:		1
**
*/
fGetBeta(const iSims, const adBeta)
{
	decl vTheta, vThetaStart, vThetaStar, dFunc, iA;

	//initialise startparameter(s)
	vTheta = zeros(1,1);
	vTheta = <0.9>;			// dBeta
	vThetaStart = vTheta;

	//transform startparameter(s)
	fTransform2(&vThetaStar, vTheta);

	//maximise
	iA=MaxBFGS(fExpectation, &vThetaStar, &dFunc, 0, TRUE);
	
	//Transform thetasStar back
   	fTransformBack2(&vTheta, vThetaStar);

	print("\nStart & Optimal parameter(s) with alpha fixed at ",dALPHA," and ",iSIMS," simulations such that we get I(1). \n",
          "%r", { "dBeta"},
          "%c", {"thetaStart","theta"}, vThetaStart~vTheta);
	adBeta[0] = vTheta[0];

	return 1;
}



/*
**  Function:	Simulate GARCH returns for given parameters
**
**  Input:		dAlpha, dBeta, dOmega, dGamma, avReturns, iIteration [to get different Zt's]
**
**  Output:		1
**
*/

fSimALGAS(const dAlpha, const dBeta, const dOmega, const dGamma, const dP, const avReturns, const iIteration){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_ALD[(i + (iIteration*iSIZE))];
		//vH[i+1] = dOmega+ dBeta*vH[i] + dAlpha*sqr(vTemp[i]) ;
		vH[i+1] = dOmega + dBeta*vH[i] + dAlpha*((dP-findicator(vTemp[i]))*(2*sqrt(1-2*dP+2*sqr(dP))/((1-dP)*dP))*vTemp[i]*sqrt(vH[i])-2*vH[i]);
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
	return 1;
}

/*
**  Function:	Transform (start)parameters	  Alpha, Beta, Omega, Gamma Startingvalues
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0] = vTheta;

	avThetaStar[0][0] = log(vTheta[0]);
	avThetaStar[0][1] = log(vTheta[1]);
	avThetaStar[0][2] = log(vTheta[2]);
	avThetaStar[0][3] = log(vTheta[3]);
	avThetaStar[0][4] = log((vTheta[4])/(1-vTheta[4]));
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adGamma,, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega, const adGamma, const adP, const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);
	adOmega[0] = exp(vTheta[2]);
	adGamma[0] = exp(vTheta[3]);
	adP[0] = exp(vTheta[4])/(1+exp(vTheta[4]));

	return 1;
}

/*
**  Function:	Calculates average value loglikelihood for GARCH given parameter values
**
**  Input: 		vTheta [parameter values], adFunc [adress functievalue], avScore [the score], amHessian [hessian matrix]
**
**  Output:		1
**
*/

fLogLike_alGAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dAlpha, dBeta, dOmega, dGamma, dP;
	fGetPars( &dAlpha,  &dBeta, &dOmega,  &dGamma, &dP, vTheta);

	decl dS2 = dGamma;	//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = 0.5*log(1-2*dP+2*sqr(dP)) - 0.5*log(dS2) -sqrt((1-2*dP+2*sqr(dP)))/(sqrt(dS2)*(1-dP)*dP)*s_vY[i]*(dP-findicator(s_vY[i])); //Gaussian
						
			//GARCH recursion
			dS2 = dOmega + dBeta* dS2 +  dAlpha*((dP-findicator(s_vY[i]))*(2*sqrt(1-2*dP+2*sqr(dP))/((1-dP)*dP))*s_vY[i]*sqrt(dS2)-2*dS2); //s_vY[i]^2;
	}		                           //dAlpha*((dP-findicator(vTemp[i]))*(2*sqrt(1-2*dP+2*sqr(dP))/((1-dP)*dP))*vTemp[i]*sqrt(vH[i])-2*vH[i]);
	
	adFunc[0] = sumc(vLogEta)/sizerc(s_vY); //sumc(vLogEta)/(-2*sizerc(s_vY)); //Average
	return 1;
}

/*
**  Function:	Transform parameters back	Alpha, Beta, Omega, Gamma Startingvalues
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/

fTransformBack(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	avTheta[0][1] = exp(vThetaStar[1]);
	avTheta[0][2] = exp(vThetaStar[2]);
	avTheta[0][3] = exp(vThetaStar[3]);
	avTheta[0][4] = exp(vThetaStar[4])/(1+exp(vThetaStar[4]));
	return 1;
}

/*
**  Function:	calculate standard errors
**
**  Input: 		vThetaStar
**
**  Output: 	vStdErrors
*/
fSigmaStdError(const vThetaStar){

 		decl iN, mHessian, mHess, mJacobian, vStdErrors, vP;

		iN 			= sizerc(s_vY);
		Num2Derivative(fLogLike_alGAS, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

		return 	vStdErrors;
}

/*
**  Function:	Estimate Garch parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adGamma_hat
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimateGarch(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adGamma_hat, const adP_hat){

	//initialise parameter values
	decl vTheta = zeros(iPARS,1);
	vTheta = <0.1 ; 0.9 ; 0.05 ; 0.1; 0.5>; // Alpha, Beta, Omega, Gamma, P Startingvalues
	decl vThetaStart = vTheta;

	//globalalize returns and vectorize true pars
	s_vY = vReturns;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_alGAS, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return alpha, beta, omega and gamma
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];
	adOmega_hat[0] = vTheta[2];
	adGamma_hat[0] = vTheta[3];
	adP_hat[0] = vTheta[4];

	if(bSTD_ERROR){		//only do this for fMonteCarlo2
		decl vSigmaStdError = fSigmaStdError(vThetaStar);
		return vSigmaStdError;
	}else{
		return 1;
	}
}

/*
**  Function:	Simulates and Estimates Garch data and parameters many times
**				to illustrate Asymptotic normality
**
**  Input: 		amMonteCarlo [matrix of many estimated parameters];
**
**  Output: 	1
*/

fMonteCarlo(const amMonteCarlo){
	decl mTemp;
	mTemp = zeros(iB,iPARS);

	for(decl i = 0; i<iB ; i++){
		decl vReturns;
		fSimALGAS(dALPHA, dBETA, dOMEGA, dGAMMA, dPROB, &vReturns, i);

		decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat, dP_hat;
		fEstimateGarch(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &dP_hat);	 //Omega and Gamma also estimated

		mTemp[i][0] =  dAlpha_hat - dALPHA;
		mTemp[i][1]	=  dBeta_hat - dBETA;
		mTemp[i][2] =  dOmega_hat - dOMEGA;
		mTemp[i][3]	=  dGamma_hat - dGAMMA;
		mTemp[i][4]	=  dP_hat - dPROB;
	}
	amMonteCarlo[0] = mTemp;
	return 1;
}

/*
**  Function:	Simulated and Estimates Garch data and parameters many times
**				to illustrate consistency it returns minimum, mean and maximum values for the estimated parameters
**
**  Input: 		amAlpha [matrix containing the min, max and mean of estimated alpha],
**				amBeta [matrix containing the min, max and mean of estimated beta], 
**				amOmega [matrix containing the min, max and mean of estimated omega],
**				amGamma [matrix containing the min, max and mean of estimated gamma]
**
**  Output: 	1
*/

fMonteCarlo2(const amAlpha, const amBeta, const amOmega, const amGamma, const amP, const amAlpha2, const amBeta2, const amOmega2, const amGamma2, const amP2){

	decl mTemp, mTempAlpha, mTempBeta, mTempOmega, mTempGamma, mTempP;
	decl mTemp2, mTempAlpha2, mTempBeta2, mTempOmega2, mTempGamma2, mTempP2;
	mTempAlpha = mTempBeta = mTempOmega = mTempGamma= mTempP = zeros((iSIZE/iSTEPS),3);
	mTempAlpha2 = mTempBeta2 = mTempOmega2 = mTempGamma2 = mTempP2 = zeros((iSIZE/iSTEPS),3);
	mTemp = mTemp2 =zeros(iB,iPARS);

	decl iSize = iSIZE;

	for(decl j = 0; j<(iSize/iSTEPS) ; j++){
		iSIZE = ((iSTEPS)*(j+1));
		for(decl i = 0; i<iB ; i++){
			decl vReturns;
			fSimALGAS(dALPHA, dBETA, dOMEGA, dGAMMA, dPROB, &vReturns, i);
	
			decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat, dP_hat, vSE;
			vSE = fEstimateGarch(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &dP_hat);	 //Omega and Gamma also estimated
			
			mTemp[i][0] =  sqrt(iSIZE)*(dAlpha_hat-dALPHA);				//SQRT(T)*(\hat_\alpha_T - \alpha_0) ~ N(0, \SIGMA)
			mTemp[i][1]	=  sqrt(iSIZE)*(dBeta_hat-dBETA);
			mTemp[i][2]	=  sqrt(iSIZE)*(dOmega_hat-dOMEGA);
			mTemp[i][3]	=  sqrt(iSIZE)*(dGamma_hat-dGAMMA);
			mTemp[i][4]	=  sqrt(iSIZE)*(dP_hat-dPROB);
			
			mTemp2[i][0] 	=  (dAlpha_hat-dALPHA)/vSE[0];				//(\hat_\alpha_T - \alpha_0)/SE(\hat_\alpha) ~ N(0, 1)
			mTemp2[i][1]	=  (dBeta_hat-dBETA)/vSE[1];
			mTemp2[i][2]	=  (dOmega_hat-dOMEGA)/vSE[2];
			mTemp2[i][3]	=  (dGamma_hat-dGAMMA)/vSE[3];
			mTemp2[i][4]	=  (dP_hat-dPROB)/vSE[4];

		}
		// v0.025_quantile, vMean, v0.975_quantile;				We get 95%-intervals
		mTempAlpha[j][0] = quantilec(mTemp[][],0.025)'[0];
		mTempAlpha[j][1] = meanc(mTemp[][])'[0];
		mTempAlpha[j][2] = quantilec(mTemp[][],0.975)'[0];
	
		mTempBeta[j][0] = quantilec(mTemp[][],0.025)'[1];
		mTempBeta[j][1] = meanc(mTemp[][])'[1];
		mTempBeta[j][2] = quantilec(mTemp[][],0.975)'[1];

		mTempOmega[j][0] = quantilec(mTemp[][],0.025)'[2];
		mTempOmega[j][1] = meanc(mTemp[][])'[2];
		mTempOmega[j][2] = quantilec(mTemp[][],0.975)'[2];
	
		mTempGamma[j][0] = quantilec(mTemp[][],0.025)'[3];
		mTempGamma[j][1] = meanc(mTemp[][])'[3];
		mTempGamma[j][2] = quantilec(mTemp[][],0.975)'[3];

		mTempP[j][0] = quantilec(mTemp[][],0.025)'[4];
		mTempP[j][1] = meanc(mTemp[][])'[4];
		mTempP[j][2] = quantilec(mTemp[][],0.975)'[4];
		

		mTempAlpha2[j][0] = quantilec(mTemp2[][],0.025)'[0];
		mTempAlpha2[j][1] = quantilec(mTemp2[][],0.5)'[0];	  //deletec()
		mTempAlpha2[j][2] = quantilec(mTemp2[][],0.975)'[0];
	
		mTempBeta2[j][0] = quantilec(mTemp2[][],0.025)'[1];
		mTempBeta2[j][1] = quantilec(mTemp2[][],0.5)'[1];
		mTempBeta2[j][2] = quantilec(mTemp2[][],0.975)'[1];

		mTempOmega2[j][0] = quantilec(mTemp2[][],0.025)'[2];
		mTempOmega2[j][1] = quantilec(mTemp2[][],0.5)'[2];
		mTempOmega2[j][2] = quantilec(mTemp2[][],0.975)'[2];
	
		mTempGamma2[j][0] = quantilec(mTemp2[][],0.025)'[3];
		mTempGamma2[j][1] = quantilec(mTemp2[][],0.5)'[3];
		mTempGamma2[j][2] = quantilec(mTemp2[][],0.975)'[3];

		mTempP2[j][0] = quantilec(mTemp2[][],0.025)'[4];
		mTempP2[j][1] = quantilec(mTemp2[][],0.5)'[4];
		mTempP2[j][2] = quantilec(mTemp2[][],0.975)'[4];
	}

	amAlpha[0] = mTempAlpha;
	amBeta[0] = mTempBeta;
	amOmega[0] = mTempOmega;
	amGamma[0] = mTempGamma;
	amP[0] = mTempP;

	amAlpha2[0] = mTempAlpha2;
	amBeta2[0] = mTempBeta2;
	amOmega2[0] = mTempOmega2;
	amGamma2[0] = mTempGamma2;
	amP2[0] = mTempP2;

	return 1;
}

/*
**				MAIN PROGRAM
**
**  Purpose:	Simulate GARCH returns for alpha, omega, gamma and beta many times.
**				Estimate GARCH parameters alpha, beta, omega and gamma.
**
**  Input: 		dALPHA, dBETA, dOMEGA, dGAMMA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/
main()
{
	//SET PARAMETERS
	dALPHA = 0.1;
 	dBETA = 0.99;			//later dus schatten
	dOMEGA = 0.01;
	dGAMMA = 0.1;
	dPROB = 0.45;
//	dGAMMA = dOMEGA/(1-dALPHA-dBETA); //doesn't work becomes negative
	iPARS = 5;


/*
** ..................................................................................	
**	 		ASYMPTOTIC NORMALITY
**	Get distributions of alpha and beta (to check for asymptotic normality)
**..................................................................................
*/

	//SET # OF SIMULATIONS 
	iB = 500; 			//max 5000
	iSIZE = 500;		//max 5000
	iSIMS = iB*iSIZE;
	//vSTD_ALD = rann(iSIMS,1);
	decl vRanExp, vRanExp2, vRanALD;
	vRanExp = ranexp(iSIMS, 1, 1);
	vRanExp2 = ranexp(iSIMS, 1, 1);
	vRanALD = vRanExp/dPROB-vRanExp2/(1-dPROB);	 //~ALD(0,1,p)
	vSTD_ALD = (1-dPROB)*dPROB/(sqrt(1-2*dPROB+2*sqr(dPROB)))*vRanALD;
	bSTD_ERROR = FALSE;				 //boolean

	//GET BETA SUCH THAT WE GET EQUALITY "Elog(alpha_0 z_t^2 + beta_0) = 0".
	decl dBeta_0;
	fGetBeta(iSIMS,  &dBeta_0);
	dBETA = dBeta_0;

	
	//DO MANY SIMULATIONS AND ESITMATIONS	
	decl mMonteCarlo;
	fMonteCarlo(&mMonteCarlo);	  

	//DRAW GRAPHS
	SetDrawWindow("SIM_ald-GAS_2_AsymN");
	DrawDensity(0, (mMonteCarlo[][0])', {"(i) Density alpha"});
	DrawDensity(1, (mMonteCarlo[][1])', {"(ii) Density beta"});
	DrawDensity(2, (mMonteCarlo[][2])', {"(iii) Density omega"});
	DrawDensity(3, (mMonteCarlo[][3])', {"(iv) Density gamma"});
	DrawDensity(4, (mMonteCarlo[][4])', {"(v) Density p"});
	ShowDrawWindow();

	print("\nFirst Graph Finished at ",time(),"\n");
/*
** ..................................................................................	
**	 			CONSISTENCY
**	Check consistency for alpha and beta
** ..................................................................................
*/	

	//SET # OF SIMULATIONS 
	iB = 50;			//100
	iSIZE = 10000;		//10000
	iSIMS = iB*iSIZE;
	vRanExp = ranexp(iSIMS, 1, 1);
	vRanExp2 = ranexp(iSIMS, 1, 1);
	vRanALD = vRanExp/dPROB-vRanExp2/(1-dPROB);	 //~ALD(0,1,p)
	vSTD_ALD = (1-dPROB)*dPROB/(sqrt(1-2*dPROB+2*sqr(dPROB)))*vRanALD;
	bSTD_ERROR = TRUE;

//	//GET BETA SUCH THAT WE GET EQUALITY "Elog(alpha_0 z_t^2 + beta_0) = 0".
	fGetBeta(iSIMS,  &dBeta_0);
	dBETA = dBeta_0;
	
	//DO MANY SIMULATIONS AND ESITMATIONS
	decl mAlpha, mBeta, mOmega, mGamma, mP, mAlpha2, mBeta2, mOmega2, mGamma2, mP2;
	iSTEPS = iSIZE/10;				 	//steps of iSIZE/100 takes a while (steps of iSIZE/10 is faster)
	fMonteCarlo2(&mAlpha, &mBeta, &mOmega, &mGamma, &mP, &mAlpha2, &mBeta2, &mOmega2, &mGamma2, &mP2);

	//DRAW GRAPHS
	SetDrawWindow("SIM_ald-GAS_2_Cons");
	Draw(0, (mAlpha)',iSTEPS,iSTEPS);
	Draw(1, (mBeta)',iSTEPS,iSTEPS);
	Draw(2, (mOmega)',iSTEPS,iSTEPS);
	Draw(3, (mGamma)',iSTEPS,iSTEPS);
	Draw(4, (mP)',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) alpha");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) omega");	
	DrawTitle(3,"(iv) gamma");
	DrawTitle(4,"(v) p");
	ShowDrawWindow();
	print("\nSecond Graph Finished at ",time(),"\n");

	SetDrawWindow("SIM_ald-GAS_2_NormC");
	Draw(0, (mAlpha2)',iSTEPS,iSTEPS);
	Draw(1, (mBeta2)',iSTEPS,iSTEPS);
	Draw(2, (mOmega2)',iSTEPS,iSTEPS);
	Draw(3, (mGamma2)',iSTEPS,iSTEPS);
	Draw(4, (mP2)',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) alpha");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) omega");	
	DrawTitle(3,"(iv) gamma");
	DrawTitle(4,"(v) p");
	ShowDrawWindow();
	print("\nThird Graph Finished at ",time(),"\n");
}

