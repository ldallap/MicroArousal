#ifndef Channels_ECells_H
#define Channels_ECells_H

#include<cmath>
#include "Variables.h"

/***********************************************************/
/*         FAST SODIUM CHANNEL FOR E-CELLS                 */
/***********************************************************/
double Na_steadystate_m_E(double memV){
  double alpham=4.*0.1*(memV+33.)/(1.-exp(-0.1*(memV+33.)));;
  double betam=4.*4.*exp(-0.083333333333*(memV+53.7));
  return alpham/(alpham+betam);
}
double Na_inactgate_h_E(double variable,double memV){
  double alphah=4.*0.07*exp(-0.1*(memV+50.));;
  double betah=4.*1./(1.+exp(-0.1*(memV+20.)));
  return (alphah*(1.0-variable))-(betah*variable);
}
/**/
/***************************************************************/
/*          DELAYED RECTIFIER - K CHANNEL- FOR E-CELLS         */
/***************************************************************/
double K_actgate_n_E(double variable, double memV){
  double betan=4.*0.125*exp(-0.04*(memV+44.));
  double alphan=4*0.01*(memV+34)/(1.-exp(-0.1*(memV+34.)));
  return (alphan*(1.0-variable))-(betan*variable);
}
/**/
/**********************************************************/
/*                A-TYPE K CHANNEL FOR E-CELLS            */
/**********************************************************/
double A_steadystate_m_E(double memV){
  return 1./(1.+exp(-(memV+50.)/20.));
}
double A_inactgate_h_E(double variable, double memV){
  double inactivOpenRate=1./(1.+exp((memV+80.)/6.))/15.;
  double inactivCloseRate=1./15.-1./(1.+exp((memV+80.)/6.))/15.;

  return (inactivOpenRate+inactivCloseRate)*((inactivOpenRate/(inactivOpenRate+inactivCloseRate))-variable);
}
/**/
/**********************************************************/
/* SLOWLY INACTIVATING A-TYPE K CHANNEL FOR E-CELLS - KS  */
/**********************************************************/
double KS_actgate_m_E(double variable, double memV){
  double minf=1./(1.+exp(-(memV+34.)/6.5));
  double tau=8./(exp(-(memV+55.)/30.)+exp((memV+55.)/30.));
  return (minf-variable)/tau;
}
/**/
/**********************************************************/
/*   HIGH THRESHOLD CALCIUM CHANNEL FOR E-CELLS  - Ca     */
/**********************************************************/
double Ca_act_m_E(double memV){
  return 1./(1.+exp(-0.11111111*(memV+20.)));
}
/**/
/**********************************************************/
/*          ANOMALOUS RECTIFIER K CHANNEL FOR E-CELLS     */
/**********************************************************/
double AR_act_h_E(double memV){
  return 1./(1.+exp((memV+75.0)*0.25));
}
/**/
/**********************************************************/
/*          PERSISTENT NA CHANNEL FOR E-CELLS - NaP       */
/**********************************************************/
double NaP_act_m_E(double memV){
  return 1./(1.+exp(-(memV+55.7)/7.7));
}
/**/
/**********************************************************/
/*         NA DEPENDENT K  CHANNEL FOR E-CELLS * KNa       */
/**********************************************************/
double KNa_winf_E(double variable){
  double fac=pow(38.7/variable,3.5);
  return 0.37/(1.+fac);
}
/**/



/**********************************************************/
/*                The Dynamics start here                 */
/**********************************************************/
double syn_sigmoid(double memV){
  return 1./(1.+exp(-0.5*(memV-20.)));
}
void allderivsE(double *y, double *dy, double IsynSoma, double IsynDend, int index, double gKNa_cond, double gKCa_cond,  int keyKNa, int keyKCa){

	const double phi_eq = Rpump*naEq*naEq*naEq/(naEq*naEq*naEq+kp3);

  double NaChannel=0.0;
	NaChannel=pow(Na_steadystate_m_E(y[0]),3)*y[1]*gNa_E*0.5*0.3*(y[0]-55.0);

	double KChannel=0.0;
	KChannel=y[2]*y[2]*y[2]*y[2]*gK_E*0.5*0.3*(y[0]+100.0);

	double LeakChannelS=0.0;
	LeakChannelS=glE[index]*0.3*(y[0]-VlE[index]);

  double LeakChannelD=0.0;
	LeakChannelD=glE[index]*0.7*(y[5]-VlE[index]);

  double AChannel=0.0;
  AChannel=A_steadystate_m_E(y[0])*A_steadystate_m_E(y[0])*A_steadystate_m_E(y[0]);
  AChannel*=gA*0.5*0.3*y[3]*(y[0]+100.0);

  double KSChannel=0.0;
  KSChannel=gKS*0.5*0.3*y[4]*(y[0]+100.0);

  /*Dendrite below - Soma above*/

  double CaChannel=0.0;
  CaChannel=gCa*0.5*0.7*Ca_act_m_E(y[5])*Ca_act_m_E(y[5])*(y[5]-120.0);

  double NaPChannel=0.0;
  NaPChannel=NaP_act_m_E(y[5])*NaP_act_m_E(y[5])*NaP_act_m_E(y[5]);
  NaPChannel*=gNaP*0.5*0.7*(y[5]-55.0);

  double ARChannel=0.0;
  ARChannel=gAR*0.5*0.7*AR_act_h_E(y[5])*(y[5]+100.0);

  double KCaChannel=0.0;
  KCaChannel=gKCa_cond*0.5*0.7*y[6]/(y[6]+30.)*(y[5]+100.0);
  KCaChannel=KCaChannel*keyKCa;

  double KNaChannel=0.0;
  KNaChannel=gKNa_cond*0.5*0.3*KNa_winf_E(y[7])*(y[0]+100.0);
  KNaChannel=KNaChannel*keyKNa;

  // if(keyKCa==1) cout << "WITH " << KCaChannel << endl;
  // if(keyKCa==0) cout << "WITHOUT " << KCaChannel << endl;

  const double iCsoma = 1./0.5/0.3;
  const double iCdend = 1./0.5/0.7;

  double iSoma = -LeakChannelS-NaChannel-KChannel-AChannel-KSChannel-KNaChannel-gsd[index]*(y[0]-y[5])+IsynSoma;

  double iDend = -LeakChannelD-CaChannel-NaPChannel-ARChannel-KCaChannel-gsd[index]*(y[5]-y[0])+IsynDend;

  dy[0] = iSoma*iCsoma;
  dy[1] = Na_inactgate_h_E(y[1],y[0]);
  dy[2] = K_actgate_n_E(y[2],y[0]);
  dy[3] = A_inactgate_h_E(y[3],y[0]);
  dy[4] = KS_actgate_m_E(y[4],y[0]);

  dy[5] = iDend*iCdend;

  dy[6] = -alphaCa*CaChannel-y[6]/tauCa;
  double na3=y[7]*y[7]*y[7];
  dy[7] = -naAlpha*(NaChannel+NaPChannel)-Rpump*na3/(na3+kp3) + phi_eq;

  dy[8] = 3.48*syn_sigmoid(y[0]) - y[8]/2.0;
  dy[9] = 0.5*y[10]*(1.-y[9])-y[9]/100.;
  dy[10] = 3.48*syn_sigmoid(y[0]) - y[10]/2.;
}

/**/
void rk4E(double curr, double cur, int ind, double gKNa_cond, double gKCa_cond,  int keyKNa, int keyKCa){
	int i;
	double hh=dt*0.5;
	double h6=dt*0.166666666666666667;
	for (i=0;i<nVarsE;i++)
		yt[i]=varsE[i]+hh*dVarsE[i];
	allderivsE(yt,dyt,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa);
	for (i=0;i<nVarsE;i++)
		yt[i]=varsE[i]+hh*dyt[i];
	allderivsE(yt,dym,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa);
	for (i=0;i<nVarsE;i++) {
		yt[i]=varsE[i]+dt*dym[i];
		dym[i] += dyt[i];
	}
	allderivsE(yt,dyt,curr,cur,ind, gKNa_cond, gKCa_cond, keyKNa, keyKCa);
	for (i=0;i<nVarsE;i++)
		varsE[i] += h6*(dVarsE[i]+dyt[i]+2.0*dym[i]);
}
/**/
void rungeKutta4E(double timesim, int index, double currExt1, double currExt2, double gKNa_cond, double gKCa_cond, int keyKNa, int keyKCa){
	varsE[0]=y[0];
	varsE[1]=y[1];
	varsE[2]=y[2];
  varsE[3]=y[3];
  varsE[4]=y[4];
  varsE[5]=y[5];
  varsE[6]=y[6];
  varsE[7]=y[7];
  varsE[8]=y[8];
  varsE[9]=y[9];
  varsE[10]=y[10];
	allderivsE(varsE,dVarsE,currExt1,currExt2,index, gKNa_cond, gKCa_cond, keyKNa, keyKCa);
  rk4E(currExt1,currExt2,index, gKNa_cond, gKCa_cond, keyKNa, keyKCa);
	y[0]=varsE[0];
	y[1]=varsE[1];
	y[2]=varsE[2];
  y[3]=varsE[3];
  y[4]=varsE[4];
  y[5]=varsE[5];
  y[6]=varsE[6];
  y[7]=varsE[7];
  y[8]=varsE[8];
  y[9]=varsE[9];
  y[10]=varsE[10];
}
/**/





#endif
