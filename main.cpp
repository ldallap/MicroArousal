#include<iostream>
#include "Randomnumber.hpp"
#include "Channels_ECells.hpp"
#include "Channels_ICells.hpp"
#include "Variables.hpp"
#include "Network.hpp"
#include "NeuronsStructure.hpp"

using namespace std;




int main (int argc, char *argv[]){

  double gKNa_conductance=0.0135; /*1.33*/
  double gKCa_conductance=0.70; /*0.57*/
  double block_gKCa=0.0;

  /*original*/
  // double gKNa_conductance=0.0135; /*1.33*/
  // double gKCa_conductance=0.70; /*0.57*/
  // double block_gKCa=0.0;
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  if(argc != 3){
    cout << "#########-ERROR-########" << endl;
    cout << "# Give 2 paramenters as input:" << endl;
    cout << "# 1) Block gKCa" << endl;
    cout << "# 2) Seed" << endl;
    cout << "########################" << endl;
    return 1;
  }
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  block_gKCa= (double) atof(argv[1]);
  gKCa_conductance=gKCa_conductance*block_gKCa;
  
  seed=(long) atof(argv[2]);
  long seed1=seed;

  startpointers();
  construct_network();



  /*SAVING OUTFILE*/
  FILE *outf, *outf2;
  char *Outfile;
  char *Outfile2;
  Outfile=new char[500];
  Outfile2=new char[500];
  sprintf(Outfile,"Raster_Exc_gKCa%.3lf_blockgKCa%.3lf_seed%ld.dat",gKCa_conductance,block_gKCa,seed);
  sprintf(Outfile2,"Raster_Inh_gKCa%.3lf_blockgKCa%.3lf_seed%ld.dat",gKCa_conductance,block_gKCa,seed);

  outf = fopen(Outfile,"w");
  outf2 = fopen(Outfile2,"w");

  int *electrode=new int [5];
  electrode[0]=200+int(ran(&seed)*100);
  electrode[1]=400+int(ran(&seed)*100);
  electrode[2]=600+int(ran(&seed)*100);
  electrode[3]=800+int(ran(&seed)*100);
  electrode[4]=1000+int(ran(&seed)*100);

// timesimulation
  for(int t=0;t<700000;t++){
    double currtimesimu=t*dt*0.001; /*current time simulation*/
    double auxsumLFP=0.0;
    double auxsumLFP1=0.0;
    double auxsumLFP2=0.0;
    double auxsumLFP3=0.0;
    double auxsumLFP4=0.0;

    if(currtimesimu>=40) return 1;

    ComputeSynpases();

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /*Excitatory looping over time - EXC DYNAMICS*/
    for(int j=0;j<numEcells;j++){
      double eCondAmpa=gEE_A*EESynInputAMPA[j];
      double eCondNMDA=gEE_N*EESynInputNMDA[j]/(1.0 + 0.2801*exp(-0.062*memV_DendE[j]));
      double iCondGABA=gIE*IESynInput[j];
      double currtodend=(eCondAmpa+eCondNMDA)*(vSynGlu-memV_DendE[j]);
      double currtosoma=iCondGABA*(vSynGABA-memV_SomaE[j]);

      y[0]=memV_SomaE[j];y[1]=h_Na_E[j];
      y[2]=n_K_E[j];y[3]=h_A_E[j];
      y[4]=m_KS_E[j];y[5]=memV_DendE[j];
      y[6]=conc_Ca_E[j];y[7]=conc_Na_E[j];
      y[8]= sAmpaE[j];y[9]= sNmdaE[j];
      y[10]= xNmdaE[j];
      for(int z=0;z<nVarsE;z++){
        dVarsE[z]=aux_dvarsE[j][z];
        dym[z]=aux_dym[j][z];
        dyt[z]=aux_dyt[j][z];
        yt[z]=aux_yt[j][z];
      }


      if(j<=electrode[0]+25 and j>=electrode[0]-25) auxsumLFP+=abs(currtosoma)+abs(currtodend);
      if(j<=electrode[1]+25 and j>=electrode[1]-25) auxsumLFP1+=abs(currtosoma)+abs(currtodend);
      if(j<=electrode[2]+25 and j>=electrode[2]-25) auxsumLFP2+=abs(currtosoma)+abs(currtodend);
      if(j<=electrode[3]+25 and j>=electrode[3]-25) auxsumLFP3+=abs(currtosoma)+abs(currtodend);
      if(j<=electrode[4]+25 and j>=electrode[4]-25) auxsumLFP4+=abs(currtosoma)+abs(currtodend);

      rungeKutta4E(currtimesimu,j,currtosoma,currtodend, gKNa_conductance, gKCa_conductance, 1., 1.);



      // if(j==585) cout << currtimesimu << "\t" << memV_SomaE[j] << "\t" << gKCa_conductance*0.5*0.7*y[6]/(y[6]+30.) << "\t" << y[6];
      // if(j==604) cout << "\t" << memV_SomaE[j] << "\t" << gKCa_conductance*0.5*0.7*y[6]/(y[6]+30.) << "\t" << y[6];
      // if(j==606) cout << "\t" << memV_SomaE[j] << "\t" << gKCa_conductance*0.5*0.7*y[6]/(y[6]+30.) << "\t" << y[6];
      if(((memV_SomaE[j]-vth)*(y[0]-vth)<0.0) && (y[0]>vth) && t>10000){
        fprintf(outf,"%lf\t%d\n", currtimesimu-(10000*0.001*dt), j);
      }

      memV_SomaE[j]=y[0];h_Na_E[j]=y[1];
      n_K_E[j]=y[2];h_A_E[j]=y[3];
      m_KS_E[j]=y[4];memV_DendE[j]=y[5];
      conc_Ca_E[j]=y[6];conc_Na_E[j]=y[7];
      sAmpaE[j]=y[8];sNmdaE[j]=y[9];
      xNmdaE[j]=y[10];
      for(int z=0;z<nVarsE;z++){
        aux_dvarsE[j][z]=dVarsE[z];
        aux_dym[j][z]=dym[z];
        aux_dyt[j][z]=dyt[z];
        aux_yt[j][z]=yt[z];
      }
    }/*End looping EXC*/
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

    if(t%10==0) cout << currtimesimu << "\t" << auxsumLFP << "\t" << auxsumLFP1 << "\t" << auxsumLFP2 << "\t" << auxsumLFP3 << "\t" << auxsumLFP4 << endl;



    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /*Inhibitory looping over time - EXC DYNAMICS*/
    for(int j=0;j<numIcells;j++){
      double eCondAmpa=gEI_A*EISynInputAMPA[j];
      double eCondNMDA=gEI_N*EISynInputNMDA[j]/(1.0 + 0.2801*exp(-0.062*memV_SomaI[j]));
      double iCondGABA=gII*IISynInput[j];
      double currtodend=(eCondAmpa+eCondNMDA)*(vSynGlu-memV_SomaI[j]);
      double currtosoma=iCondGABA*(vSynGABA-memV_SomaI[j]);

      x[0]=memV_SomaI[j];
      x[1]=h_Na_I[j];
      x[2]=n_K_I[j];
      x[3]=sGabaI[j];
      for(int z=0;z<nVarsI;z++){
        aux_dvarsI[j][z]=dVarsI[z];
        aux_dxm[j][z]=dxm[z];
        aux_dxt[j][z]=dxt[z];
        aux_xt[j][z]=xt[z];
      }

      rungeKutta4I(currtimesimu,j, currtodend+currtosoma);

      // if(j==145) cout << "\t" << memV_SomaI[j];
      // if(j==150) cout << "\t" << memV_SomaI[j] ;
      // if(j==155) cout << "\t" << memV_SomaI[j] << endl;
      // if(((memV_SomaI[j]-vth)*(x[0]-vth)<0.0) && t>10000){
      //   fprintf(outf2,"%lf\t%d\n", currtimesimu-(10000*0.001*dt), j);
      // }

      memV_SomaI[j]=x[0];h_Na_I[j]=x[1];
      n_K_I[j]=x[2];sGabaI[j]=x[3];
      for(int z=0;z<nVarsI;z++){
        aux_dvarsI[j][z]=dVarsI[z];
        aux_dxm[j][z]=dxm[z];
        aux_dxt[j][z]=dxt[z];
        aux_xt[j][z]=xt[z];
      }
    }/*End Looping inh*/
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

  }/*END LOOP TIME*/
  fclose(outf);
  fclose(outf2);

//  /*SAVING OUTFILE*/
//  FILE *outf, *outf2, *outf3,*outf5;
// char *Outfile;
//  char *Outfile2;
//  char *Outfile3;
// Outfile=new char[500];
//  Outfile2=new char[500];
//  Outfile3=new char[500];
//  sprintf(Outfile,"Raster_ConcentrationKNa_%.2lf_gKNa%.3lf_ConcentrationKCa_%.2lf_gKCa%.3lf_seed%ld.ts",KNa_percentage,gKNa_conductance,KCa_percentage,gKCa_conductance,seed1);
//  sprintf(Outfile2,"MembranePotential_IonChannels_ConcentrationKNa_%.2lf_gKNa%.3lf_ConcentrationKCa_%.2lf_gKCa%.3lf_seed%ld_keyChannel%d.ts",KNa_percentage,gKNa_conductance,KCa_percentage,gKCa_conductance,seed1,keychannel);
//  sprintf(Outfile3,"LFP_IonChannels_ConcentrationKNa_%.2lf_gKNa%.3lf_ConcentrationKCa_%.2lf_gKCa%.3lf_seed%ld_keyChannel%d.ts",KNa_percentage,gKNa_conductance,KCa_percentage,gKCa_conductance,seed1,keychannel);


//  outf = fopen(Outfile,"w");
// for(int t=0;t<timesimulation;t++){
//    if(timeindexes[t]==-1 or spikesindexes[t]==-1) break;
//    else fprintf(outf,"%lf\t%d\n", timeindexes[t], spikesindexes[t]);
//  }
//  fclose(outf);



  deletepointers();

  return 0;
}
