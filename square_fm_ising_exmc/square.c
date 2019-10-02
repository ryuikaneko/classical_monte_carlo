/*-------------------------------------------------------------------*
 * exchange monte carlo
 *-------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../dSFMT/dSFMT.h"

#define min(a,b) (((a)<(b)) ? (a):(b))
dsfmt_t dsfmt;

/*-------------------------------------------------------------------*
 * get time
 *-------------------------------------------------------------------*/
double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

/*-------------------------------------------------------------------*
 * energy = -J sum_i,j spins[i][j] * (spins[i+1][j] + spins[i][j+1])
 *-------------------------------------------------------------------*/
int calc_energy
(
  int n,
  int m,
  double J,
  int ***spins,
  double *ene
)
{
  int i,ipp,j,jpp;
  double tmp_ene;
  tmp_ene = 0.0;
  for(i=0; i<n; i++){
    ipp = (i==n-1) ? 0 : i+1;
    for(j=0; j<n; j++){
      jpp = (j==n-1) ? 0 : j+1;
      tmp_ene -= J * spins[m][i][j]
        * (spins[m][ipp][j] + spins[m][i][jpp]);
    }
  }
  ene[0] = tmp_ene;
  return 0;
}

/*-------------------------------------------------------------------*
 * magnetization = sum_i,j spins[i][j]
 *-------------------------------------------------------------------*/
double calc_mag
(
  int n,
  int m,
  int ***spins,
  double *mag
)
{
  int i,j;
  double tmp_mag;
  tmp_mag = 0.0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      tmp_mag += spins[m][i][j];
    }
  }
  mag[0] = tmp_mag;
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step
 *-------------------------------------------------------------------*/
int one_flip
(
  int n,
  int m,
  double **BoltzFac,
  int *list_temp,
  int ***spins
)
{
  int i,j,ipp,imm,jpp,jmm;
  int egap;
  /* spin flip candidate */
  i = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ipp = (i==n-1) ? 0 : i+1;
  imm = (i==0) ? n-1 : i-1;
  jpp = (j==n-1) ? 0 : j+1;
  jmm = (j==0) ? n-1 : j-1;
  egap = 2 * spins[m][i][j]
    *(spins[m][ipp][j] + spins[m][i][jpp]
    + spins[m][imm][j] + spins[m][i][jmm]);
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][egap+12]){
    spins[m][i][j] *= -1;
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step
 * (calculate energy and magnetization differences)
 *-------------------------------------------------------------------*/
int one_flip_ene_mag
(
  int n,
  int m,
  double J,
  double **BoltzFac,
  int *list_temp,
  int ***spins,
  double *ene,
  double *mag
)
{
  int i,j,ipp,imm,jpp,jmm;
  int egap;
  /* spin flip candidate */
  i = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ipp = (i==n-1) ? 0 : i+1;
  imm = (i==0) ? n-1 : i-1;
  jpp = (j==n-1) ? 0 : j+1;
  jmm = (j==0) ? n-1 : j-1;
  egap = 2 * spins[m][i][j]
    *(spins[m][ipp][j] + spins[m][i][jpp]
    + spins[m][imm][j] + spins[m][i][jmm]);
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][egap+12]){
    spins[m][i][j] *= -1;
    ene[m] += J*egap;
    mag[m] += 2 * spins[m][i][j];
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * temperature exchange
 *-------------------------------------------------------------------*/
int exchange_states
(
  int n,
  int m,
  double J,
  double *T,
  int *list_temp,
  int *list_conf,
  int ***spins
)
{
  int conf0,conf1;
  double e0,e1;
  double egap,w;
  /* exclude edge temperature */
  /* calculate gap */
  conf0 = list_conf[m];
  conf1 = list_conf[m+1];
  calc_energy(n,conf0,J,spins,&e0);
  calc_energy(n,conf1,J,spins,&e1);
  egap = (1.0/T[m+1] - 1.0/T[m]) * (e0 - e1);
  /* metropolis method */
  w = exp(-egap);
  if(dsfmt_genrand_close_open(&dsfmt) <= w){
    list_temp[conf0] = m+1;
    list_temp[conf1] = m;
    list_conf[m]     = conf1;
    list_conf[m+1]   = conf0;
  }
  return 0;
}

int exchange_states_ene_mag
(
  int n,
  int m,
  double J,
  double *T,
  int *list_temp,
  int *list_conf,
  int ***spins,
  double *ene,
  double *mag
)
{
  int conf0,conf1;
  double egap,w;
  /* calculate gap */
  conf0 = list_conf[m];
  conf1 = list_conf[m+1];
  egap = (1.0/T[m+1] - 1.0/T[m]) * (ene[conf0] - ene[conf1]);
  /* metropolis method */
  w = exp(-egap);
  if(dsfmt_genrand_close_open(&dsfmt) <= w){
    list_temp[conf0] = m+1;
    list_temp[conf1] = m;
    list_conf[m]     = conf1;
    list_conf[m+1]   = conf0;
  }
  return 0;
}

int main
(
  int argc,
  char *argv[]
)
{
  int i,j,m,l;
  char opt;
  int N = 8;
  int seed = N * 10000;
  int NN = N*N;
  int N_t = 61;/* # of temperature */
  int N_equ = 10000;
  int N_smp = 10000;
  int N_cor = 10000;
  int N_blk = 5;
  double T_i = 2.10;
  double T_f = 2.405;
  double J = +1.0;/* ferro */
  double timei,timef;

  /* prameter settings */
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 'n'){
        N = atoi(argv[i+1]);
        seed = N * 10000;
        NN = N*N;
      }else if(opt == 'e'){
        N_equ = atoi(argv[i+1]);
      }else if(opt == 's'){
        N_smp = atoi(argv[i+1]);
      }else if(opt == 'c'){
        N_cor = atoi(argv[i+1]);
      }else if(opt == 'b'){
        N_blk = atoi(argv[i+1]);
      }else if(opt == 'i'){
        T_i = atof(argv[i+1]);
      }else if(opt == 'f'){
        T_f = atof(argv[i+1]);
      }else if(opt == 't'){
        N_t = atoi(argv[i+1]);
      }else if(opt == 'j'){
        J = atof(argv[i+1]);
      }else if(opt == 'r'){
        seed = atoi(argv[i+1]);
      }
    }
  }

  double T_width = (T_f - T_i)/(double)N_t;
  double Etemp,Mtemp;
  double Ene[N_t],Ene2[N_t];
  double Mag[N_t],Mag2[N_t],Mag4[N_t];
  double C[N_t],X[N_t],BP[N_t];
  double absMag[N_t],absX[N_t];
  double Ea[N_t],Ma[N_t],Ca[N_t],Xa[N_t],Ba[N_t];
  double Ev[N_t],Mv[N_t],Cv[N_t],Xv[N_t],Bv[N_t];
  double Ee[N_t],Me[N_t],Ce[N_t],Xe[N_t],Be[N_t];
  double absMa[N_t],absXa[N_t],Mag2a[N_t];
  double absMv[N_t],absXv[N_t],Mag2v[N_t];
  double absMe[N_t],absXe[N_t],Mag2e[N_t];
  double Eadd[N_t],Madd[N_t];
  int list_temp[N_t];
  int list_conf[N_t];
  double T[N_t];

  int ***spins;
  spins = (int***)malloc(N_t*sizeof(int**));
  spins[0] = (int**)malloc(N_t*N*sizeof(int*));
  spins[0][0] = (int*)malloc(N_t*N*N*sizeof(int));
  for(m=0; m<N_t; m++){
    spins[m] = spins[0] + m*N;
    for(i=0; i<N; i++){
      spins[m][i] = spins[0][0] + (m*N + i)*N;
    }
  }
  double **BoltzFac;
  BoltzFac = (double**)malloc(N_t*sizeof(double*));
  BoltzFac[0] = (double*)malloc(N_t*25*sizeof(double));
  for(m=0; m<N_t; m++){
    BoltzFac[m] = BoltzFac[0] + m*25;
  }
  for(m=0; m<N_t; m++){
    list_temp[m] = m;
    list_conf[m] = m;
    T[m] = T_i + T_width * m;
  }

  /* table of Boltzmann factors */
  for(m=0; m<N_t; m++){
    for(i=-12; i<=12; i++){
      BoltzFac[m][i+12] = exp(-J*i/T[m]);
    }
  }
  /* end of table of Boltzmann factors */

  printf("# N=%d NN=%d N_equ=%d N_smp=%d N_cor=%d N_blk=%d\n",
    N,NN,N_equ,N_smp,N_cor,N_blk);
  printf("# seed=%d T_i=%13.6e T_f=%13.6e T_width=%13.6e J=%13.6e\n",
    seed,T_i,T_f,T_width,J);
  printf("# T E dE C dC M dM X dX BP dBP Mag2 dMag2 absM dabsM absX dabsX\n");

  timei = gettimeofday_sec();

  /* initialization */
  dsfmt_init_gen_rand(&dsfmt,seed);
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      spins[N_t-1][i][j] = (dsfmt_genrand_close_open(&dsfmt) > 0.5) ? +1 : -1;
//      spins[N_t-1][i][j] = 1;
    }
  }
  for(i=0; i<N_equ; i++){
    for(j=0; j<NN; j++){
      one_flip(N,N_t-1,BoltzFac,list_temp,spins);
    }
  }
  for(m=N_t-2; m>=0; m--){
    /* m: parameter of environment or temperature */
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
        spins[m][i][j] = spins[m+1][i][j];
      }
    }
    for(i=0; i<N_equ; i++){
      for(j=0; j<NN; j++){
        one_flip(N,m,BoltzFac,list_temp,spins);
      }
    }
  }

  /* initialization */
/*
  for(m=0; m<N_t; m++){// loop for each temperature
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
        spins[m][i][j] = (dsfmt_genrand_close_open(&dsfmt) > 0.5) ? +1 : -1;
      }
    }
  }
  // loop until equilibrium is achieved
  for(i=0; i<N_equ; i++){
    for(m=0; m<N_t; m++){// loop for each temperature
      for(j=0; j<NN; j++){
        one_flip(N,m,BoltzFac,list_temp,spins);
      }
    }
  }
*/



  /* parameter initialization */
  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of environment or temperature */
    Ea[m] = 0.0;  Ma[m] = 0.0;  Ca[m] = 0.0;  Xa[m] = 0.0;  Ba[m] = 0.0;
    Ev[m] = 0.0;  Mv[m] = 0.0;  Cv[m] = 0.0;  Xv[m] = 0.0;  Bv[m] = 0.0;
    Ee[m] = 0.0;  Me[m] = 0.0;  Ce[m] = 0.0;  Xe[m] = 0.0;  Be[m] = 0.0;
    Mag2a[m] = 0.0;  Mag2v[m] = 0.0;  Mag2e[m] = 0.0;
    absMa[m] = 0.0;  absXa[m] = 0.0;
    absMv[m] = 0.0;  absXv[m] = 0.0;
    absMe[m] = 0.0;  absXe[m] = 0.0;
  }
  /* end of parameter initialization */

  /* per each block */
  for(i=0; i<N_blk; i++){
    /* neglecting and exchanging */
    for(j=0; j<N_cor; j++){
      for(m=0; m<N_t-1; m++){/* loop for each temperature */
        /* m: parameter of environment */
        exchange_states(N,m,J,T,list_temp,list_conf,spins);
      }
      for(m=0; m<N_t; m++){/* loop for each temperature */
        /* m: parameter of environment */
        for(l=0; l<NN; l++){
          one_flip(N,m,BoltzFac,list_temp,spins);
        }
      }
    }
    /* end of neglecting and exchanging */

    /* parameter initialization */
    Etemp = 0.0;
    Mtemp = 0.0;
    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of environment */
      calc_energy(N,m,J,spins,&(Eadd[m]));
      calc_mag(N,m,spins,&(Madd[m]));
      Ene[m]  = 0.0;
      Ene2[m] = 0.0;
      Mag[m]  = 0.0;
      Mag2[m] = 0.0;
      Mag4[m] = 0.0;
      C[m]    = 0.0;
      X[m]    = 0.0;
      BP[m]   = 0.0;
      absMag[m] = 0.0;
      absX[m]   = 0.0;
    }
    /* end of parameter initialization */

    /* measurement loop */
    for(j=0; j<N_smp; j++){
      for(m=0; m<N_t-1; m++){/* loop for each temperature */
        /* m: parameter of environment */
        exchange_states_ene_mag
          (N,m,J,T,list_temp,list_conf,spins,Eadd,Madd);
      }
      for(m=0; m<N_t; m++){/* loop for each temperature */
        /* m: parameter of environment */
        for(l=0; l<NN; l++){
          one_flip_ene_mag
            (N,m,J,BoltzFac,list_temp,spins,Eadd,Madd);
        }
      }
      for(m=0; m<N_t; m++){/* loop for each temperature */
        /* m:            parameter of environment */
        /* list_temp[m]: parameter of temperature */
        Etemp    = Eadd[m];/* E */
        Ene[list_temp[m]]  += Etemp;
        Etemp    = Etemp*Etemp;/* E^2 */
        Ene2[list_temp[m]] += Etemp;
        Mtemp    = Madd[m];/* M */
        Mag[list_temp[m]]  += Mtemp;
        Mtemp    = Mtemp*Mtemp;/* M^2 */
        Mag2[list_temp[m]] += Mtemp;
        Mtemp    = Mtemp*Mtemp;/* M^4 */
        Mag4[list_temp[m]] += Mtemp;
        Mtemp    = fabs(Madd[m]);/* |M| */
        absMag[list_temp[m]]  += Mtemp;
      }
    }
    /* end of measurement loop */

    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of temperature */
      /*                 ^^^^^^^^^^^ */
      Ene[m]  /= (double)N_smp;
      Ene2[m] /= (double)N_smp;
      Mag[m]  /= (double)N_smp;
      Mag2[m] /= (double)N_smp;
      Mag4[m] /= (double)N_smp;
      C[m]  = (Ene2[m] - Ene[m]*Ene[m]) /T[m]/T[m];
      X[m]  = (Mag2[m] - Mag[m]*Mag[m]) /T[m];
      BP[m] = 1.0 - Mag4[m]/Mag2[m]/Mag2[m]/3.0;
      absMag[m]  /= (double)N_smp;
      absX[m]  = (Mag2[m] - absMag[m]*absMag[m]) /T[m];

      Ea[m] += Ene[m];
      Ev[m] += Ene[m]*Ene[m];
      Ma[m] += Mag[m];
      Mv[m] += Mag[m]*Mag[m];
      Ca[m] += C[m];
      Cv[m] += C[m]*C[m];
      Xa[m] += X[m];
      Xv[m] += X[m]*X[m];
      Ba[m] += BP[m];
      Bv[m] += BP[m]*BP[m];
      Mag2a[m] += Mag2[m];
      Mag2v[m] += Mag2[m]*Mag2[m];
      absMa[m] += absMag[m];
      absMv[m] += absMag[m]*absMag[m];
      absXa[m] += absX[m];
      absXv[m] += absX[m]*absX[m];
    }
  }
  /* end of per each block */

  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of environment or temperature */
    Ea[m] /= (double)N_blk;
    Ev[m] /= (double)N_blk;
    Ee[m]  = sqrt(fabs(Ev[m]-Ea[m]*Ea[m])/(double)(N_blk-1));
    Ma[m] /= (double)N_blk;
    Mv[m] /= (double)N_blk;
    Me[m]  = sqrt(fabs(Mv[m]-Ma[m]*Ma[m])/(double)(N_blk-1));
    Ca[m] /= (double)N_blk;
    Cv[m] /= (double)N_blk;
    Ce[m]  = sqrt(fabs(Cv[m]-Ca[m]*Ca[m])/(double)(N_blk-1));
    Xa[m] /= (double)N_blk;
    Xv[m] /= (double)N_blk;
    Xe[m]  = sqrt(fabs(Xv[m]-Xa[m]*Xa[m])/(double)(N_blk-1));
    Ba[m] /= (double)N_blk;
    Bv[m] /= (double)N_blk;
    Be[m]  = sqrt(fabs(Bv[m]-Ba[m]*Ba[m])/(double)(N_blk-1));
    Mag2a[m] /= (double)N_blk;
    Mag2v[m] /= (double)N_blk;
    Mag2e[m]  = sqrt(fabs(Mag2v[m]-Mag2a[m]*Mag2a[m])/(double)(N_blk-1));
    absMa[m] /= (double)N_blk;
    absMv[m] /= (double)N_blk;
    absMe[m]  = sqrt(fabs(absMv[m]-absMa[m]*absMa[m])/(double)(N_blk-1));
    absXa[m] /= (double)N_blk;
    absXv[m] /= (double)N_blk;
    absXe[m]  = sqrt(fabs(absXv[m]-absXa[m]*absXa[m])/(double)(N_blk-1));
  }

  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of temperature */
    /*                 ^^^^^^^^^^^ */
    printf("%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e \
      %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
      T[m],Ea[m]/NN,Ee[m]/NN,Ca[m]/NN,Ce[m]/NN,Ma[m]/NN,Me[m]/NN,Xa[m]/NN,Xe[m]/NN,Ba[m],Be[m],
      Mag2a[m]/NN/NN,Mag2e[m]/NN/NN,absMa[m]/NN,absMe[m]/NN,absXa[m]/NN,absXe[m]/NN);
  }

  timef = gettimeofday_sec();
  printf("# %13.6e\n", timef-timei);

  free(spins[0][0]);
  free(spins[0]);
  free(spins);
  free(BoltzFac[0]);
  free(BoltzFac);

  return 0;
}
