/*-------------------------------------------------------------------*
 * Ising (square lattice)
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
 * show spin configuration
 *-------------------------------------------------------------------*/
int show_spins(int n, int **spins)
{
  int i,j;
  
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      printf("%2d ", spins[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  
  return 0;
}

/*-------------------------------------------------------------------*
 * energy = -J sum_i,j spins[i][j] * (spins[i+1][j] + spins[i][j+1])
 *-------------------------------------------------------------------*/
double calc_energy(int n, double J, int **spins)
{
  int i,ipp,j,jpp;
  int ene;
  
  ene = 0.0;
  for(i=0; i<n; i++){
    ipp = (i==n-1) ? 0 : i+1;
    for(j=0; j<n; j++){
      jpp = (j==n-1) ? 0 : j+1;
      ene -= spins[i][j] * (spins[ipp][j] + spins[i][jpp]);
    }
  }
  return(J*ene);
}

/*-------------------------------------------------------------------*
 * magnetization = sum_i,j spins[i][j]
 *-------------------------------------------------------------------*/
double calc_mag(int n, int **spins)
{
  int i,j;
  double mag;
  
  mag = 0.0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      mag += spins[i][j];
    }
  }
  return(mag);
}

/*-------------------------------------------------------------------*
 * grow cluster
 *-------------------------------------------------------------------*/
void grow_cluster(int n, double *BoltzFac, int **spins, int i, int j, int mag_ij)
{
  int ipp,imm,jpp,jmm;
  int new_i,new_j;
  
  spins[i][j] *= -1;
  
  ipp = (i==n-1) ? 0 : i+1;
  imm = (i==0) ? n-1 : i-1;
  jpp = (j==n-1) ? 0 : j+1;
  jmm = (j==0) ? n-1 : j-1;
  
  new_i = ipp;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = imm;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = i;
  new_j = jpp;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = i;
  new_j = jmm;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
}

/*-------------------------------------------------------------------*
 * one monte carlo step
 *-------------------------------------------------------------------*/
int one_flip(int n, double *BoltzFac, int **spins)
{
  int i,j;
  int mag_ij;
  
  /* spin flip candidate */
  i = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  mag_ij = spins[i][j];
  
  /* grow and flip a cluster */
  grow_cluster(n,BoltzFac,spins,i,j,mag_ij);
  
  return 0;
}

/*-------------------------------------------------------------------*
 * grow cluster
 * (calculate energy and magnetization differences)
 *-------------------------------------------------------------------*/
void grow_cluster_ene_mag(int n, double J, double *BoltzFac, int **spins, int i, int j, int mag_ij,
  double *ene, double *mag)
{
  int egap;
  int ipp,imm,jpp,jmm;
  int new_i,new_j;
  
  ipp = (i==n-1) ? 0 : i+1;
  imm = (i==0) ? n-1 : i-1;
  jpp = (j==n-1) ? 0 : j+1;
  jmm = (j==0) ? n-1 : j-1;
  
  spins[i][j] *= -1;
  
  /* calculate energy and magnetization */
  egap = 2 * spins[i][j]
    *(spins[ipp][j] + spins[i][jpp]
    + spins[imm][j] + spins[i][jmm]);
  ene[0] += -J * egap;
  mag[0] += 2 * spins[i][j];
  /* end of calculate energy and magnetization */
  
/*
  new_i = ipp;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = imm;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = i;
  new_j = jpp;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
  new_i = i;
  new_j = jmm;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster(n,BoltzFac,spins,new_i,new_j,mag_ij);
  }
*/
  
  new_i = ipp;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster_ene_mag(n,J,BoltzFac,spins,new_i,new_j,mag_ij,ene,mag);
  }
  new_i = imm;
  new_j = j;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster_ene_mag(n,J,BoltzFac,spins,new_i,new_j,mag_ij,ene,mag);
  }
  new_i = i;
  new_j = jpp;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster_ene_mag(n,J,BoltzFac,spins,new_i,new_j,mag_ij,ene,mag);
  }
  new_i = i;
  new_j = jmm;
  if(spins[new_i][new_j]==mag_ij && dsfmt_genrand_close_open(&dsfmt)<BoltzFac[0]){
    grow_cluster_ene_mag(n,J,BoltzFac,spins,new_i,new_j,mag_ij,ene,mag);
  }

}

/*-------------------------------------------------------------------*
 * one monte carlo step
 * (calculate energy and magnetization differences)
 *-------------------------------------------------------------------*/
int one_flip_ene_mag(int n, double J, double *BoltzFac, int **spins,
  double *ene, double *mag)
{
  int i,j;
  int mag_ij;
  
  /* spin flip candidate */
  i = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (n * dsfmt_genrand_close_open(&dsfmt));
  mag_ij = spins[i][j];
  
  /* grow and flip a cluster */
  grow_cluster_ene_mag(n,J,BoltzFac,spins,i,j,mag_ij,ene,mag);

  /* calculate physical quantity */
/*
  int ipp,imm,jpp,jmm;
  ipp = (i==n-1) ? 0 : i+1;
  imm = (i==0) ? n-1 : i-1;
  jpp = (j==n-1) ? 0 : j+1;
  jmm = (j==0) ? n-1 : j-1;
  double tmp_ene;
  tmp_ene = 0.0;
  for(i=0; i<n; i++){
    ipp = (i==n-1) ? 0 : i+1;
    for(j=0; j<n; j++){
      jpp = (j==n-1) ? 0 : j+1;
      tmp_ene -= spins[i][j] * (spins[ipp][j] + spins[i][jpp]);
    }
  }
  ene[0] = tmp_ene;
  double tmp_mag;
  tmp_mag = 0.0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      tmp_mag += spins[i][j];
    }
  }
  mag[0] = tmp_mag;
*/
  /* end of calculate physical quantity */
  
  return 0;
}

/* specific heat */
/*---------------------------------------------------------------*
 * C = (<E^2> - <E>^2) / (kT^2)
 * <E> = sum_i E_i / N_smp
 * <E^2> = sum_i E_i^2 / N_smp
 *---------------------------------------------------------------*/

/* magnetic susceptibility */
/*---------------------------------------------------------------*
 * X = (<M^2> - <M>^2) / (kT)
 * <M> = sum_i M_i / N_smp
 * <M^2> = sum_i M_i^2 / N_smp
 *---------------------------------------------------------------*/

/* Binder parameter */
/*---------------------------------------------------------------*
 * U_4 = 1 - <M^4>/(3 * <M^2>^2)
 *---------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  int i,j;
  char opt;
  int N = 8;
  int seed = N * 10000;
  int NN = N*N;
  int N_eq_per_spin  = 10000;
  int N_smp_per_spin = 10000;
  int N_cor_per_spin = 10000;
  int N_blk = 5;
  double T_i = 2.10;
  double T_f = 2.40;
  double T_width = 0.005;
  double J = +1.0;/* ferromagnetic */
//  double J = -1.0;/* antiferromagnetic */
  
  /* prameter settings */
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 'n'){
        N = atoi(argv[i+1]);
        seed = N * 10000;
        NN = N*N;
      }else if(opt == 'e'){
        N_eq_per_spin  = atoi(argv[i+1]);
      }else if(opt == 's'){
        N_smp_per_spin = atoi(argv[i+1]);
      }else if(opt == 'c'){
        N_cor_per_spin = atoi(argv[i+1]);
      }else if(opt == 'b'){
        N_blk = atoi(argv[i+1]);
      }else if(opt == 'i'){
        T_i = atof(argv[i+1]);
      }else if(opt == 'f'){
        T_f = atof(argv[i+1]);
      }else if(opt == 'w'){
        T_width = atof(argv[i+1]);
      }else if(opt == 'j'){
        J = atof(argv[i+1]);
      }else if(opt == 'r'){
        seed = atoi(argv[i+1]);
      }
    }
  }
  
  double N_eq  = N_eq_per_spin  * NN;
  double N_smp = N_smp_per_spin * NN;
  double N_cor = N_cor_per_spin * NN;
  double BoltzFac[1];
  double T;
  double Etemp,Mtemp;
  double Ene,Ene2;
  double Mag,Mag2,Mag4;
  double C,X,BP;
  double absMag,absX;
  double Ea,Ma,Ca,Xa,Ba;
  double Ev,Mv,Cv,Xv,Bv;
  double Ee,Me,Ce,Xe,Be;
  double absMa,absXa,Mag2a;
  double absMv,absXv,Mag2v;
  double absMe,absXe,Mag2e;
  double Eadd[1],Madd[1];
  double timei,timef;
  int **spins;
  
  spins = (int**)malloc(N*sizeof(int*));
  for(i=0; i<N; i++) spins[i] = (int*)malloc(N*sizeof(int));
  
  printf("# N=%d NN=%d ",
    N,NN);
  printf("N_eq_per_spin=%d N_smp_per_spin=%d N_cor_per_spin=%d N_blk=%d\n",
    N_eq_per_spin,N_smp_per_spin,N_cor_per_spin,N_blk);
  printf("# seed=%d T_i=%9.2e T_f=%9.2e T_width=%9.2e J=%9.2e\n",
    seed,T_i,T_f,T_width,J);
  printf("# T E dE C dC M dM X dX BP dBP rnd time\n");
  
  /* random number seed initialization */
  dsfmt_init_gen_rand(&dsfmt,seed);
  /* end of random number seed initialization */
  
  /* hot or cold start */
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      spins[i][j] = (dsfmt_genrand_close_open(&dsfmt) > 0.5) ? +1 : -1;
/*
      spins[i][j] = 1;
*/
    }
  }
  /* end of hot or cold start */
  
  /* calculate observables from temperature T_i to T_f */
  for(T=T_f; T>T_i-T_width/2.0; T-=T_width){
    
    timei = gettimeofday_sec();
    
    /* random number seed initialization */
    seed++;
    dsfmt_init_gen_rand(&dsfmt,seed);
    /* end of random number seed initialization */
    
    /* table of Boltzmann factors */
//    BoltzFac[0] = exp(-2.0*J/T);
    BoltzFac[0] = 1.0 - exp(-2.0*J/T);
    /* end of table of Boltzmann factors */
    
    /* loop until equilibrium is achieved */
    for(i=0; i<N_eq; i++){
      one_flip(N,BoltzFac,spins);
    }
    /* end of loop until equilibrium is achieved */
    
    /* parameter initialization */
    Ea = 0.0;  Ma = 0.0;  Ca = 0.0;  Xa = 0.0;  Ba = 0.0;
    Ev = 0.0;  Mv = 0.0;  Cv = 0.0;  Xv = 0.0;  Bv = 0.0;
    Ee = 0.0;  Me = 0.0;  Ce = 0.0;  Xe = 0.0;  Be = 0.0;
    Mag2a = 0.0;  Mag2v = 0.0;  Mag2e = 0.0;
    absMa = 0.0;  absXa = 0.0;
    absMv = 0.0;  absXv = 0.0;
    absMe = 0.0;  absXe = 0.0;
    /* end of parameter initialization */
    
    /* per each block */
    for(i=0; i<N_blk; i++){
      
      /* loop until correlation is neglected */
      for(j=0; j<N_cor; j++){
        one_flip(N,BoltzFac,spins);
      }
      /* end of loop until correlation is neglected */
      
      /* parameter initialization per each block*/
      Eadd[0] = calc_energy(N,J,spins);
      Madd[0] = calc_mag(N,spins);
      Etemp = 0.0;  Mtemp = 0.0;
      Ene = 0.0;  Ene2 = 0.0;
      Mag = 0.0;  Mag2 = 0.0;  Mag4 = 0.0;
      absMag = 0.0;
      C = 0.0;  X = 0.0;  BP = 0.0;
      /* end of parameter initialization per each block*/
      
      /* measurement loop */
      for(j=0; j<N_smp; j++){
        one_flip_ene_mag(N,J,BoltzFac,spins,Eadd,Madd);
        Etemp = Eadd[0];/* E */
        Mtemp = Madd[0];/* M */
        Ene  += Etemp;
        Etemp = Etemp*Etemp;/* E^2 */
        Ene2 += Etemp;
        Mag  += Mtemp;
        Mtemp = Mtemp*Mtemp;/* M^2 */
        Mag2 += Mtemp;
        Mtemp = Mtemp*Mtemp;/* M^4 */
        Mag4 += Mtemp;
        Mtemp   = fabs(Madd[0]);/* |M| */
        absMag += Mtemp;
      }
      /* end of measurement loop */
      
      Ene  /= (double)N_smp;
      Ene2 /= (double)N_smp;
      Mag  /= (double)N_smp;
      Mag2 /= (double)N_smp;
      Mag4 /= (double)N_smp;
      C = (Ene2 - Ene*Ene)/T/T;
      X = (Mag2 - Mag*Mag)/T;
      BP = 1.0 - Mag4/Mag2/Mag2/3.0;
      absMag /= (double)N_smp;
      absX = (Mag2 - absMag*absMag)/T;
      
/*
      printf("%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
        T,Ene/NN,C/NN,Mag/NN,X/NN,BP);
*/
      
      Ea += Ene;
      Ev += Ene*Ene;
      Ma += Mag;
      Mv += Mag*Mag;
      Ca += C;
      Cv += C*C;
      Xa += X;
      Xv += X*X;
      Ba += BP;
      Bv += BP*BP;
      Mag2a += Mag2;
      Mag2v += Mag2*Mag2;
      absMa += absMag;
      absMv += absMag*absMag;
      absXa += absX;
      absXv += absX*absX;
    }
    /* end of per each block */
    
    Ea /= (double)N_blk;
    Ev /= (double)N_blk;
    Ee  = sqrt(fabs(Ev-Ea*Ea)/(double)(N_blk-1));
    Ma /= (double)N_blk;
    Mv /= (double)N_blk;
    Me  = sqrt(fabs(Mv-Ma*Ma)/(double)(N_blk-1));
    Ca /= (double)N_blk;
    Cv /= (double)N_blk;
    Ce  = sqrt(fabs(Cv-Ca*Ca)/(double)(N_blk-1));
    Xa /= (double)N_blk;
    Xv /= (double)N_blk;
    Xe  = sqrt(fabs(Xv-Xa*Xa)/(double)(N_blk-1));
    Ba /= (double)N_blk;
    Bv /= (double)N_blk;
    Be  = sqrt(fabs(Bv-Ba*Ba)/(double)(N_blk-1));
    Mag2a /= (double)N_blk;
    Mag2v /= (double)N_blk;
    Mag2e  = sqrt(fabs(Mag2v-Mag2a*Mag2a)/(double)(N_blk-1));
    absMa /= (double)N_blk;
    absMv /= (double)N_blk;
    absMe  = sqrt(fabs(absMv-absMa*absMa)/(double)(N_blk-1));
    absXa /= (double)N_blk;
    absXv /= (double)N_blk;
    absXe  = sqrt(fabs(absXv-absXa*absXa)/(double)(N_blk-1));
    
    timef = gettimeofday_sec();
    
    printf("%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e \
      %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %d %13.6e\n",
      T,Ea/NN,Ee/NN,Ca/NN,Ce/NN,Ma/NN,Me/NN,Xa/NN,Xe/NN,Ba,Be,
      Mag2a/NN/NN,Mag2e/NN/NN,absMa/NN,absMe/NN,absXa/NN,absXe/NN,seed,timef-timei);
  }
  /* end of calculate observables from temperature T_i to T_f */
  
  for(i=0; i<N; i++) free(spins[i]);
  free(spins);
  
  return 0;
}

