#include <iostream>
#include <fstream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=160; //con 600 sale segmentation fault
const int Ly=100;
const int Q=9;
const double C=0.5;
const double W0=4/9.;


const double Uentrada=0.1;
const double RHOinicial=1.0;


const double tau=0.5;
const double Utau=1./tau;
const double UmUtau=1-Utau;

const double s7=Utau;
const double s6=1.64;
const double s5=1.14;
const double s4=1.92;
const double s3=3*((2-s7)/(3-s7));


class LatticeBoltzmann{
private:
  
  double w[Q];
  int V[2][Q]; //V[alpha][i] alpha=1 es x, alpha=0 es y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; //f[ix][iy][i]
  
  double zeta[Lx][Ly][Q], deltazeta[Lx][Ly][Q];
  double S[Q];
  double MporF[Q], MporFnew[Q],  M1porDeltazeta[Q], M1porZeta[Q];
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew, double sigma);
  double Jx(int ix, int iy,bool UseNew);
  double Jy(int ix, int iy,bool UseNew);
  double Ccelda(int ix, int iy);
  double feq(int i, int ix, int iy, double rho0, double Jx0, double Jy0);
  double zetaequilibrio(int i, double rho0, double Jx0, double Jy0);
  void Inicie(double rho0, double Jx0, double Jy0);
  double GetSigma(int ix, int iy, int t);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo, int t);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1/9.;
  w[5]=w[6]=w[7]=w[8]=1/36.;

  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1;  V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;   V[1][4]=-1;
  
  V[0][5]=1;  V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1;
  V[1][5]=1;  V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;

  S[0]=S[1]=S[2]=0;
  S[3]=-s3; S[4]=-s4; S[5]=-s5; S[6]=-s6; S[7]=S[8]=-s7;
  
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew, double sigma){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma+(1./2)*sigma;
}

double LatticeBoltzmann::Jx(int ix, int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[0][i]*fnew[ix][iy][i];
    else
      suma+=V[0][i]*f[ix][iy][i];
  
  return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[1][i]*fnew[ix][iy][i];
    else
      suma+=V[1][i]*f[ix][iy][i];
  
  return suma;
}

double LatticeBoltzmann::Ccelda(int ix, int iy){
  return C;
}

double LatticeBoltzmann::feq(int i, int ix, int iy, double rho0, double Jx0, double Jy0){
  if(i==0)
    return (1-3*Ccelda(ix,iy)*Ccelda(ix,iy)*(1-W0))*rho0;
  else
    return w[i]*(3*Ccelda(ix,iy)*Ccelda(ix,iy)*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}

double LatticeBoltzmann::zetaequilibrio(int i, double rho0, double Jx0, double Jy0){
  
  double qx_eq,qy_eq;
  double energia_eq,energia2_eq,A;
  double pxx_eq,pxy_eq;

  A=3*((Jx0*Jx0+Jy0*Jy0)/rho0);
  
  qx_eq=-Jx0;                  qy_eq=-Jy0;
  
  energia_eq= -2*rho0 + A;     energia2_eq= rho0- A;

  pxx_eq= (Jx0*Jx0-Jy0*Jy0)/rho0;
  pxy_eq=Jx0*Jy0/rho0;

  double zeta_eq[Q] = {rho0,Jx0,Jy0,qx_eq,qy_eq,energia2_eq,energia_eq,pxx_eq,pxy_eq};
  return zeta_eq[i];

}

void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){
  int ix,iy,i,j;
  double a=1/9.,b=1/6.,c=1/18.,d=1/36.,e=1/4.,g=1/12.;

 double M1[Q][Q] = {{a,0,0, 0,0,a, -a,0,0},    //version presentacion
		     {a,b,0,-b,0,-c,-d,e,0},
		     {a,b,b, g,g, d, c,0,e},
		     {a,0,b, 0,-b,-c,-d,-e,0},
		     {a,-b,b,-g,g,d,c,0,-e},
		     {a,-b,0,b,0,-c,-d,e,0},
		     {a,-b,-b,-g,-g,d,c,0,e},
		     {a,0,-b,0,b,-c,-d,-e,0},
		     {a,b,-b,g,-g,d,c,0,-e}};
  

  /* double M1[Q][Q] = {{a,-a, a, 0, 0, 0, 0, 0, 0},  //version articulo
		     {a,-d,-c, b,-b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0, b,-b,-e, 0},
		     {a,-d,-c,-b, b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0,-b, b,-e, 0},
		     {a, c, d, b, g, b, g, 0, e},
		     {a, c, d,-b,-g, b, g, 0,-e},
		     {a, c, d,-b,-g,-b,-g, 0, e},
		     {a, c, d, b, g,-b,-g, 0,-e}};*/
  
 int M[Q][Q]={  { 1,  1, 1, 1,  1, 1, 1, 1, 1},   //version presentacion
		 { 0,  1, 1, 0, -1,-1,-1, 0, 1},
		 { 0,  0, 1, 1,  1, 0,-1,-1, -1},
		 { 0, -2, 1, 0, -1, 2,-1, 0, 1},
		 { 0,  0, 1,-2,  1, 0,-1, 2, -1},
		 { 4, -2, 1,-2,  1,-2, 1,-2, 1},
		 {-4, -1, 2,-1,  2,-1, 2,-1, 2},
		 { 0,  1, 0,-1,  0, 1, 0,-1, 0},
		 { 0,  0, 1, 0, -1, 0, 1, 0, -1}};
  
 /*  int M[Q][Q]={{1,1,1,1,1,1,1,1,1},                  //version articulo
	       {-4,-1,-1,-1,-1,2,2,2,2},
	       {4,-2,-2,-2,-2,1,1,1,1},
	       {0,1,0,-1,0,1,-1,-1,1},
	       {0,-2,0,2,0,1,-1,-1,1},
	       {0,0,1,0,-1,1,1,-1,-1},
	       {0,0,-2,0,2,1,1,-1,-1},
	       {0,1,-1,1,-1,0,0,0,0},
	       {0,0,0,0,0,1,-1,1,-1}};*/
  
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
	//	zeta[ix][iy][i]=zetaequilibrio(i,rho0,Jx0,Jy0);
	f[ix][iy][i]=feq(i,ix,iy,rho0,Jx0,Jy0);
      }
      
      for(i=0;i<Q;i++){
	//M1porZeta[i]=0;
	MporF[i]=0;
	for(j=0;j<Q;j++){
	  //   M1porZeta[i]+=M1[i][j]*zeta[ix][iy][j];
	  MporF[i]+=M[i][j]*f[ix][iy][j];
	}
      }
      
      for(i=0;i<Q;i++){
	//f[ix][iy][i]=M1porZeta[i];
	zeta[ix][iy][i]=MporF[i];
      }
      
    }
  }
}



double LatticeBoltzmann::GetSigma(int ix, int iy, int t){
  double A=1, lambda=1000, omega=2*M_PI*Ccelda(ix,iy)/lambda;
  if(ix==0)
    return A*sin(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i,j; double rho0,Jx0,Jy0,sigma;
  double a=1/9.,b=1/6.,c=1/18.,d=1/36.,e=1/4.,g=1/12.;
  
  double M1[Q][Q] = {{a,0,0, 0,0,a, -a,0,0},    //version presentacion
		     {a,b,0,-b,0,-c,-d,e,0},
		     {a,b,b, g,g, d, c,0,e},
		     {a,0,b, 0,-b,-c,-d,-e,0},
		     {a,-b,b,-g,g,d,c,0,-e},
		     {a,-b,0,b,0,-c,-d,e,0},
		     {a,-b,-b,-g,-g,d,c,0,e},
		     {a,0,-b,0,b,-c,-d,-e,0},
		     {a,b,-b,g,-g,d,c,0,-e}};
 
 
  /*  double M1[Q][Q] = {{a,-a, a, 0, 0, 0, 0, 0, 0},  //version articulo
		     {a,-d,-c, b,-b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0, b,-b,-e, 0},
		     {a,-d,-c,-b, b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0,-b, b,-e, 0},
		     {a, c, d, b, g, b, g, 0, e},
		     {a, c, d,-b,-g, b, g, 0,-e},
		     {a, c, d,-b,-g,-b,-g, 0, e},
		     {a, c, d, b, g,-b,-g, 0,-e}};*/
  

  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){ //para cada celda      
      
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,false,sigma); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false); //calculo campos
 
 
      for(i=0;i<Q;i++){
	deltazeta[ix][iy][i]= S[i]*(zeta[ix][iy][i]-zetaequilibrio(i,rho0,Jx0,Jy0));
      }
      
      for(i=0;i<Q;i++){
	M1porDeltazeta[i]=0;
	for(j=0;j<Q;j++)
	  M1porDeltazeta[i]+=M1[i][j]*deltazeta[ix][iy][j];
      }
      
      for(i=0;i<Q;i++){
	//	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*fequilibrio(i,rho0,Jx0,Jy0);
	fnew[ix][iy][i]=f[ix][iy][i]+M1porDeltazeta[i];
      }                                                //evoluciono
      
    }
  }

}
   
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i,j;

  int M[Q][Q]={  { 1,  1, 1, 1,  1, 1, 1, 1, 1},  //version presentacion
		 { 0,  1, 1, 0, -1,-1,-1, 0, 1},
		 { 0,  0, 1, 1,  1, 0,-1,-1, -1},
		 { 0, -2, 1, 0, -1, 2,-1, 0, 1},
		 { 0,  0, 1,-2,  1, 0,-1, 2, -1},
		 { 4, -2, 1,-2,  1,-2, 1,-2, 1},
		 {-4, -1, 2,-1,  2,-1, 2,-1, 2},
		 { 0,  1, 0,-1,  0, 1, 0,-1, 0},
		 { 0,  0, 1, 0, -1, 0, 1, 0, -1}};

  /*  int M[Q][Q]={{1,1,1,1,1,1,1,1,1},                  //version articulo
	       {-4,-1,-1,-1,-1,2,2,2,2},
	       {4,-2,-2,-2,-2,1,1,1,1},
	       {0,1,0,-1,0,1,-1,-1,1},
	       {0,-2,0,2,0,1,-1,-1,1},
	       {0,0,1,0,-1,1,1,-1,-1},
	       {0,0,-2,0,2,1,1,-1,-1},
	       {0,1,-1,1,-1,0,0,0,0},
	       {0,0,0,0,0,1,-1,1,-1}};*/
  
 
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      
      for(i=0;i<Q;i++){
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
      }
      
      
      for(i=0;i<Q;i++){
	MporFnew[i]=0;
	for(j=0;j<Q;j++){
	  MporFnew[i]+=M[i][j]*fnew[ix][iy][j];
	}
      }
      
      for(i=0;i<Q;i++)
	zeta[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=MporFnew[i];
      
    }
  }
}

void LatticeBoltzmann::Imprimase(char const * NombreArchivo, int t){
  double rho0,Jx0,Jy0,sigma;
  ofstream MiArchivo(NombreArchivo); 
  for(int ix=0;ix<Lx;ix+=4)
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true,sigma); Jx0=Jx(ix,iy,true); Jy0=Jy(ix,iy,true);
      sigma=GetSigma(ix,iy,t);
      MiArchivo<<ix<<" "<<iy<<" "<<4.0/Uentrada*Jx0/rho0<<" "<<4.0/Uentrada*Jy0/rho0<<endl;
      //MiArchivo<<ix<<" "<<iy<<" "<<Jx0<<" "<<Jy0<<endl;
    }
  MiArchivo.close();
}



//-----------------------Funciones Globales-----------------------


int main(void){
  
  LatticeBoltzmann Ala;
  int t,tmax=1000;
  
  //Inicie
  Ala.Inicie(RHOinicial,Uentrada,0);
 
  //Corra
  for(t=0;t<tmax;t++){
    /*double Rho0,Jx0,Jy0,Ux0,Uy0; int ix,iy;
    for (ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
	Rho0=Ala.rho(ix,iy,false);
	Jx0=Ala.Jx(ix,iy,false);
	Jy0=Ala.Jy(ix,iy,false);
	Ala.ImponerCampos(ix,iy,Rho0,Jx0,Jy0,t);
	Ux0=Jx0/Rho0;
	Uy0=Jy0/Rho0;
	
      }if(Rho0==1)cout<<ix<<"\t"<<iy<<"\t"<<Rho0<<"\t"<<Ux0<<"\t"<<Uy0<<endl;
    }
    cout<<"--------------------Nuevo tiempo = "<<t<<" ---------------------------------------"<<endl;      */
    
    Ala.Colisione(t);
    Ala.Adveccione();
  }
  
  Ala.Imprimase("Ala_sigma.dat", t);

 
  
  return 0;
}
