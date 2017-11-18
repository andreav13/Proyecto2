#include <iostream>
#include <fstream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=1;
const int Ly=64;
const int Q=9;
const double W0=4/9.;

const double Uentrada=0.0;
const double RHOinicial=1.0;

const double tau=12.0;
const double Utau=1./tau;
const double UmUtau=1-Utau;
const double UmU2tau=1.0-1.0/(2*tau);

const double g=0.01;


class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; //V[alpha][i] alpha=1 es x, alpha=0 es y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; //f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew);
  double Ux(int ix, int iy, bool UseNew);
  double Uy(int ix, int iy, bool UseNew);
  double fequilibrio(int i, double rho0, double Ux0, double Uy0);
  double Fi(int i, double rho0, double Ux0, double Uy0,double Fx0,double Fy0);
  void Inicie(double rho0, double Ux0, double Uy0);
  void ImponerCampos(int ix, int iy, double & rho0, double & Ux0, double & Uy0, int t);
  double GetFx(int ix, int iy, int t){return g;};
  double GetFy(int ix, int iy, int t){return 0;};
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
  
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
    suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Ux(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[0][i]*fnew[ix][iy][i];
    else  
      suma+=V[0][i]*f[ix][iy][i];
  return suma/rho(ix,iy,UseNew);
}

double LatticeBoltzmann::Uy(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[1][i]*fnew[ix][iy][i];
    else
      suma+=V[1][i]*f[ix][iy][i];
  return suma/rho(ix,iy,UseNew);
}

double LatticeBoltzmann::fequilibrio(int i, double rho0, double Ux0, double Uy0){
  double U2, UdotVi;
  U2=Ux0*Ux0+Uy0*Uy0;   UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  return w[i]*rho0*(1 + 3*UdotVi + 9/2.*UdotVi*UdotVi - 3/2.*U2);
  
}

double LatticeBoltzmann::Fi(int i, double rho0, double Ux0, double Uy0,double Fx0,double Fy0){
  double A,B,C;
  A=(V[0][i]-Ux0)*Fx0+(V[1][i]-Uy0)*Fy0; B=V[0][i]*Ux0+V[1][i]*Uy0;  C=V[0][i]*Fx0+V[1][i]*Fy0;
  return UmU2tau*w[i]*(3*A+9*B*C);
}

void LatticeBoltzmann::Inicie(double rho0, double Ux0, double Uy0){
int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	f[ix][iy][i]=fequilibrio(i,rho0,Ux0,Uy0);
      }
}

void LatticeBoltzmann::ImponerCampos(int ix, int iy, double & rho0, double & Ux0, double & Uy0, int t){

  //Las paredes.
  if(iy==0 || iy==Ly-1){Ux0=Uy0=0;}
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Ux0,Uy0; double Fx0,Fy0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //para cada celda
      
      rho0=rho(ix,iy,false); Ux0=Ux(ix,iy,false); Uy0=Uy(ix,iy,false); //calculo campos
      Fx0=GetFx(ix,iy,t);  Fy0=GetFy(ix,iy,t);
      Ux0+=0.5*Fx0; Uy0+=0.5*Fy0;      
      ImponerCampos(ix,iy,rho0,Ux0,Uy0,t);
      
      for(i=0;i<Q;i++) //para cada direccion
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*fequilibrio(i,rho0,Ux0,Uy0)+Fi(i,rho0,Ux0,Uy0,Fx0,Fy0); //evoluciono
    }
}
   
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}

void LatticeBoltzmann::Imprimase(char const * NombreArchivo, int t){
  double rho0,Ux0,Uy0; double Fx0,Fy0;
  ofstream MiArchivo(NombreArchivo); 
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false); Ux0=Ux(ix,iy,false); Uy0=Uy(ix,iy,false); //calculo campos
      Fx0=GetFx(ix,iy,t);  Fy0=GetFy(ix,iy,t);
      Ux0+=0.5*Fx0; Uy0+=0.5*Fy0;      
      ImponerCampos(ix,iy,rho0,Ux0,Uy0,t);
      MiArchivo<<ix<<" "<<iy<<" "<<Ux0<<" "<<Uy0<<endl;
    }
  MiArchivo.close();
}


//-----------------------Funciones Globales-----------------------


int main(void){

  LatticeBoltzmann Poiseuille;
  int t,tmax=10000;

  //  modelo tridiminseional latticce boltzmann para potenciales electroestaticos
  
  //Inicie
  Poiseuille.Inicie(RHOinicial,Uentrada,0);

  //Corra
  for(t=0;t<tmax;t++){
    Poiseuille.Colisione(t);
    Poiseuille.Adveccione();
  }
  
  Poiseuille.Imprimase("Poiseuille.dat", t);
  
  return 0;
}
