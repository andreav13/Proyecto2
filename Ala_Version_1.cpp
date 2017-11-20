#include <iostream>
#include <fstream>
#include <cmath> 
#include "Random64.h"

using namespace std;

const int Lx=260;
const int Ly=64;
const int Q=9;
const double W0=4/9.;

const double Uentrada=0.1;
const double RHOinicial=1.0;

const double tau=0.55;
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
  double MporF[Q], M1porDeltazeta[Q];
  
public:
  LatticeBoltzmann(void);
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double fequilibrio(int i, double rho0, double Jx0, double Jy0);
  double zetaequilibrio(int i, double rho0, double Jx0, double Jy0);
  void Inicie(double rho0, double Jx0, double Jy0);
  void ImponerCampos(int ix, int iy, double & rho0, double & Jx0, double & Jy0, int t);
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

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
    suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[0][i]*fnew[ix][iy][i];
    else  
      suma+=V[0][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=V[1][i]*fnew[ix][iy][i];
    else
      suma+=V[1][i]*f[ix][iy][i];
  return suma;
}

double LatticeBoltzmann::fequilibrio(int i, double rho0, double Jx0, double Jy0){
  double U2, UdotVi;
  U2=Jx0*Jx0+Jy0*Jy0;   UdotVi=Jx0*V[0][i]+Jy0*V[1][i];
  return w[i]*rho0*(1 + 3*UdotVi + 9/2.*UdotVi*UdotVi - 3/2.*U2);
  
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

  /*int M[Q][Q]={  { 1,  1, 1, 1,  1, 1, 1, 1, 1},   //version presentacion
		 { 0,  1, 1, 0, -1,-1,-1, 0, 1},
		 { 0,  0, 1, 1,  1, 0,-1,-1, -1},
		 { 0, -2, 1, 0,  1, 2,-1, 0, 1},
		 { 0,  0, 1,-2,  1, 0,-1, 2, -1},
		 { 4, -2, 1,-2,  1,-2, 1,-2, 1},
		 {-4, -1, 2,-1,  2,-1, 2,-1, 2},
		 { 0,  1, 0,-1,  0, 1, 0,-1, 0},
		 { 0,  0, 1, 0, -1, 0, 1, 0, -1}};*/

  int M[Q][Q]={{1,1,1,1,1,1,1,1,1},                  //version articulo
	       {-4,-1,-1,-1,-1,2,2,2,2},
	       {4,-2,-2,-2,-2,1,1,1,1},
	       {0,1,0,-1,0,1,-1,-1,1},
	       {0,-2,0,2,0,1,-1,-1,1},
	       {0,0,1,0,-1,1,1,-1,-1},
	       {0,0,-2,0,2,1,1,-1,-1},
	       {0,1,-1,1,-1,0,0,0,0},
	       {0,0,0,0,0,1,-1,1,-1}};
  
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){

      for(i=0;i<Q;i++)
	f[ix][iy][i]=fequilibrio(i,rho0,Jx0,Jy0);
      //zeta[ix][iy][i]=zetaequilibrio(i,rho0,Jx0,Jy0);
      
      for(i=0;i<Q;i++){
	MporF[i]=0;
	for(j=0;j<Q;j++)
	  MporF[i]+=M[i][j]*f[ix][iy][j];
      }
      
      for(i=0;i<Q;i++)
	zeta[ix][iy][i]=MporF[i];
      
    }
}

void LatticeBoltzmann::ImponerCampos(int ix, int iy, double & rho0, double & Jx0, double & Jy0, int t){
  double ixc=Lx/8, iyc=Ly/2, R=Ly/5, R2=R*R;
  //El obstaculo
  if((ix-ixc)*(ix-ixc) + (iy-iyc)*(iy-iyc) <= R2)
    Jx0=Jy0=0;
  //Un puntito extra
  if(ix==ixc && iy==iyc+R+1)
    Jx0=Jy0=0;
  //El ventilador
  if(ix==0){
    Jx0=Uentrada*RHOinicial;
    Jy0=0;
  }
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i,j; double rho0,Jx0,Jy0;
  double a=1/9.,b=1/6.,c=1/18.,d=1/36.,e=1/4.,g=1/12.;

  /*double M1[Q][Q] = {{a,0,0, 0,0,a, -a,0,0},    //version presentacion
		     {a,b,0,-b,0,-c,-d,e,0},
		     {a,b,b, g,g, d, c,0,e},
		     {a,0,b, 0,-b,c,d,-e,0},
		     {a,-b,b,-g,g,d,c,0,-e},
		     {a,-b,0,b,0,-c,-d,e,0},
		     {a,-b,-b,-g,-g,d,c,0,e},
		     {a,0,-b,0,b,-c,-d,-e,0},
		     {a,b,-b,g,-g,d,c,0,-e}};	*/

  double M1[Q][Q] = {{a,-a, a, 0, 0, 0, 0, 0, 0},  //version articulo
		     {a,-d,-c, b,-b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0, b,-b,-e, 0},
		     {a,-d,-c,-b, b, 0, 0, e, 0},
		     {a,-d,-c, 0, 0,-b, b,-e, 0},
		     {a, c, d, b, g, b, g, 0, e},
		     {a, c, d,-b,-g, b, g, 0,-e},
		     {a, c, d,-b,-g,-b,-g, 0, e},
		     {a, c, d, b, g,-b,-g, 0,-e}};

  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //para cada celda      
      
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false); //calculo campos
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      for(i=0;i<Q;i++){
	deltazeta[ix][iy][i]= S[i]*(zeta[ix][iy][i]-zetaequilibrio(i,rho0,Jx0,Jy0));  
      }

      for(i=0;i<Q;i++){
	M1porDeltazeta[i]=0;
	for(j=0;j<Q;j++)
	  M1porDeltazeta[i]+=M1[i][j]*deltazeta[ix][iy][j];
      }
      
      for(i=0;i<Q;i++)
	fnew[ix][iy][i]=f[ix][iy][i]+M1porDeltazeta[i]; //evoluciono
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
  double rho0,Jx0,Jy0;
  ofstream MiArchivo(NombreArchivo); 
  for(int ix=0;ix<Lx;ix+=4)
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); Jx0=Jx(ix,iy,true); Jy0=Jy(ix,iy,true);
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      // MiArchivo<<ix<<" "<<iy<<" "<<4.0/Uentrada*Jx0/rho0<<" "<<4.0/Uentrada*Jy0/rho0<<endl;
      //MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      MiArchivo<<ix<<" "<<iy<<" "<<Jx0<<" "<<Jy0<<endl;
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
    Ala.Colisione(t);
    Ala.Adveccione();
  }
  
  Ala.Imprimase("Ala.dat", t);

  
  return 0;
}
