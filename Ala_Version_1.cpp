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
  double zeta[Lx][Ly][Q], zetanew[Lx][Ly][Q];
  
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

 int  M[9][9]={  { 1,  1, 1, 1,  1, 1, 1, 1, 1},
		 { 0,  1, 1, 0, -1,-1,-1, 0, 1},
		 { 0,  0, 1, 1,  1, 0,-1,-1, -1},
		 { 0, -2, 1, 0,  1, 2,-1, 0, 1},
		 { 0,  0, 1,-2,  1, 0,-1, 2, -1},
		 { 4, -2, 1,-2,  1,-2, 1,-2, 1},
		 {-4, -1, 2,-1,  2,-1, 2,-1, 2},
		 { 0,  1, 0,-1,  0, 1, 0,-1, 0},
		 { 0,  0, 1, 0, -1, 0, 1, 0, -1}};

 
double  S[9][9]={{ 0,  0, 0, 0,  0, 0, 0, 0, 0},
		 { 0,  0, 0, 0,  0, 0, 0, 0, 0},
     		 { 0,  0, 0, 0,  0, 0, 0, 0, 0},
		 { 0,  0, 0,s3,  0, 0, 0, 0, 0},
		 { 0,  0, 0, 0, s4, 0, 0, 0, 0},
		 { 0,  0, 0, 0,  0,s5, 0, 0, 0},
		 { 0,  0, 0, 0,  0, 0,s6, 0, 0},
		 { 0,  0, 0, 0,  0, 0, 0,s7, 0},
		 { 0,  0, 0, 0,  0, 0, 0, 0,s7}};
/*
 for (int i=0;i<9;i++){
   for(int j=0;j<9;j++){
     std::cout << S[i][j] << "  ";       //Visualizar las matrices
   }
   std::cout<<endl;
 }
*/

 
}
/*
double** S(void){                                         // Intento de una Función para la matriz S, no se como hacer que devuelva la matriz
  double S[9][9]={{ 0,  0, 0, 0,  0, 0, 0, 0, 0},
		 { 0,  0, 0, 0,  0, 0, 0, 0, 0},
		 { 0,  0, 0, 0,  0, 0, 0, 0, 0},
		 { 0,  0, 0,s3,  0, 0, 0, 0, 0},
		 { 0,  0, 0, 0, s4, 0, 0, 0, 0},
		 { 0,  0, 0, 0,  0,s5, 0, 0, 0},
		 { 0,  0, 0, 0,  0, 0,s6, 0, 0},
		 { 0,  0, 0, 0,  0, 0, 0,s7, 0},
		 { 0,  0, 0, 0,  0, 0, 0, 0,s7}};

   for(int i = 0; i < 9; ++i){
     for (int j=0; j<9;j++){
       std::cout << S[i][j] << ' ';
     }
     std::cout << std::endl;
   }
   
   return S;

}
*/
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

double LatticeBoltzmann::zetaequilibrio(int i, double rho0, double Jx0, double Jy0){
  
  double qx_eq,qy_eq;
  double energia_eq,energia2_eq,A;
  double pxx_eq,pxy_eq;

  A=3*((Jx0*Jx0*+Jy0*Jy0)/rho0);
  
  qx_eq=-Jx0;                  qy_eq=-Jy0;
  
  energia_eq= -2*rho0 + A;     energia2_eq= rho0- A;

  pxx_eq= (Jx0*Jx0-Jy0*Jy0)/rho0;
  pxy_eq=Jx0*Jy0/rho0;

  double zeta_eq[9][1]= {{rho0},{Jx0},{Jy0},{qx_eq},{qy_eq},{energia2_eq},{energia_eq},{pxx_eq},{pxy_eq}};
  return zeta_eq[i][1];


  
}

double LatticeBoltzmann::fequilibrio(int i, double rho0, double Jx0, double Jy0){
  double U2, UdotVi;
  U2=Jx0*Jx0+Jy0*Jy0;   UdotVi=Jx0*V[0][i]+Jy0*V[1][i];
  return w[i]*rho0*(1 + 3*UdotVi + 9/2.*UdotVi*UdotVi - 3/2.*U2);
  
}

void LatticeBoltzmann::Inicie(double rho0, double Jx0, double Jy0){
int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	f[ix][iy][i]=fequilibrio(i,rho0,Jx0,Jy0);
	zeta[ix][iy][i]=zetaequilibrio(i,rho0,Jx0,Jy0);           //Posiblemente hay que inicializar zeta pero no estoy seguro
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
    Jx0=Uentrada;
    Jy0=0;
  }
}

void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //para cada celda
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false); //calculo campos
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      for(i=0;i<Q;i++){ //para cada direccion
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*fequilibrio(i,rho0,Jx0,Jy0); //evoluciono
       	zetanew[ix][iy][i]= S*(zeta[ix][iy][i]+zetaequilibrio(i,rho0,Jx0,Jy0));    //No se como introducir S en esta función.
      }
    }
}
   
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
  //	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=(M^⁻1)*zetanew[ix][iy][i];   //No se como introducir M en esta funcion.
}

void LatticeBoltzmann::Imprimase(char const * NombreArchivo, int t){
  double rho0,Jx0,Jy0;
  ofstream MiArchivo(NombreArchivo); 
  for(int ix=0;ix<Lx;ix+=4)
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); Jx0=Jx(ix,iy,true); Jy0=Jy(ix,iy,true);
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      MiArchivo<<ix<<" "<<iy<<" "<<4.0/Uentrada*Jx0<<" "<<4.0/Uentrada*Jy0<<endl;
    }
  MiArchivo.close();
}


//-----------------------Funciones Globales-----------------------


int main(void){

  LatticeBoltzmann Ang;
  int t,tmax=1000;
  
  //Inicie
  Ang.Inicie(RHOinicial,Uentrada,0);

  //Corra
  for(t=0;t<tmax;t++){
    Ang.Colisione(t);
    Ang.Adveccione();
  }
  
  Ang.Imprimase("Ang.dat", t);
  
  return 0;
}
