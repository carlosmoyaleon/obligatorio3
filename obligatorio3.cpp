#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <algorithm>


using namespace std;

#define h 0.1
#define nciclos N/8
#define pi 3.141592
#define npasos 1000
#define N 200

int main(void){

    
    int j, n;
    double  norma, lambda=0.3, s, k0, V[N];

    complex<double> i, phi[N], alfa[N], beta[N], chi[N];
    
    ofstream fich, schrodinger_fich;
     

    //Inicializamos los valores

    
i=complex<double> (0.0,1.0);
    

//Calculamos k0

k0=nciclos*2*pi/(N);
   

    //Calculamos s

    s=(1./(4.*k0*k0));


    //Iniciamos el potencial
    

    for ( j = 0; j < N; j++)
    {
        if((j>=(2*N/5))&&(j<=(3*N/5))){
            V[j]=lambda*k0*k0;
        }
        else{
            V[j]=0.;
        }
    } 

        //Generamos los paquetes de onda

    phi[0]=0.;
    phi[N-1]=0.;

  for ( j = 1; j < N-1; j++)
    {

        phi[j]= exp(j*1.0*k0*i)*exp(-8.*(4.*j-N)*(4.*j-N)/(N*N));
    } 






    
    

    //Calculamos alfa

    alfa[N-2]=0.;

      for ( j = N-2; j > 0; j--)
    {
        alfa[j-1]=-1./(-2.+2.*(i/s)+alfa[j]-V[j]);
    }
    


fich.open ("norma.dat");
schrodinger_fich.open("schrodinger_data.dat");


for ( n = 0; n < npasos; n++)
{
    norma=0;

    for ( j = 0; j < N; j++)
    {
           norma=norma + norm(phi[j]);
            schrodinger_fich << j*h << ", " << norm(phi[j]) << ", " << V[j] << endl;
        

    }

//Calculamos beta

    beta[N-2]=0.;

    for ( j = N-2; j > 0; j--){
          
           beta[j-1]=(1./(-2.+2.*(i/s)+alfa[j]-V[j]))*(4.*i*(phi[j]/s)-beta[j]);

    }

    
          

    //Calculamos chi

       chi[0]=0.;    
        for ( j = 0; j < N-2 ; j++)
        {
            chi[j+1]=alfa[j]*chi[j]+beta[j];
        }

    //Calculamos funcion de onda en j, n+1

    for ( j = 0; j < N; j++)
    {
        phi[j]=chi[j]-phi[j];
    }

    //Guardamos para cada salto de tiempo

    fich << norma << endl;

    schrodinger_fich << endl;

  
}

fich.close();
schrodinger_fich.close();

    return 0;
}