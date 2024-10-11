#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#define nx 512
#define ny 512

using namespace std;

double altura,largura;
double f2[nx][ny] = {0}, f1[nx][ny] = {0}, f4[nx][ny] = {0}, f3[nx][ny] = {0};
double f2aux[nx][ny]= {0}, f1aux[nx][ny]= {0}, f4aux[nx][ny]= {0}, f3aux[nx][ny]= {0};
double f0[nx][ny]= {0};
double f0aux[nx][ny]= {0}; //f1eções diagonais(matriz diagonal)
double f0eq[nx][ny]= {0}, f2eq[nx][ny]= {0}, f4eq[nx][ny]= {0},f1eq[nx][ny]= {0},f3eq[nx][ny]= {0};


double T[nx][ny]= {0};
double Ux[nx][ny]= {0}, Uy[nx][ny]= {0};

double D=2, dt=0.1, dx=0.01, alpha = 0.003;
double inv_omega   = (dt*D*alpha/(dx*dx) + 0.5);
double tau = dt*inv_omega;
double omega = 1./inv_omega;



void colisao()
{
    int i, j;
    double t;// ux2uy2;
    //Pesos
    double w0 = (2.0/6.0), w = (1.0/6.0);



    // COLISÃO
    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            t = f0[i][j] + f2[i][j] + f4[i][j] + f1[i][j] + f3[i][j];


            T[i][j] = t;

            f0eq[i][j]= w0*(t);
            f1eq[i][j]= w*(t);
            f2eq[i][j]= w*(t);
            f3eq[i][j]= w*(t);
            f4eq[i][j]= w*(t);


            f0aux[i][j]=(f0[i][j]*(1.-omega))+ omega*f0eq[i][j] ;
            f1aux[i][j]=(f1[i][j]*(1.-omega))+ omega*f1eq[i][j];
            f2aux[i][j]=(f2[i][j]*(1.-omega))+ omega*f2eq[i][j];
            f3aux[i][j]=(f3[i][j]*(1.-omega))+ omega*f3eq[i][j];
            f4aux[i][j]=(f4[i][j]*(1.-omega))+ omega*f4eq[i][j];

        }
    }
}

void paredes()
{
     for(int i=0; i<nx; i++)
    {
        // cima
        f4[i][ny-1] = f2[i][ny-1];

        // baixo
        f2[i][0] = f4[i][0];
    }

    for(int j=0; j<ny; j++)
    {
        //Esquerda
        f1[0][j] = f3[0][j];


        //Direita
        f3[nx-1][j] = f1[nx-1][j];

    }
}
void propagacao()
{
    int i=0,j=0, im, ip, jm, jp;
    for(i=0; i<nx; i++)
    {
        im=i-1;
        ip=i+1;
        if(i==0)
        {
            im=0;
        }
        else if(i==nx-1)
        {
            ip=nx-1;
        }
        for(j=0; j<ny; j++)
        {
            jm=j-1;
            jp=j+1;
            if(j==0)
            {
                jm=0;
            }
            else if(j==ny-1)
            {
                jp=ny-1;
            }

            f0[i][j]=f0aux[i][j];
            f1[i][j]=f1aux[im][j];
            f2[i][j]=f2aux[i][jm];
            f3[i][j]=f3aux[ip][j];
            f4[i][j]=f4aux[i][jp];

        }
    }


}



void inicializa()
{
    int i, j;
    double v = 1.0;
    //PESOS
    double w0 = (2.0/6.0), w = (1.0/6.0);


    for(i=(nx/4); i<(nx-(nx/4)); i++)
    {
        for(j=(ny/4); j<(ny-(ny/4)); j++)
        {

            f0[i][j]= w0*(v); //parado VERIFICAR!!!

            f1[i][j]= w*(v);
            f2[i][j]= w*(v);
            f3[i][j]= w*(v);
            f4[i][j]= w*(v);



        }
    }
}

void save_vtk(char *nome, double t[][ny])
{
    ofstream arqvtk;
    arqvtk.open(nome);
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " double\n";
    for(int i=0; i<nx; i++)
        arqvtk << i << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " double\n";
    for(int j=0; j<ny; j++)
        arqvtk << j << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 double\n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx*ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
    arqvtk << "VelocityMagnitude 1 " << nx*ny << " double\n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
        {
            arqvtk << t[i][j] << " ";
        }

    arqvtk << "\n";


    arqvtk.close();
}

void velocidade(char *salvando){
    int j;
    ofstream VELOCIDADE;
    VELOCIDADE.open(salvando);
    cout<< "Valores da Velocidade" << endl;

    for(j=0;j<ny;j++){

        VELOCIDADE << T[nx/2][j] << endl;

    }
    VELOCIDADE.close();

}


int main(int argc, char** argv)
{
    setlocale(LC_ALL, "Portuguese");
    int t=0;
    inicializa();


    char nome[150];



    while(t<1000)
    {
    	cout << "t " << t << std::endl;

    	cout << "Colisão" << std::endl;
	    colisao();
	    cout << "Propagação \n" << std::endl;
	    propagacao();
	    cout << "Contorno" << std::endl;
	    paredes();

	    if(t%10==0)
	    {
	    	sprintf(nome, "resultado_%d.vtk", t/10);
	    	save_vtk(nome,T);

	    }

	    t++;

    }
  char texto[25] = "VELOCIDADE.txt";
  velocidade(texto);

}
