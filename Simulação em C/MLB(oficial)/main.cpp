#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#define nx 256
#define ny 64

using namespace std;

double altura,largura;
double f2[nx][ny] = {0}, f1[nx][ny] = {0}, f4[nx][ny] = {0}, f3[nx][ny] = {0};
double f2aux[nx][ny]= {0}, f1aux[nx][ny]= {0}, f4aux[nx][ny]= {0}, f3aux[nx][ny]= {0};
double f5[nx][ny]= {0},f8[nx][ny]= {0},f6[nx][ny]= {0},f7[nx][ny]= {0}, f0[nx][ny]= {0};
double f5aux[nx][ny]= {0}, f8aux[nx][ny]= {0}, f6aux[nx][ny]= {0},f7aux[nx][ny]= {0},f0aux[nx][ny]= {0}; //f1eções diagonais(matriz diagonal)
double f0eq[nx][ny]= {0}, f2eq[nx][ny]= {0}, f4eq[nx][ny]= {0},f1eq[nx][ny]= {0},f3eq[nx][ny]= {0},f5eq[nx][ny]= {0},f8eq[nx][ny]= {0},f7eq[nx][ny]= {0},f6eq[nx][ny]= {0};


double V[nx][ny]= {0};
double Ux[nx][ny]= {0}, Uy[nx][ny]= {0};

double tau    = 1.0; //3.*kvisc + 0.5;
double umax = 0.1;
double lbD = ny-1;
double Re = 100;
double kvisc  = (1.0/3.0)*(tau-0.5); // viscosidade??
double forcex = (8.0 * umax * kvisc)/((lbD)*(lbD));
double forcey = 0.0;


void colisao()
{
    int i, j;
    double densidade, x,y, ux2uy2;
    //Pesos
    double w0 = (4.0/9.0), w = (1.0/9.0), wdiagonal = (1.0/36.0);



    // COLISÃO
    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            densidade = f0[i][j] + f2[i][j]+ f4[i][j]+f1[i][j]+f3[i][j]+f6[i][j]+f8[i][j]+f5[i][j]+f7[i][j];

            //decompor as forças ????
            // onde acrecentar a f1eção???

            x = (f1[i][j]-f3[i][j]+f5[i][j]+f8[i][j]-f6[i][j]-f7[i][j])/(densidade);
            y = (f2[i][j]-f4[i][j]+f5[i][j]+f6[i][j]-f7[i][j]-f8[i][j])/(densidade);

            Ux[i][j]= x;
            Uy[i][j]= y;

            V[i][j] = sqrt(x*x + y*y);


            ux2uy2 = x*x+y*y;

            f0eq[i][j]= w0*(densidade)*(1.0-(3.0/2.0)*(ux2uy2)); //parado VERIFICAR!!!

            f1eq[i][j]= w*(densidade)*(1.0+3.0*(x)+(9.0/2.0)*(x*x)-(3.0/2.0)*ux2uy2);
            f2eq[i][j]= w*(densidade)*(1.0+3.0*(y)+(9.0/2.0)*(y*y)-(3.0/2.0)*ux2uy2);
            f3eq[i][j]= w*(densidade)*(1.0+3.0*(-x)+(9.0/2.0)*(x*x)-(3.0/2.0)*ux2uy2);
            f4eq[i][j]= w*(densidade)*(1.0+3.0*(-y)+(9.0/2.0)*(y*y)-(3.0/2.0)*ux2uy2);


            f5eq[i][j]= wdiagonal*(densidade)*(1.0+3.0*(x+y)+(9.0/2.0)*((x+y)*(x+y))-(3.0/2.0)*ux2uy2);
            f6eq[i][j]= wdiagonal*(densidade)*(1.0+3.0*(-x+y)+(9.0/2.0)*((-x+y)*(-x+y))-(3.0/2.0)*ux2uy2);
            f7eq[i][j]= wdiagonal*(densidade)*(1.0+3.0*(-x-y)+(9.0/2.0)*((-x-y)*(-x-y))-(3.0/2.0)*ux2uy2);
            f8eq[i][j]= wdiagonal*(densidade)*(1.0+3.0*(x-y)+(9.0/2.0)*((x-y)*(x-y))-(3.0/2.0)*ux2uy2);


            f0aux[i][j]=(1./tau)* f0eq[i][j] + (1.-1./tau)*f0[i][j];
            f1aux[i][j]=(1./tau)* f1eq[i][j] + (1.-1./tau)*f1[i][j] + forcex;
            f2aux[i][j]=(1./tau)* f2eq[i][j] + (1.-1./tau)*f2[i][j];
            f3aux[i][j]=(1./tau)* f3eq[i][j] + (1.-1./tau)*f3[i][j] - forcex;
            f4aux[i][j]=(1./tau)* f4eq[i][j] + (1.-1./tau)*f4[i][j];
            f5aux[i][j]=(1./tau)* f5eq[i][j] + (1.-1./tau)*f5[i][j] + forcex;
            f6aux[i][j]=(1./tau)* f6eq[i][j] + (1.-1./tau)*f6[i][j] - forcex;
            f7aux[i][j]=(1./tau)* f7eq[i][j] + (1.-1./tau)*f7[i][j] - forcex;
            f8aux[i][j]=(1./tau)* f8eq[i][j] + (1.-1./tau)*f8[i][j] + forcex;

        }


    }
}

void paredes()
{
     for(int i=0; i<nx; i++)
    {
        // cima
        f4[i][ny-1] = f2[i][ny-1];
        f7[i][ny-1] = f5[i][ny-1];
        f8[i][ny-1] = f6[i][ny-1];

        // baixo
        f2[i][0] = f4[i][0];
        f5[i][0] = f7[i][0];
        f6[i][0] = f8[i][0];
    }

    for(int j=0; j<ny; j++)
    {
        //Esquerda
        f1[0][j] = f1[nx-1][j];
        f5[0][j] = f5[nx-1][j];
        f8[0][j] = f8[nx-1][j];

        //Direita
        f3[nx-1][j] = f3[0][j];
        f6[nx-1][j] = f6[0][j];
        f7[nx-1][j] = f7[0][j];
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
            f5[i][j]=f5aux[im][jm];
            f6[i][j]=f6aux[ip][jm];
            f7[i][j]=f7aux[ip][jp];
            f8[i][j]=f8aux[im][jp];
        }
    }


}



void inicializa()
{
    int i, j;
    double densidade = 1.0, ux2uy2;
    //PESOS
    double w0 = (4.0/9.0), w = (1.0/9.0), wdiagonal = (1.0/36.0);
    double Ux = 0.0, Uy = 0.0;

    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            ux2uy2 = Ux*Ux+Uy*Uy;

            f0[i][j]= w0*(densidade)*(1.0-(3.0/2.0)*(ux2uy2)); //parado VERIFICAR!!!

            f1[i][j]= w*(densidade)*(1.0+3.0*(Ux)+(9.0/2.0)*(Ux*Ux)-(3.0/2.0)*ux2uy2);
            f2[i][j]= w*(densidade)*(1.0+3.0*(Uy)+(9.0/2.0)*(Uy*Uy)-(3.0/2.0)*ux2uy2);
            f3[i][j]= w*(densidade)*(1.0+3.0*(-Ux)+(9.0/2.0)*(Ux*Ux)-(3.0/2.0)*ux2uy2);
            f4[i][j]= w*(densidade)*(1.0+3.0*(-Uy)+(9.0/2.0)*(Uy*Uy)-(3.0/2.0)*ux2uy2);


            f5[i][j]= wdiagonal*(densidade)*(1.0+3.0*(Ux+Uy)+(9.0/2.0)*((Ux+Uy)*(Ux+Uy))-(3.0/2.0)*ux2uy2);
            f6[i][j]= wdiagonal*(densidade)*(1.0+3.0*(-Ux+Uy)+(9.0/2.0)*((-Ux+Uy)*(-Ux+Uy))-(3.0/2.0)*ux2uy2);
            f7[i][j]= wdiagonal*(densidade)*(1.0+3.0*(-Ux-Uy)+(9.0/2.0)*((-Ux-Uy)*(-Ux-Uy))-(3.0/2.0)*ux2uy2);
            f8[i][j]= wdiagonal*(densidade)*(1.0+3.0*(Ux-Uy)+(9.0/2.0)*((Ux-Uy)*(Ux-Uy))-(3.0/2.0)*ux2uy2);

        }
    }
}

void save_vtk(char *nome, double v[][ny], double vx[][ny], double vy[][ny])
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
            arqvtk << v[i][j] << " ";
        }

    arqvtk << "\n";

    arqvtk << "VECTORS vectors double\n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
        {
            arqvtk << vx[i][j] << " "<< vy[i][j] << " "<< 0.0 << "    ";
        }


    arqvtk << "\n";
    arqvtk.close();
}

void velocidade(char *salvando){
    int j;
    ofstream velocity;
    velocity.open(salvando);
    cout<< "Valores da Velocidade" << endl;

    for(j=0;j<ny;j++){

        velocity << V[nx/2][j] << endl;

    }
    velocity.close();

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
	    	save_vtk(nome, V, Ux, Uy);

	    }

	    t++;

    }

    velocidade("velocity.txt");

}
