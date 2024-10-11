#include<iostream>
#include<fstream>
#define nx 512
#define ny 512

using namespace std;
double f0[nx][ny]= {0},f1[nx][ny]= {0},f2[nx][ny]= {0},f3[nx][ny]= {0},f4[nx][ny]= {0};
double f0aux[nx][ny]= {0},f1aux[nx][ny]= {0},f2aux[nx][ny]= {0},f3aux[nx][ny]= {0},f4aux[nx][ny]= {0};
double f0eq[nx][ny]= {0},f1eq[nx][ny]= {0},f2eq[nx][ny]= {0},f3eq[nx][ny]= {0},f4eq[nx][ny]= {0};

double V[nx][ny]= {0};

double D=2, dt=0.1, dx=0.02, alpha = 0.0003;
double inv_omega   = (dt*D*alpha/(dx*dx) + 0.5);
double tau = dt*inv_omega;
double omega = 1./inv_omega;


double h[nx][ny]= {0};
double t=0;


double dvdt(int i, int j,double t, double v, double h)
{
    double Tin=0.3, Tout = 6.0;
    double ti=0.1;
    double tf=1.1;
    double Istim;

    double ti1 = 380;
    double tf1 = 383;


    if(((t>=ti)&&(t<=tf)&& (i>0) && (i<20) && (j>0) && (j<512)) || ((t>=ti1)&&(t<=tf1)&&(i>0)&&(i<512)&&(j>0)&&(j<100)))
        Istim = 0.1;
    else
        Istim=0;

    double C = v*v*(1.0 -v);
    double Iin= h*C/Tin;
    double Iout = -v/Tout;

    return Iin + Iout + Istim;

}

double dhdt(int i, int j,double t, double v, double h)
{
    double Topen=120.0, Tclose = 150.0, Vgate = 0.13;

    if(v<Vgate)
        return(1.0 - h)/Topen;
    else
        return -h/Tclose;

}

void colisao()
{
    int i, j;
    double v;
    double H;


    double w0 = (2./6.), w = (1./6.);

    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {

            v = f0[i][j] + f2[i][j] + f4[i][j] + f1[i][j] + f3[i][j];
            H = h[i][j];


            V[i][j] = v;

            f0eq[i][j]= w0*(v);

            f1eq[i][j]= w*(v);
            f2eq[i][j]= w*(v);
            f3eq[i][j]= w*(v);
            f4eq[i][j]= w*(v);



            f0aux[i][j]=(f0[i][j]*(1.-omega))+ omega*f0eq[i][j]+ w0*dt*dvdt(i,j,t,v,H);
            f1aux[i][j]=(f1[i][j]*(1.-omega))+ omega*f1eq[i][j]+ w*dt*dvdt(i,j,t,v,H);
            f2aux[i][j]=(f2[i][j]*(1.-omega))+ omega*f2eq[i][j]+ w*dt*dvdt(i,j,t,v,H);
            f3aux[i][j]=(f3[i][j]*(1.-omega))+ omega*f3eq[i][j]+ w*dt*dvdt(i,j,t,v,H);
            f4aux[i][j]=(f4[i][j]*(1.-omega))+ omega*f4eq[i][j]+ w*dt*dvdt(i,j,t,v,H);


            h[i][j] = H + dt*dhdt(i,j,t,v,H);



        }



    }
    t+=dt;


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
    //PESOS
    double w0 = (2./6.), w = (1./6.);
    double v=0;


    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {

            h[i][j] = 1.;


            f0[i][j]= w0*(v); //parado VERIFICAR!!!

            f1[i][j]= w*(v);
            f2[i][j]= w*(v);
            f3[i][j]= w*(v);
            f4[i][j]= w*(v);



        }
    }
}

void save_vtk(char *nome, double d[][ny]) //d=densidade
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
    arqvtk << "Potential 1 " << nx*ny << " double\n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
        {
            arqvtk << d[i][j] << " ";
        }

    arqvtk << "\n";


    arqvtk.close();
}

int main(int argc, char** argv)
{

    int TempoSim=0;
    inicializa();


    char nome[150];



    while(TempoSim<1000/dt) // 10 000 interações
    {
        cout << "tempo" << TempoSim << std::endl;

        cout << "Colisão" << std::endl;
        colisao();
        cout << "Propagação \n" << std::endl;
        propagacao();
        cout << "Contorno" << std::endl;
        paredes();

        if(TempoSim%100==0)// relativo ao "skip = 100" do WebGL
        {
            sprintf(nome, "resultado_%d.vtk", TempoSim/100);
            save_vtk(nome,V);

        }

        TempoSim++;

    }
}
