#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#define nx 50
#define ny 50
#define NEQ 8

using namespace std;

typedef struct
{
    float sv[8];
} celula;

float f0[nx][ny] = {0},f1[nx][ny] = {0},f2[nx][ny] = {0},f3[nx][ny] = {0},f4[nx][ny] = {0};
float f0aux[nx][ny] = {0},f1aux[nx][ny] = {0},f2aux[nx][ny] = {0},f3aux[nx][ny] = {0},f4aux[nx][ny] = {0};
float f0eq[nx][ny] = {0},f1eq[nx][ny] = {0},f2eq[nx][ny] = {0},f3eq[nx][ny] = {0},f4eq[nx][ny] = {0};

float D=3., dt=0.01, dx=0.02, alpha = 0.0003;
float inv_omega   = (dt*D*alpha/(dx*dx) + 0.5);
float tau = dt*inv_omega;
float omega = 1./inv_omega;

//float h[nx][ny];
celula V[nx][ny]= {0};

float t=0;

void RHS_cpu(float time_, float *sv_, float *rDY_, int i, int j)
{
    // State variables
    const float V_old_ = sv_[0];	 // initial value = -83.853 millivolt
    const float m_old_ = sv_[1];	 // initial value = 0.00187018 dimensionless
    const float h_old_ = sv_[2];	 // initial value = 0.9804713 dimensionless
    const float j_old_ = sv_[3];	 // initial value = 0.98767124 dimensionless
    const float d_old_ = sv_[4];	 // initial value = 0.00316354 dimensionless
    const float f_old_ = sv_[5];	 // initial value = 0.99427859 dimensionless
    const float X_old_ = sv_[6];	 // initial value = 0.16647703 dimensionless
    const float Cai_old_ = sv_[7];	 // initial value = 0.0002 millimolar


	// Parameters
	//const float stim_amplitude = -25.5;	 // microA_per_cm2

	const float C = 1;	 // microF_per_cm2
	const float R = 8314;	 // joule_per_kilomole_kelvin
	const float T = 310;	 // kelvin
	const float F = 96484.6;	 // coulomb_per_mole
	const float Nao = 140;	 // millimolar
	const float Nai = 18;	 // millimolar
	const float g_Na = 23;	 // milliS_per_cm2
	const float Ko = 5.4;	 // millimolar
	const float PR_NaK = 0.01833;	 // dimensionless
	const float Ki = 145;	 // millimolar
	const float g_Kp = 0.0183;	 // milliS_per_cm2
	const float g_b = 0.03921;	 // milliS_per_cm2
	const float E_b = -59.87;	 // millivolt
    float calc_I_stim=0.0;
	// Independent Variable
	float time_new = time_;

//     for(int i=0; i<1; i++)
//      {
    /*
        istart=estimulos[mapij3(0,i,7)];
        isend=estimulos[mapij3(1,i,7)];
        istimamp=estimulos[mapij3(2,i,7)];
        xini=estimulos[mapij3(3,i,7)];
        xfim=estimulos[mapij3(4,i,7)];
        yini=estimulos[mapij3(5,i,7)];
        yfim=estimulos[mapij3(6,i,7)];
    */
    //istart[0]=0.1; isend[0]=4.0; xini[0]=0; xfim[0]=10; yini[0]=0; yfim[0]=200; istimamp[0]=-25.5;
	const float stim_amplitude = -25.5;	 // microA_per_cm2
	const float stim_start = 100;	 // millisecond
	const float stim_end = 100000000000;	 // millisecond
	const float stim_period = 1000;	 // millisecond
	const float stim_duration = 2;	 // millisecond
      //  if((cellID/neq)%nx>=xini[i] && (cellID/neq)%nx<=xfim[i] && (cellID/neq)/nx>=yini[i] && (cellID/neq)/nx <= yfim[i] && time_ >= istart[i] && time_ <= isend[i])
      float iistart=0.1, iisend=4.0, ixini=0, ixfim=10, iyini=0, iyfim=(ny-1), iistimamp= -25.5;

        if(i>=ixini && i<=ixfim && j>=iyini && j<= iyfim && time_ >= iistart && time_ <= iisend)
        {
          calc_I_stim = iistimamp;
//          break;
        }
        else
          calc_I_stim = 0.0;
//      }

float calc_alpha_h, calc_beta_h, calc_alpha_j, calc_beta_j, calc_Xi;

    float calc_E_Na = (((R*T)/F)*log((Nao/Nai)));	//2
	float calc_alpha_m = (V_old_ !=  -47.13)? ((0.32*(V_old_+47.13))/(1.0-exp(((-0.1)*(V_old_+47.13))))): 3.2;	//4
	float calc_beta_m = (0.08*exp(((-V_old_)/11.0)));	//5

    if((V_old_<(-4.0e+01))) calc_alpha_h = 1.35e-01*exp(((8.0e+01+V_old_)/(-6.8e+00)));
    else calc_alpha_h = 0.0e+00;

//	float calc_alpha_h = ((V_old_<(-40.0)))? (0.135*exp(((80.0+V_old_)/(-6.8)))): 0.0;	//7
    if((V_old_<(-4.0e+01))) calc_beta_h = (3.56e+00*exp((7.9e-02*V_old_)))+(3.1e+05*exp((3.5e-01*V_old_)));
    else calc_beta_h = 1.0e+00/(1.3e-01*(1.0e+00+exp(((V_old_+1.066e+01)/(-1.11e+01)))));

//	float calc_beta_h = ((V_old_<(-40.0)))? ((3.56*exp((0.079*V_old_)))+(310000.0*exp((0.35*V_old_)))): (1.0/(0.13*(1.0+exp(((V_old_+10.66)/(-11.1))))));	//8
	if((V_old_<(-4.0e+01))) calc_alpha_j = ((((-1.2714e+05)*exp((2.444e-01*V_old_)))-(3.474e-05*exp(((-4.391e-02)*V_old_))))*(V_old_+3.778e+01))/(1.0e+00+exp((3.11e-01*(V_old_+7.923e+01))));
    else calc_alpha_j = 0.0e+00;

//	float calc_alpha_j = ((V_old_<(-40.0)))? (((((-127140.0)*exp((0.2444*V_old_)))-(0.00003474*exp(((-0.04391)*V_old_))))*(V_old_+37.78))/(1.0+exp((0.311*(V_old_+79.23))))): 0.0;	//10
    if((V_old_<(-4.0e+01))) calc_beta_j = (1.212e-01*exp(((-1.052e-02)*V_old_)))/(1.0e+00+exp(((-1.378e-01)*(V_old_+4.014e+01))));
    else calc_beta_j = (((3.0e-01*exp(((-2.535e-07)*V_old_)))/(1.0e+00+exp(((-1.0e-01)*(V_old_+3.2e+01))))));
//	float calc_beta_j = ((V_old_<(-40.0)))? ((0.1212*exp(((-0.01052)*V_old_)))/(1.0+exp(((-0.1378)*(V_old_+40.14))))): ((0.3*exp(((-0.0000002535)*V_old_)))/(1.0+exp(((-0.1)*(V_old_+32.0)))));	//11
	float calc_E_si = (7.7-(13.0287*log((Cai_old_/1.0))));	//13
	float calc_alpha_d = ((0.095*exp(((-0.01)*(V_old_-5.0))))/(1.0+exp(((-0.072)*(V_old_-5.0)))));	//15
	float calc_beta_d = ((0.07*exp(((-0.017)*(V_old_+44.0))))/(1.0+exp((0.05*(V_old_+44.0)))));	//16
	float calc_alpha_f = ((0.012*exp(((-0.008)*(V_old_+28.0))))/(1.0+exp((0.15*(V_old_+28.0)))));	//18
	float calc_beta_f = ((0.0065*exp(((-0.02)*(V_old_+30.0))))/(1.0+exp(((-0.2)*(V_old_+30.0)))));	//19
	float calc_g_K = (0.282*pow((Ko/5.4),1.0/2.0));	//21
	float calc_E_K = (((R*T)/F)*log(((Ko+(PR_NaK*Nao))/(Ki+(PR_NaK*Nai)))));	//22
	float calc_alpha_X = ((0.0005*exp((0.083*(V_old_+50.0))))/(1.0+exp((0.057*(V_old_+50.0)))));	//24
	float calc_beta_X = ((0.0013*exp(((-0.06)*(V_old_+20.0))))/(1.0+exp(((-0.04)*(V_old_+20.0)))));	//25

    if((V_old_>(-1.0e+02))) calc_Xi = (((2.837e+00*(exp((4.0e-02*(V_old_+7.7e+01)))-1.0e+00))/((V_old_+7.7e+01)*exp((4.0e-02*(V_old_+3.5e+01))))));
     else calc_Xi = 1.0e+00;
//	float calc_Xi = ((V_old_>(-100.0)))? ((V_old_ != -77.0)?((2.837*(exp((0.04*(V_old_+77.0)))-1.0))/((V_old_+77.0)*exp((0.04*(V_old_+35.0))))):0.608883292 ): 1.0;	//27
	float calc_g_K1 = (0.6047*pow((Ko/5.4),1.0/2.0));	//28
	float calc_E_K1 = (((R*T)/F)*log((Ko/Ki)));	//29
	float calc_Kp = (1.0/(1.0+exp(((7.488-V_old_)/5.98))));	//35
	float calc_i_b = (g_b*(V_old_-E_b));	//37
	float calc_i_Na = (g_Na*pow(m_old_,3.0)*h_old_*j_old_*(V_old_-calc_E_Na));	//3
	float calc_i_si = (0.09*d_old_*f_old_*(V_old_-calc_E_si));	//14
	float calc_alpha_K1 = (1.02/(1.0+exp((0.2385*((V_old_-calc_E_K1)-59.215)))));	//31
	float calc_beta_K1 = (((0.49124*exp((0.08032*((V_old_+5.476)-calc_E_K1))))+(1.0*exp((0.06175*(V_old_-(calc_E_K1+594.31))))))/(1.0+exp(((-0.5143)*((V_old_-calc_E_K1)+4.753)))));	//32
	float calc_E_Kp = calc_E_K1;	//34
	float calc_i_K = (calc_g_K*X_old_*calc_Xi*(V_old_-calc_E_K));	//23
	float calc_K1_infinity = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));	//33
	float calc_i_Kp = (g_Kp*calc_Kp*(V_old_-calc_E_Kp));	//36
	float calc_i_K1 = (calc_g_K1*calc_K1_infinity*(V_old_-calc_E_K1));	//30

	// Differential Equations
	float d_dt_V = (((-1.0)/C)*(calc_I_stim+calc_i_Na+calc_i_si+calc_i_K+calc_i_K1+calc_i_Kp+calc_i_b));	// 1
	float d_dt_m = ((calc_alpha_m*(1.0-m_old_))-(calc_beta_m*m_old_));	// 6
	float d_dt_h = ((calc_alpha_h*(1.0-h_old_))-(calc_beta_h*h_old_));	// 9
	float d_dt_j = ((calc_alpha_j*(1.0-j_old_))-(calc_beta_j*j_old_));	// 12
	float d_dt_d = ((calc_alpha_d*(1.0-d_old_))-(calc_beta_d*d_old_));	// 17
	float d_dt_f = ((calc_alpha_f*(1.0-f_old_))-(calc_beta_f*f_old_));	// 20
	float d_dt_X = ((calc_alpha_X*(1.0-X_old_))-(calc_beta_X*X_old_));	// 26
	float d_dt_Cai = ((((-0.0001)/1.0)*calc_i_si)+(0.07*(0.0001-Cai_old_)));	// 38

	rDY_[0] = d_dt_V;
	rDY_[1] = d_dt_m;
	rDY_[2] = d_dt_h;
	rDY_[3] = d_dt_j;
	rDY_[4] = d_dt_d;
	rDY_[5] = d_dt_f;
	rDY_[6] = d_dt_X;
	rDY_[7] = d_dt_Cai;

}

//-----------------------------------------------------------------------------------------------------------
void solve_Forward_Euler_cpu(float time, celula *cell, float dt, int i,int j)
{
	// sv contains the membrane state variables (NEQ) of all cells in the 'tissue' (ni*nj)
	// cellID defines where the state variable starts in sv
 //printf("safe 2\n");
	float rY[NEQ], rDY[NEQ];

	for(int a = 0; a < NEQ; a++)
		rY[a] = cell->sv[a];
 //printf("safe 3\n");
	RHS_cpu(time, rY, rDY, i ,j);

	for(int a = 0; a < NEQ; a++)
		cell->sv[a] = dt*rDY[a] + rY[a];

    //cout<<"V  "<<sv[0]<<endl;

}


void solve_ode_cpu(float time, float dt, celula *cell, int i, int j)
{
    //printf("safe 1\n");
    solve_Forward_Euler_cpu(time, cell, dt, i, j);
}
//-----------------------------------------------------------------------------------------------------------

void colisao()
{
    int i, j;
    float v, v_old;
    //Pesos
    float w0 = (2.0/6.0), w = (1.0/6.0);

    float passo = 0.01;


    // COLISÃO
    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {
            v = f0[i][j] + f2[i][j] + f4[i][j] + f1[i][j] + f3[i][j];
           // h_1 = h[i][j];



            V[i][j].sv[0] = v;



            f0eq[i][j]= w0*(v);

            f1eq[i][j]= w*(v);
            f2eq[i][j]= w*(v);
            f3eq[i][j]= w*(v);
            f4eq[i][j]= w*(v);


            v_old = v;

            solve_ode_cpu(t, passo, &V[i][j] ,i,j);

            v = V[i][j].sv[0];


            f0aux[i][j]=(f0[i][j]*(1.-omega))+ omega*f0eq[i][j] + w0*(v - v_old);// w0*passo*dvdt(i,j,t,v,h_1);
            f1aux[i][j]=(f1[i][j]*(1.-omega))+ omega*f1eq[i][j] + w*(v - v_old);//w*passo*dvdt(i,j,t,v,h_1);
            f2aux[i][j]=(f2[i][j]*(1.-omega))+ omega*f2eq[i][j] + w*(v - v_old);//w*passo*dvdt(i,j,t,v,h_1);
            f3aux[i][j]=(f3[i][j]*(1.-omega))+ omega*f3eq[i][j] + w*(v - v_old);//w*passo*dvdt(i,j,t,v,h_1);
            f4aux[i][j]=(f4[i][j]*(1.-omega))+ omega*f4eq[i][j] + w*(v - v_old);//w*passo*dvdt(i,j,t,v,h_1);


          //  h[i][j] =  h_1+ passo*dhdt(i,j,t,v,h_1);


        }
    }
    t += passo;
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
//-------------------------------------------------------------------------------------------------------------------------------------

void init_st_var(celula *v_)
{

    float IC[8];

    IC[0] = -83.853;	 // V millivolt
    IC[1] = 0.00187018;	 // m dimensionless
    IC[2] = 0.9804713;	 // h dimensionless
    IC[3] = 0.98767124;	 // j dimensionless
    IC[4] = 0.00316354;	 // d dimensionless
    IC[5] = 0.99427859;	 // f dimensionless
    IC[6] = 0.16647703;	 // X dimensionless
    IC[7] = 0.0002;	 // Cai millimolar

    for(int i=0; i<8; i++)
    {
        v_->sv[i] = IC[i];
    }
}


void inicializa()
{
    int i, j;
    //PESOS
    float w0 = (2.0/6.0), w = (1.0/6.0);


    for(i=0; i<nx; i++)
    {
        for(j=0; j<ny; j++)
        {

            // h[i][j] = 1.;
            init_st_var(&V[i][j]);


            f0[i][j]= w0 * V[i][j].sv[0];

            f1[i][j]= w * V[i][j].sv[0];
            f2[i][j]= w * V[i][j].sv[0];
            f3[i][j]= w * V[i][j].sv[0];
            f4[i][j]= w * V[i][j].sv[0];


        }
    }
}

//------------------------------------------------------------------------------------------------------------------------]

void save_vtk(char *nome, celula d[nx][ny])
{
    ofstream arqvtk;
    arqvtk.open(nome);
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " float\n";
    for(int i=0; i<nx; i++)
        arqvtk << i << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " float\n";
    for(int j=0; j<ny; j++)
        arqvtk << j << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 float\n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx*ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
    arqvtk << "VelocityMagnitude 1 " << nx*ny << " float\n";

    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++)
        {
            arqvtk << d[i][j].sv[0] << " ";
        }

    arqvtk << "\n";


    arqvtk.close();
}

/*void velocidade(char *salvando)
{
    int j;
    ofstream VELOCIDADE;
    VELOCIDADE.open(salvando);
    cout<< "Valores da Velocidade" << endl;

    for(j=0; j<ny; j++)
    {

        VELOCIDADE << V[nx/2][j] << endl;

    }
    VELOCIDADE.close();

}*/
//-------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "Portuguese");

    int TempoSim =0;


    inicializa();


    char nome[150];



    while(TempoSim<10000)
    {
        cout << "t " << t << std::endl;
printf("Ponto 1 :");
printf("valor de V: %g", V[0][0].sv[0]);
printf("valor de m: %g", V[0][0].sv[1]);
printf("valor de h: %g", V[0][0].sv[2]);
printf("valor de j: %g\n", V[0][0].sv[3]);


        cout << "Colisão" << std::endl;
        colisao();
        cout << "Propagação \n" << std::endl;
        propagacao();
        cout << "Contorno" << std::endl;
        paredes();

        if(TempoSim%100==0)
        {
            sprintf(nome, "resultado_%d.vtk", TempoSim/10);
            save_vtk(nome,V);

        }

        TempoSim++;

    }
   /* char texto[25] = "VELOCIDADE.txt";
    velocidade(texto);*/

}
