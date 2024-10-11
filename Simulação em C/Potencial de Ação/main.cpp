#include <iostream>
#include <fstream>
#define TAM 100
using namespace std;

float dvdt(float t, float v, float h){
    float Tin=0.3, Tout = 6.0;
    float ti=0.1;
    float tf=1;
    float Istim;

    if((t>=ti)&&(t<=tf))
            Istim = 0.1;
        else
            Istim=0;

            float C = v*v*(1.0 -v);
            float Iin= h*C/Tin;
            float Iout = -v/Tout;

            return Iin + Iout + Istim;

}

float dhdt(float t, float v, float h){
    float Topen=120.0, Tclose = 150.0, Vgate = 0.13;

    if(v<Vgate)
        return(1.0 - h)/Topen;
    else
        return -h/Tclose;

}

int main()
{

    ofstream arqSaida;
    arqSaida.open("saida.txt");
    float tempo_sim = 500; //tempo de simulação em segundos
    float passo=0.1;
    int total_passos = int(tempo_sim/passo);
    float u[total_passos];
    u[0]=0;
    float t[total_passos];
    t[0]=0;
    float h[total_passos]; //abertura e fechamento dos canais ionicos;
    h[0]=1;

    arqSaida<< t[0] << " " << u[0] << " " << h[0] << " " << endl;


    for(int i=1; i<total_passos; i++){
        t[i]=t[i-1]+passo;
        u[i]=u[i-1]+passo*dvdt(t[i-1],u[i-1],h[i-1]);
        h[i]=h[i-1]+passo*dhdt(t[i-1],u[i-1],h[i-1]);
        cout<< "Valor de tn " << t[i]<<endl;
        cout<< "Valor de v " << u[i] << endl;
        cout<< "Valor de h " << h[i]<<endl;

        arqSaida << t[i] << " " << u[i] << " " << h[i] << " " << endl ;

    }



    arqSaida.close();



    return 0;
}





