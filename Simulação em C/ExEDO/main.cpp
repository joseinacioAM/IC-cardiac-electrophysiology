#include <iostream>
#include<math.h>

using namespace std;

int main()
{
    float u[3];
    u[0]=2;
    float to=0;
    float tf=0.3;
    float h=(tf-to)/3;
    float t[3];
    t[0]=0;

    for(int i=0;i<3;i++){
        t[i+1]=t[i]+h;
        u[i+1]=u[i]+h*(-u[i]+t[i]+2);
        cout<< "Valores de t " << t[i] <<endl;
        cout<< "Valores de u " << u[i] << endl;
    }
    return 0;
}
