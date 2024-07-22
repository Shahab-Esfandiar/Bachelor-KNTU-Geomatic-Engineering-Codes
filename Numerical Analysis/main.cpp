#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#define EPSILON 0.001

double h(double x,double y)
{
    return ((-2*x*x*y - y*y + 22*y + x - 7)/(4*x*y - 1));
}
double k(double x,double y)
{
    return ((-(x*x) - 2*x*y*y + 14*x + y - 11)/(4*x*y - 1));
}
void fixed(double a, double b)
{
    double X0=2.5;
    double Y0=2.5;
    double X1,Y1,z,w;

    do
    {
        X1 = X0 + h(X0,Y0);
        Y1 = Y0 + k(X0,Y0);

        z=fabs(X1-X0);
        w=fabs(Y1-Y0);

            cout<<setprecision(4);
            cout<<std::fixed;
            cout<<"\t"<<"\t"<<"\t"<<X0<<"\t"<<"|"<<"    "<<X1<<"\t"<<"|"<<"\t"<<z<<"\t"<<"        |"<<"\t"<<Y0<<"\t"<<"  |"<<"    "<<Y1<<"\t"<<"|"<<"\t"<<w<<endl;


            X0 = X1;
            Y0 = Y1;



    }
    while (!((z < EPSILON) && (w < EPSILON)));

    cout<<"                                                                                                                            -------> the ANSWER is |----> Xn = "<<X1<<endl;
    cout<<"                                                                                                                                                   |"<<endl;
    cout<<"                                                                                                                                                   |----> Yn = "<<Y1<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;

}

int main()
{
    cout<<endl;
    cout<<"  Question 6 ---->(Taylor method)"<<endl;
    cout<<endl;
    cout<<"  EPSILON = 0.001"<<endl;
    cout<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<endl;

    double a, b;
    for(int i1=1 ; i1 <=1 ; i1++)
    {
        cout<<"  for first root |----> X0 = 2.5"<<endl;
        cout<<"                 |----> Y0 = 2.5"<<endl;
        cout<<endl;
        cout<<endl;


        a=1, b=1;

        cout<<"\t"<<"\t"<<"\t"<<" X(n)"<<"\t"<<"|"<<"    "<<"X(n-1)"<<"\t"<<"|"<<"    |X(n)-X(n-1)|"<<"\t"<<"|"<<"    "<<"    Y(n)"<<"\t"<<"  |"<<"    "<<"Y(n-1)"<<"\t"<<"|"<<"    |Y(n)-Y(n-1)|"<<endl;
        cout<<"                     -----------------------------------------------    |    ----------------------------------------------"<<endl;

        fixed(a,b);
    }

    cout<<endl;

    return 0;
}

