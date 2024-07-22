#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#define EPSILON 0.0001

double Func(double x)
{
    return (((pow(x,3)) + (4*pow(x,2)) - (3*x)) - 7);
}

double Fprime(double x)
{
    return (3*(pow(x,2)) + (8*x) - 3);
}

double Fzegond(double x)
{
    return ((6*x) + 8);
}
void newton(double a, double b)
{

    double X0, X1;

    if (Func(a)*Func(b) >= 0)
    {
        cout << "  You haven't assumed right a and b\n";
        return;
    }

    else if (Func(a)*Fprime(a) > 0)
    {
        X0 = a;
    }
    else
    {
        X0 = b;
    }

    do
    {

        X1 = X0 - (Func(X0)/Fprime(X0));

        cout<<setprecision(4);
        cout<<std::fixed;
        cout<<"\t"<<"\t"<<"\t"<<X1<<"\t"<<"|"<<"\t"<<X0<<"\t"<<"|"<<"\t"<<Func(X1)<<endl;
        if (Func(X1) == 0.0)
        {
            break;
        }
        X0=X1;


    }
    while (!(abs(Func(X1)) < EPSILON));

    cout<<"                                                                -------> the ANSWER is "<<X1<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;

}
int main()
{
    cout<<endl;
    cout<<"  Question 9"<<endl;
    cout<<endl;
    cout<<"  F(x)= x^3 + 4(x^2) - 3x - 7 , EPSILON = 0.0001"<<endl;
    cout<<endl;
    cout<<"  F'(x)= 3(x^2) + 8x - 3"<<endl;
    cout<<endl;
    cout<<"  F''(x)= 6x + 8"<<endl;
    cout<<endl;
    cout<<"  Root ranges ----> (a1,b1)=[-5,-4]    (a2,b2)=[-2,-1]    (a3,b3)=[1,2]"<<endl;
    cout<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<endl;

    double a, b;
    for(int i1=1 ; i1 <=1 ; i1++)
    {
        cout<<"  in (a1,b1)=[-4,-3.5] |----> 1. the function is continuous"<<endl;
        cout<<"                       |----> 2. [F(a1).F(b1) < 0] ----> (correct)"<<endl;
        cout<<"                       |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<"                       |----> 4. [F\"(x) !=0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;


        a=-5, b=-4;

        cout<<"\t"<<"\t"<<"\t"<<"X(n)"<<"\t"<<"|"<<"\t"<<"X(n-1)"<<"\t"<<"|"<<"\t"<<"Func(X(n))"<<endl;
        cout<<"                     ------------------------------------------------"<<endl;

        newton(a,b);
    }

    cout<<endl;

    for(int i2=1 ; i2 <=1 ; i2++)
    {
        cout<<"  in (a2,b2)=[-1,-0.5] |----> 1. the function is continuous"<<endl;
        cout<<"                       |----> 2. [F(a2).F(b2) < 0] ----> (correct)"<<endl;
        cout<<"                       |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;

        a=-2, b=-1;

        cout<<"\t"<<"\t"<<"\t"<<"X(n)"<<"\t"<<"|"<<"\t"<<"X(n-1)"<<"\t"<<"|"<<"\t"<<"Func(X(n))"<<endl;
        cout<<"                      ----------------------------------------------"<<endl;

        newton(a,b);
    }

    cout<<endl;

    for(int i3=1 ; i3 <=1 ; i3++)
    {
        cout<<"  in (a3,b3)=[0.5,1]  |----> 1. the function is continuous"<<endl;
        cout<<"                      |----> 2. [F(a3).F(b3) < 0] ----> (correct)"<<endl;
        cout<<"                      |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;

        a=1, b=2;

        cout<<"\t"<<"\t"<<"\t"<<"X(n)"<<"\t"<<"|"<<"\t"<<"X(n-1)"<<"\t"<<"|"<<"\t"<<"Func(X(n))"<<endl;
        cout<<"                      -----------------------------------------------"<<endl;

        newton(a,b);
    }
    return 0;
}
