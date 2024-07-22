#include <iomanip>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#define EPSILON 0.001

double func(double x)
{
    return (pow(x,4) + 2*pow(x,3) - 7*pow(x,2) + 3);
}

void bisection(double a, double b)
{
    if (func(a)*func(b) >= 0)
    {
        cout << "  You haven't assumed right a and b\n";
        return;
    }

    double c;

    do
    {
        c = (a+b)/2;

        if (func(c) == 0.0)
            break;

        else if (func(c)*func(a) < 0)
        {
            cout<<setprecision(4);
            cout<<std::fixed;
            cout<<"\t"<<"\t"<<"\t"<<a<<"\t"<<"|"<<"\t"<<b<<"\t"<<"|"<<"\t"<<c<<"\t"<<"|"<<"\t"<<func(c)<<endl;
            b = c;

        }
        else
        {
            cout<<setprecision(4);
            cout<<std::fixed;
            cout<<"\t"<<"\t"<<"\t"<<a<<"\t"<<"|"<<"\t"<<b<<"\t"<<"|"<<"\t"<<c<<"\t"<<"|"<<"\t"<<func(c)<<endl;
            a = c;

        }


    }
    while (!(abs(func(c)) < EPSILON));

    cout<<"                                                                               -------> the ANSWER is "<<c<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;

}
int main()
{
    cout<<endl;
    cout<<"  Question 3"<<endl;
    cout<<endl;
    cout<<"  F(x)= x^4 + 2(x^3) - 7(x^2) + 3 , EPSILON = 0.001"<<endl;
    cout<<endl;
    cout<<"  F'(x)= 4(x^3) + 6(x^2) - 14x"<<endl;
    cout<<endl;
    cout<<"  Root ranges ----> (a1,b1)=[-4,-3.5]    (a2,b2)=[-1,-0.5]    (a3,b3)=[0.5,1]    (a4,b4)=[1.5,2]"<<endl;
    cout<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<endl;

    double a, b;
    for(int i1=1 ; i1 <=1 ; i1++)
    {
        cout<<"  in (a1,b1)=[-4,-3.5] |----> 1. the function is continuous"<<endl;
        cout<<"                       |----> 2. [F(a1).F(b1) < 0] ----> (correct)"<<endl;
        cout<<"                       |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;


        a=-4, b=-3.5;

        cout<<"\t"<<"\t"<<"\t"<<"a"<<"\t"<<"|"<<"\t"<<"b"<<"\t"<<"|"<<"\t"<<"c"<<"\t"<<"|"<<"\t"<<"func(c)"<<endl;
        cout<<"                     -----------------------------------------------------------"<<endl;

        bisection(a,b);
    }

    cout<<endl;

    for(int i2=1 ; i2 <=1 ; i2++)
    {
        cout<<"  in (a2,b2)=[-1,-0.5] |----> 1. the function is continuous"<<endl;
        cout<<"                       |----> 2. [F(a2).F(b2) < 0] ----> (correct)"<<endl;
        cout<<"                       |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;

        a=-1, b=-0.5;

        cout<<"\t"<<"\t"<<"\t"<<"a"<<"\t"<<"|"<<"\t"<<"b"<<"\t"<<"|"<<"\t"<<"c"<<"\t"<<"|"<<"\t"<<"func(c)"<<endl;
        cout<<"                      -----------------------------------------------------------"<<endl;

        bisection(a,b);
    }

    cout<<endl;

    for(int i3=1 ; i3 <=1 ; i3++)
    {
        cout<<"  in (a3,b3)=[0.5,1]  |----> 1. the function is continuous"<<endl;
        cout<<"                      |----> 2. [F(a3).F(b3) < 0] ----> (correct)"<<endl;
        cout<<"                      |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;

        a=0.5, b=1;

        cout<<"\t"<<"\t"<<"\t"<<"a"<<"\t"<<"|"<<"\t"<<"b"<<"\t"<<"|"<<"\t"<<"c"<<"\t"<<"|"<<"\t"<<"func(c)"<<endl;
        cout<<"                      -----------------------------------------------------------"<<endl;

        bisection(a,b);
    }

    cout<<endl;

    for(int i4=1 ; i4 <=1 ; i4++)
    {
        cout<<"  in (a4,b4)=[1.5,2]  |----> 1. the function is continuous"<<endl;
        cout<<"                      |----> 2. [F(a1).F(b1) < 0] ----> (correct)"<<endl;
        cout<<"                      |----> 3. [F'(x) != 0] ----> (correct)"<<endl;
        cout<<endl;
        cout<<endl;

        a=1.5, b=2;

        cout<<"\t"<<"\t"<<"\t"<<"a"<<"\t"<<"|"<<"\t"<<"b"<<"\t"<<"|"<<"\t"<<"c"<<"\t"<<"|"<<"\t"<<"func(c)"<<endl;
        cout<<"                      -----------------------------------------------------------"<<endl;

        bisection(a,b);
    }
    return 0;
}
