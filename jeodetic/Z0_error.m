clc;
clear all;
close all;
format long g;

%ShahabEsfandiar_9819373

%Observed Data

l=[7.9985;15.9903;23.9929;31.9897;39.9973;
   8.0037;15.9957;23.9952;31.9922;
   7.9949;15.9971;23.9893;
   8.0003;15.9939
   7.9899];
   
%--------------------------------------------------------------------------

%compute z0 error with linear parametric method
syms s12 s23 s34 s45 s56 z01;

L1(1)=s12-z01;
L1(2)=s12+s23-z01;
L1(3)=s12+s23+s34-z01;
L1(4)=s12+s23+s34+s45-z01;
L1(5)=s12+s23+s34+s45+s56-z01;
L1(6)=s23-z01;
L1(7)=s23+s34-z01;
L1(8)=s23+s34+s45-z01;
L1(9)=s23+s34+s45+s56-z01;
L1(10)=s34-z01;
L1(11)=s34+s45-z01;
L1(12)=s34+s45+s56-z01;
L1(13)=s45-z01;
L1(14)=s45+s56-z01;
L1(15)=s56-z01;

A1=jacobian(L1,[s12 s23 s34 s45 s56 z01]);

l1=l;
xcap=inv(A1'*A1)*A1'*l1;
xcap1=[7034731103099988899/879609302220800000,
      42226506570225681149/5277655813324800000,
      21102522051625746277/2638827906662400000,
      28152235833173563851/3518437208883200000,
      84398829207434569583/10555311626649600000,
      8760908650121453/3518437208883200000];

vcap1=A1*xcap1-l1;
lcap1=l1+vcap1;
df=15-6;
cov=(vcap1'*eye*vcap1)/df;

%--------------------------------------------------------------------------

%compute z0 error with non_linear parametric method
syms X2 X3 X4 X5 X6 Y2 Y3 Y4 Y5 Y6 z02;

m=0;
L2(1)=sqrt((X2-0)^2+(Y2-0)^2)-z02;
L2(2)=sqrt((X3-0)^2+(Y3-0)^2)-z02;
L2(3)=sqrt((X4-0)^2+(Y4-0)^2)-z02;
L2(4)=sqrt((X5-0)^2+(Y5-0)^2)-z02;
L2(5)=sqrt((X6-0)^2+(Y6-0)^2)-z02;
L2(6)=sqrt((X3-X2)^2+(Y3-Y2)^2)-z02;
L2(7)=sqrt((X4-X2)^2+(Y4-Y2)^2)-z02;
L2(8)=sqrt((X5-X2)^2+(Y5-Y2)^2)-z02;
L2(9)=sqrt((X6-X2)^2+(Y6-Y2)^2)-z02;
L2(10)=sqrt((X4-X3)^2+(Y4-Y3)^2)-z02;
L2(11)=sqrt((X5-X3)^2+(Y5-Y3)^2)-z02;
L2(12)=sqrt((X6-X3)^2+(Y6-Y3)^2)-z02;
L2(13)=sqrt((X5-X4)^2+(Y5-Y4)^2)-z02;
L2(14)=sqrt((X6-X4)^2+(Y6-Y4)^2)-z02;
L2(15)=sqrt((X6-X5)^2+(Y6-Y5)^2)-z02;
L2(16)=m*X2+Y2;
L2(17)=m*X3+Y3;
L2(18)=m*X4+Y4;
L2(19)=m*X5+Y5;
L2(20)=m*X6+Y6;

A2=jacobian(L2,[X2 X3 X4 X5 X6 Y2 Y3 Y4 Y5 Y6 z02]);

X2=8;
X3=16;
X4=24;
X5=32;
X6=40;
Y2=0;
Y3=0;
Y4=0;
Y5=0;
Y6=0;
z02=0;

ln=[7.9985;15.9903;23.9929;31.9897;39.9973;
   8.0037;15.9957;23.9952;31.9922;
   7.9949;15.9971;23.9893;
   8.0003;15.9939
   7.9899;
   0;
   0;
   0;
   0;
   0;];

for i=1:100
    
    l2=eval(L2);
    dl2=ln-l2';
    
    A2n=eval(A2);
    dx2=inv(A2n'*A2n)*A2n'*dl2;
    
    X2=X2+dx2(1);
    X3=X3+dx2(2);
    X4=X4+dx2(3);
    X5=X5+dx2(4);
    X6=X6+dx2(5);
    Y2=Y2+dx2(6);
    Y3=Y3+dx2(7);
    Y4=Y4+dx2(8);
    Y5=Y5+dx2(9);
    Y6=Y6+dx2(10);
    z02=z02+dx2(11);
    
    if abs(dx2)<10^-3
        break
    end
    
end

xcap2=[X2;X3;X4;X5;X6;z02];
 
%--------------------------------------------------------------------------

%compute periodic error with non_linear parametric method
syms sn12 sn23 sn34 sn45 sn56 z03 A T phi;

L3(1)=sn12-z03-(A*cosd((2*pi*l(1)/T)+phi));
L3(2)=sn12+sn23-z03-(A*cosd((2*pi*l(2)/T)+phi));
L3(3)=sn12+sn23+sn34-z03-(A*cosd((2*pi*l(3)/T)+phi));
L3(4)=sn12+sn23+sn34+sn45-z03-(A*cosd((2*pi*l(4)/T)+phi));
L3(5)=sn12+sn23+sn34+sn45+sn56-z03-(A*cosd((2*pi*l(5)/T)+phi));
L3(6)=sn23-z03-(A*cosd((2*pi*l(6)/T)+phi));
L3(7)=sn23+sn34-z03-(A*cosd((2*pi*l(7)/T)+phi));
L3(8)=sn23+sn34+sn45-z03-(A*cosd((2*pi*l(8)/T)+phi));
L3(9)=sn23+sn34+sn45+sn56-z03-(A*cosd((2*pi*l(9)/T)+phi));
L3(10)=sn34-z03-(A*cosd((2*pi*l(10)/T)+phi));
L3(11)=sn34+sn45-z03-(A*cosd((2*pi*l(11)/T)+phi));
L3(12)=sn34+sn45+sn56-z03-(A*cosd((2*pi*l(12)/T)+phi));
L3(13)=sn45-z03-(A*cosd((2*pi*l(13)/T)+phi));
L3(14)=sn45+sn56-z03-(A*cosd((2*pi*l(14)/T)+phi));
L3(15)=sn56-z03-(A*cosd((2*pi*l(15)/T)+phi));

A3=jacobian(L3,[sn12 sn23 sn34 sn45 sn56 A T phi z03]);

sn12=8;
sn23=8;
sn34=8;
sn45=8;
sn56=8;
A=0.03;
T=11;
phi=0;
z03=0;

for i=1:100
    
    l3=eval(L3);
    dl3=l-l3';
    
    A3n=eval(A3);
    dx3=pinv(A3n'*A3n)*A3n'*dl3;
    
    sn12=sn12+dx3(1);
    sn23=sn23+dx3(2);
    sn34=sn34+dx3(3);
    sn45=sn45+dx3(4);
    sn56=sn56+dx3(5);
    A=A+dx3(6);
    T=T+dx3(7);
    phi=phi+dx3(8);
    z03=z03+dx3(9);
    
    if abs(dx3)<10^-3
        break
    end
end

xcap3=[sn12;sn23;sn34;sn45;sn56;A;T;phi;z03];

%--------------------------------------------------------------------------

%compute scale error with non_linear parametric method
syms K3 An Tn phin z04;

L4(1)=K3*sn12-z03-(An*cosd(((2*pi*l(1)/Tn)+phin)));
L4(2)=K3*(sn12+sn23)-z03-(An*cosd(((2*pi*l(2)/Tn)+phin)));
L4(3)=K3*(sn12+sn23+sn34)-z03-(An*cosd(((2*pi*l(3)/Tn)+phin)));
L4(4)=K3*(sn12+sn23+sn34+sn45)-z03-(An*cosd(((2*pi*l(4)/Tn)+phin)));
L4(5)=K3*(sn12+sn23+sn34+sn45+sn56)-z03-(An*cosd(((2*pi*l(5)/Tn)+phin)));
L4(6)=K3*sn23-z03-(An*cosd(((2*pi*l(6)/Tn)+phin)));
L4(7)=K3*(sn23+sn34)-z03-(An*cosd(((2*pi*l(7)/Tn)+phin)));
L4(8)=K3*(sn23+sn34+sn45)-z03-(An*cosd(((2*pi*l(8)/Tn)+phin)));
L4(9)=K3*(sn23+sn34+sn45+sn56)-z03-(An*cosd(((2*pi*l(9)/Tn)+phin)));
L4(10)=K3*sn34-z03-(An*cosd(((2*pi*l(10)/Tn)+phin)));
L4(11)=K3*(sn34+sn45)-z03-(An*cosd(((2*pi*l(11)/Tn)+phin)));
L4(12)=K3*(sn34+sn45+sn56)-z03-(An*cosd(((2*pi*l(12)/Tn)+phin)));
L4(13)=K3*sn45-z03-(An*cosd(((2*pi*l(13)/Tn)+phin)));
L4(14)=K3*(sn45+sn56)-z03-(An*cosd(((2*pi*l(14)/Tn)+phin)));
L4(15)=K3*sn56-z03-(An*cosd(((2*pi*l(15)/Tn)+phin)));

A4=jacobian(L4,[K3 An Tn phin z04]);

K3=1;
An=0.03;
Tn=11;
phin=0;
z04=0;

for i=1:100
    
    l4=eval(L4);
    dl4=l-l4';
    
    A4n=eval(A4);
    dx4=pinv(A4n'*A4n)*A4n'*dl4;
    
    K3=K3+dx4(1);
    An=An+dx4(2);
    Tn=Tn+dx4(3);
    phin=phin+dx4(4);
    z04=z04+dx4(5);
    
    if abs(dx4)<10^-3
        break
    end
end

xcap4=[K3;A;T;phi;z03];
vcap4=A4n*dx4-dl4;
lcap4=l+vcap4;

%--------------------------------------------------------------------------

%dispaly unknown parameteres

disp('distance between staitions and z0 :')
disp(xcap1)
disp('----------------------------------------')

disp('X coordinate staitions and z0 :')
disp(xcap2)
disp('----------------------------------------')

disp('distance between staitions, z0 and periodic error parameteres :')
disp(xcap3)
disp('----------------------------------------')

disp('scale error, distance between staitions, z0 and periodic error parameteres :')
disp(xcap4)

%end

