clc;
clear all;
close all;
format long g;

%Project1_Geometric Geodazy
%ShahabEsfandiar_9819373
%HomaGangali_9929953

%Exer1
%--------------------------------------------------------------------------

%Data
a=6378137;
f=1/298.2572;
long=51.338325;
lat=35.699737;
h=1187.196;
%--------------------------------------------------------------------------

%Convert geodetic coordinate to Cartesian
b=a*(1-f);
n=a/(sqrt((cosd(lat)^2)+((b^2/a^2)*sind(lat)^2)));

x=(n+h)*cosd(lat)*cosd(long);
y=(n+h)*cosd(lat)*sind(long);
z=(n*(b^2/a^2)+h)*sind(lat);
r=[x y z]';
%--------------------------------------------------------------------------

%Convert computed Cartesian coordinate to Geodetic
landa=atand(y/x);
p=sqrt(x^2+y^2);
e=sqrt((a^2-b^2)/a^2);
dN=1;dH=1;dphi=1;i=2;

N(1)=a;
H(1)=(sqrt((p^2+z^2)))-(sqrt(a*b));
phi(1)=atand((z/p)*((1-(N(1)*e^2)/(N(1)+H(1)))^-1));
Z=(N(1)*(b^2/a^2)+H(1))*sind(phi(1));

while abs(dN)>10^-9 && abs(dH)>10^-9 && abs(dphi)>10^-9

     N(i)=a/(sqrt((cosd(phi(i-1))^2)+((b^2/a^2)*sind(phi(i-1))^2)));
     H(i)=(p/cosd(phi(i-1)))-N(i);
     phi(i)=atand(((z/p)*((1-((N(i)*e^2)/(N(i)+H(i))))^-1)));
     
     dN=N(i)-N(i-1);
     dH=H(i)-H(i-1);
     dphi=phi(i)-phi(i-1);
     i=i+1;
end

eL=landa-long;
eP=phi(end)-lat;
eH=H(end)-h;

R=[landa phi(end) H(end)]';
E=[eL eP eH]';
%--------------------------------------------------------------------------

%Dispaly
disp(' The computed Cartesian coordinate is :')
disp("   ")
disp(r')
disp('----------------------------------------')

disp(' The computed Geodetic coordinate is :')
disp("   ");
disp(R')
disp('----------------------------------------')

disp(' The computed error is :')
disp("   ");
disp(E')
disp('----------------------------------------')

disp('Answer of Question 3 :')
disp("   ");
disp('The calculated values for geodetic longitude and latitude are equal to the initial coordinates with appropriate accuracy.')
disp('we have just order difference in the computed coordinates.')
disp('which is generally caused by the initial values with different accuracy that used in the calculations.')

%end






