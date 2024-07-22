clc;
clear all;
close all;
format long g;

%Project1_Geometric Geodazy
%ShahabEsfandiar_9819373
%HomaGangali_9929953

%Exer2
%--------------------------------------------------------------------------

%Data
a1=6378137;
f1=1/298.2572;
x01=-25.8;
y01=-168.1;
z01=167.8;

a2=6378388;
f2=1/297;
x02=-64.6;
y02=-154.8;
z02=-46.2;

lat1=rad2deg(0.779865469);
long1=rad2deg(1.110238844);
h1=37.46;
%--------------------------------------------------------------------------

%Convert Geodetic coordinate to Cartesian
b1=a1*(1-f1);
N1=a1/(sqrt((cosd(lat1)^2)+((b1^2/a1^2)*sind(lat1)^2)));

x1=x01+(N1*cosd(lat1)*cosd(long1))+(h1*cosd(lat1)*cosd(long1));
y1=y01+(N1*cosd(lat1)*sind(long1))+(h1*cosd(lat1)*sind(long1));
z1=z01+(N1*(b1^2/a1^2)+h1)*sind(lat1);
r1=[x1 y1 z1]';
%--------------------------------------------------------------------------

%Compute Geodetic coordinate in second datum
x2=x1-x02;
y2=y1-y02;
z2=z1-z02;
landa2=atand(y2/x2);
b2=a2*(1-f2);
p2=sqrt(x2^2+y2^2);
e2=sqrt((a2^2-b2^2)/a2^2);
dN2=1;dh2=1;dphi2=1;i=2;

N2(1)=a2;
h2(1)=(sqrt((p2^2+z2^2)))-(sqrt(a2*b2));
phi2(1)=atand((z2/p2)*((1-(N2(1)*e2^2)/(N2(1)+h2(1)))^-1));
Z2=(N2(1)*(b2^2/a2^2)+h2(1))*sind(phi2(1));


while abs(dN2)>10^-9 && abs(dh2)>10^-9 && abs(dphi2)>10^-9

     N2(i)=a2/(sqrt((cosd(phi2(i-1))^2)+((b2^2/a2^2)*sind(phi2(i-1))^2)));
     h2(i)=(p2/cosd(phi2(i-1)))-N2(i);
     phi2(i)=atand(((z2/p2)*((1-((N2(i)*e2^2)/(N2(i)+h2(i))))^-1)));
     
     dN2=N2(i)-N2(i-1);
     dh2=h2(i)-h2(i-1);
     dphi2=phi2(i)-phi2(i-1);
     i=i+1;
end

eL2=landa2-long1;
eP2=phi2(end)-lat1;
eH2=h2(end)-h1;

r2=[landa2 phi2(end) h2(end)]';
E2=[eL2 eP2 eH2]';
%--------------------------------------------------------------------------

%Dispaly
disp(' The computed Cartesian coordinate for first datum is :')
disp("   ")
disp(r1')
disp('----------------------------------------')

disp(' The computed Geodetic coordinate for second datum is :')
disp("   ");
disp(r2')
disp('----------------------------------------')

disp(' The geodetic coordinates difference between two datum is :')
disp("   ");
disp(E2')

%end







