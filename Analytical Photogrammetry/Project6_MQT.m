clc;
clear all;
close all;
format long g;

%project_6:MQ_transformation
%ShahabEsfandair_9819373

%InputData
[X,Y]=textread('DigitData1.txt','%f%f');
[E,N]=textread('UTMData1.txt','%f%f');
[Xc,Yc]=textread('DigitCheckData.txt','%f%f');
[Ec,Nc]=textread('UTMCheckData.txt','%f%f');
[Xp,Yp]=textread('DigitDataPlot.txt','%f%f');
[Eg,Ng]=textread('UTMdataPlot.txt','%f%f');

%LinearTerm

%MatrixCalc
j=1;
for i=1:2:length(X)*2
    A(i,:)=[X(j),Y(j),0,0,1,0];
    A(i+1,:)=[0,0,X(j),Y(j),0,1];
    L(i,1)=E(j);
    L(i+1,1)=N(j);
    j=j+1;
end

%UnknownsCalc
Xcap=inv(A'*A)*A'*L;
Lcap=A*Xcap;

%ResidueCalc
for i=1:length(X)
    Xr(i,1)=Lcap(i*2-1);
    Yr(i,1)=Lcap(i*2);
end

    dx=Xr-E;
    dy=Yr-N;
    
%GCP_DistanceCalc
 for i=1:length(X)
     for j=1:length(X)
         D(i,j)=sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
     end
 end
 
 Xa=inv(D'*D)*D'*dx;
 Xb=inv(D'*D)*D'*dy;
 
%CheckPointCalc
j=1;
for i=1:2:length(Xc)*2
    Ach(i,:)=[Xc(j),Yc(j),0,0,1,0];
    Ach(i+1,:)=[0,0,Xc(j),Yc(j),0,1];
    Lc(i,1)=Ec(j);
    Lc(i+1,1)=Nc(j);
    j=j+1;
end

Lch=Ach*Xcap;

%CheckPointDistanceCalc
for i=1:length(X)
     for j=1:length(Xc)
         Dc(i,j)=sqrt((X(i)-Xc(j))^2+(Y(i)-Yc(j))^2);
     end
end

for i=1:length(Xc)
dx_ch(i)=Dc(:,i)'*Xa;
dy_ch(i)=Dc(:,i)'*Xb;
end

%FinalCP_Coordinates
for i=1:length(Xc)
    Er(i,1)=Lch(i*2-1);
    Nr(i,1)=Lch(i*2);
    
    Xf(i)=Er(i)+dx_ch(i);
    Yf(i)=Nr(i)+dy_ch(i);
end

for i=1:3
ex(i)=Ec(i)-Xf(i);
ey(i)=Nc(i)-Yf(i);
end

et=[ex;ey];

%RMSECalc
Vx=Xf'-Ec;
Vy=Yf'-Nc;
for i=1:length(Vx)
    dr(i)=sqrt(Vx(i)^2+Vy(i)^2);
end

RMSE=sqrt(sum(dr.^2)/(length(Vx)-1));

%Disp
disp(" Coordinates of points in Image coordinate system :")
disp("    ")
disp([X,Y])
disp("-----------------------------")
disp(" Coordinates of points in Ground coordinate system :")
disp("    ")
disp([E,N])
disp("-----------------------------")
disp(" Coordinates of Checkpoints in Image coordinate system :")
disp("    ")
disp([Xc,Yc])
disp("-----------------------------")
disp(" Coordinates of Checkpoints in Ground coordinate system :")
disp("    ")
disp([Ec,Nc])
disp("-----------------------------")
disp(" Coefficients matrix A :")
disp("    ")
disp(A)
disp("-----------------------------")
disp(" Unknown Coefficients matrix Xcap :")
disp("    ")
disp(Xcap)
disp("-----------------------------")
disp(" GCP_GP calculated  :")
disp("    ")
disp(Lcap)
disp("-----------------------------")
disp(" Coefficinents matrix of distance :")
disp("    ")
disp([Xa,Xb])
disp("-----------------------------")
disp(" Final calculated coordinates :")
disp("    ")
disp([Xf;Yf])
disp("-----------------------------")
disp(" Error of calculated data and input data :")
disp("    ")
disp(et)
disp("-----------------------------")
disp(" RMSE Measure =")
disp("    ")
disp(RMSE)


%plots
subplot(2,2,1);
plot(Xp,Yp,'b->');
title('image system');

subplot(2,2,2);
plot(Eg,Ng,'b->');
title('ground system');

%end

xc_f(i,1)=(Xcap(1)*XC(i)+Xcap(2)*YC(i)+Xcap(3)*ZC(i)+Xcap(4))/(Xcap(9)*XC(i)+Xcap(10)*YC(i)+Xcap(11)*ZC(i)+1);
    yc_f(i,1)=(Xcap(5)*XC(i)+Xcap(6)*YC(i)+Xcap(7)*ZC(i)+Xcap(8))/(Xcap(9)*XC(i)+Xcap(10)*YC(i)+Xcap(11)*ZC(i)+1); 
end