clc;
clear all;
close all;
format long g;

%project_5:GP transformation_Quadric
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
for i=1:length(X)
    Ax(i,:)=[1,X(i),Y(i),X(i)*Y(i),X(i)^2,Y(i)^2];
    Lx(i,1)=E(i);
end

for i=1:length(Y)
    Ly(i,1)=N(i);
end

Ay=Ax;
%UnknownsCalc
Xcap=inv(Ax'*Ax)*Ax'*Lx;
Ycap=inv(Ax'*Ax)*Ax'*Ly;

%RMSWCalc
for i=1:length(Xc)
    Ax_ch(i,:)=[1,Xc(i),Yc(i),Xc(i)*Yc(i),Xc(i)^2,Yc(i)^2];
    Lx_c(i,1)=Ec(i);
end

Lx_ch=Ax_ch*Xcap;
dx=Lx_ch-Lx_c;

for i=1:length(Yc)
    Ly_c(i,1)=Nc(i);
end

Ay_ch=Ax_ch;
Ly_ch=Ay_ch*Ycap;
dy=Ly_ch-Ly_c;

for i=1:length(Xc)
    dr(i)=sqrt(dx(i,1)^2+dy(i,1)^2);
end

RMSE(1)=sqrt(sum(dr.^2)/(length(X)-1));

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
disp(Ax)
disp("-----------------------------")
disp(" Unknown Coefficients matrix Xcap :")
disp("    ")
disp(Xcap)
disp("-----------------------------")
disp(" Unknown Coefficients matrix Ycap :")
disp("    ")
disp(Ycap)
disp("-----------------------------")
disp(" RMSE Measure =")
disp("    ")
disp(RMSE(1))


%plots
subplot(2,2,1);
plot(Xp,Yp,'b->');
title('image system');

subplot(2,2,2);
plot(Eg,Ng,'b->');
title('ground system');

%end
