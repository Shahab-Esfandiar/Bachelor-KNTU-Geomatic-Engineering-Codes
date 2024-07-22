clc;
clear all;
close all;

%project_5:GP transformation_Linear
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

%RMSWCalc
j=1;
for i=1:2:length(Xc)*2
    Ach(i,:)=[Xc(j),Yc(j),0,0,1,0];
    Ach(i+1,:)=[0,0,Xc(j),Yc(j),0,1];
    Lc(i,1)=Ec(j);
    Lc(i+1,1)=Nc(j);
    j=j+1;
end

Lch=Ac*Xcap;
d=Lch-Lc;

j=1;
for i=1:2:length(Xc)*2
    dr(i)=sqrt(d(i,1)^2+d(i+1,1)^2)
    j=j+1;
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
disp(" Coefficients matrix :")
disp("    ")
disp(A)
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
