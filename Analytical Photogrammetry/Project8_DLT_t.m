clc;
clear all;
close all;
format long g;

%project_8: DLT_transformation
%ShahabEsfandair_9819373

%InputDatas
img=imread('image1.tif');
[ID c r X Y Z]=textread('AllPoints_Data.txt','%f%f%f%f%f%f');
[C_ID cC rC XC YC ZC]=textread('CPs_Data.txt','%f%f%f%f%f%f');
[G_ID cG rG XG YG ZG]=textread('GCPs_Data.txt','%f%f%f%f%f%f');

%Consts
r0=514;
c0=1894;  
Ps=0.0001;

%Convert Pixel coordinate system to photogrammetry coordinate system
for i=1:length(C_ID)
    
    xc_p(i,1)=(cC(i)-c0)*Ps;
    yc_p(i,1)=(rC(i)-r0)*Ps;
end
for i=1:length(G_ID)
    
    xg_p(i,1)=(cG(i)-c0)*Ps;
    yg_p(i,1)=(rG(i)-r0)*Ps;
end

%MatrixCalc A & L
for i=1:length(G_ID)
    A(2*i-1,:)=[XG(i) YG(i) ZG(i) 1 0 0 0 0 -xg_p(i)*XG(i) -xg_p(i)*YG(i) -xg_p(i)*ZG(i)];
    A(2*i,:)=[0 0 0 0 XG(i) YG(i) ZG(i) 1 -yg_p(i)*XG(i) -yg_p(i)*YG(i) -yg_p(i)*ZG(i)];
    L(2*i-1,1)=xg_p(i);
    L(2*i,1)=yg_p(i);
end

%UnknownsCalc
Xcap=inv(A'*A)*A'*L;

%Final CP_Calc with DLT method
for i=1:length(C_ID)
    xc_f(i,1)=(Xcap(1)*XC(i)+Xcap(2)*YC(i)+Xcap(3)*ZC(i)+Xcap(4))/(Xcap(9)*XC(i)+Xcap(10)*YC(i)+Xcap(11)*ZC(i)+1);
    yc_f(i,1)=(Xcap(5)*XC(i)+Xcap(6)*YC(i)+Xcap(7)*ZC(i)+Xcap(8))/(Xcap(9)*XC(i)+Xcap(10)*YC(i)+Xcap(11)*ZC(i)+1); 
end

%RMSECalc
for i=1:length(C_ID)
    xr(i,1)=xc_p(i)-xc_f(i);
    yr(i,1)=yc_p(i)-yc_f(i);
end

for i=1:length(C_ID)
    dr(i)=sqrt(xr(i)^2+yr(i)^2);
    theta(i)=atand(yr(i)/xr(i));
end

RMSE=sqrt(sum(dr.^2)/(length(C_ID)-1));

%Disp
disp(" Coordinates of points in Image coordinate system :")
disp("    ")
disp([r,c])
disp("-----------------------------")
disp(" Coordinates of points in Ground coordinate system :")
disp("    ")
disp([X,Y,Z])
disp("-----------------------------")
disp(" Coordinates of Checkpoints in Photogrammetry coordinate system :")
disp("    ")
disp([xc_p,yc_p])
disp("-----------------------------")
disp(" Coordinates of GCP in Photogrammetry coordinate system :")
disp("    ")
disp([xg_p,yg_p])
disp("-----------------------------")
disp(" Coefficients matrix A :")
disp("    ")
disp(A)
disp("-----------------------------")
disp(" Unknown Coefficients matrix Xcap :")
disp("    ")
disp(Xcap)
disp("-----------------------------")
disp(" Final CP calculated coordinates with DLT :")
disp("    ")
disp([xc_f;yc_f])
disp("-----------------------------")
disp(" Error of calculated data and input data :")
disp("    ")
disp([dr;theta])
disp("-----------------------------")
disp(" RMSE Measure =")
disp("    ")
disp(RMSE)

%Plots
img_s = imshow('image1.tif');
hold on
for i=1:length(ID)
    
    plot(c(i),r(i),'r.','LineWidth',20,'MarkerSize',16);
    text(c(i),r(i),num2str(ID(i)),'color','white','FontSize',14);
end

%end



