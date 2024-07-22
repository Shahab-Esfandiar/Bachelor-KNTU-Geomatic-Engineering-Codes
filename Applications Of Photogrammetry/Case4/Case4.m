clc;
clear all;
close all;
format long g;

%ProjectCase4_ApplicationalPhotogrametery
%ShahabEsfandiar_9819373
%HomaGangali_9929953
%--------------------------------------------------------------------------

%(0)Data
DataPoints=xlsread('data.xlsx');
TieA=xlsread('tie.xlsx');
EOPA=xlsread('EOP.xlsx');
Full=importdata('full.txt');
Full_comp=xlsread('full_com.xlsx');

DataPoints(:,4)=DataPoints(:,4)/1000;
DataPoints(:,5)=DataPoints(:,5)/1000;
EOP=EOPA;
Tie=TieA;
f=0.153692;
%--------------------------------------------------------------------------

%(1)Collinearity Equation
syms X Y Z X0 Y0 Z0 omega phi kapa x y

M=[cos(kapa) sin(kapa) 0;-sin(kapa) cos(kapa) 0;0 0 1]*[cos(phi) 0 -sin(phi);0 1 0;sin(phi) 0 cos(phi)]*[1 0 0;0 cos(omega) sin(omega);0 -sin(omega) cos(omega)];

    Fx=x+f*((M(1,1)*(X-X0) + M(1,2)*(Y-Y0) + M(1,3)*(Z-Z0))/(M(3,1)*(X-X0) + M(3,2)*(Y-Y0) + M(3,3)*(Z-Z0))); 
    Fy=y+f*((M(2,1)*(X-X0) + M(2,2)*(Y-Y0) + M(2,3)*(Z-Z0))/(M(3,1)*(X-X0) + M(3,2)*(Y-Y0) + M(3,3)*(Z-Z0)));
   
FAeo=[jacobian(Fx,[X0 Y0 Z0 omega phi kapa]);
      jacobian(Fy,[X0 Y0 Z0 omega phi kapa])];
   
FAg=[jacobian(Fx,[X Y Z]);
     jacobian(Fy,[X Y Z])];
%--------------------------------------------------------------------------

%(2)Compute the residuals and weight verctor  
A=zeros(424,246);
Acc(1:328)=(1/(7*(10^-6))^2);
P=diag(Acc);
normU=2;
c=1;

while normU >10^-10
    
    for i=1:length(DataPoints)
        img_num=DataPoints(i,2);

        X0=EOPA(img_num,3);
        Y0=EOPA(img_num,4);
        Z0=EOPA(img_num,5);
        omega=EOPA(img_num,6);
        phi=EOPA(img_num,7);
        kapa=EOPA(img_num,8);
        
        % Design matrix A
        if DataPoints(i,6)==0
           X=TieA(DataPoints(i,7),2);
           Y=TieA(DataPoints(i,7),3);
           Z=TieA(DataPoints(i,7),4);
    
           A(2*i-1:2*i,84+3*DataPoints(i,7)-2:84+3*DataPoints(i,7))=[(2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2);
                                                                    -(2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2)];

        else
           X=Full_comp(DataPoints(i,8),3);
           Y=Full_comp(DataPoints(i,8),4);
           Z=Full_comp(DataPoints(i,8),5);      

           A(2*i-1:2*i,234+3*DataPoints(i,8)-2:234+3*DataPoints(i,8))=[(2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2);
                                                                      -(2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2)];
        end
        
    A(2*i-1:2*i,6*img_num-5:6*img_num)=[(2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(omega)*cos(phi)*(Y - Y0) + cos(phi)*sin(omega)*(Z - Z0))*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*((sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Y - Y0) - (cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Z - Z0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(kapa)*sin(phi)*(X - X0) + cos(kapa)*cos(omega)*cos(phi)*(Z - Z0) - cos(kapa)*cos(phi)*sin(omega)*(Y - Y0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*(cos(phi)*(X - X0) - cos(omega)*sin(phi)*(Z - Z0) + sin(omega)*sin(phi)*(Y - Y0))*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2),  (2768668935719305*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
                                        (2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) + (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(omega)*cos(phi)*(Y - Y0) + cos(phi)*sin(omega)*(Z - Z0))*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*((cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Y - Y0) - (cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Z - Z0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))),   (2768668935719305*(sin(kapa)*sin(phi)*(X - X0) + cos(omega)*cos(phi)*sin(kapa)*(Z - Z0) - cos(phi)*sin(kapa)*sin(omega)*(Y - Y0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*(cos(phi)*(X - X0) - cos(omega)*sin(phi)*(Z - Z0) + sin(omega)*sin(phi)*(Y - Y0))*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), -(2768668935719305*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)))];
                                      
    % Claculate the misclosure verctor(w)                                  
    W(2*i-1,1)=DataPoints(i,4)+(2768668935719305*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
    W(2*i-0,1)=DataPoints(i,5)+(2768668935719305*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
    end
   
    % Claculate the Covariance matrix(P)
    for i=1:4
        A(328+3*i-2,234+3*i-2)=-1;
        A(328+3*i-1,234+3*i-1)=-1;
        A(328+3*i-0,234+3*i-0)=-1;

        P(328+3*i-2,328+3*i-2)=(1/(15*(10^-2))^2);
        P(328+3*i-1,328+3*i-1)=(1/(15*(10^-2))^2);
        P(328+3*i-0,328+3*i-0)=(1/(20*(10^-2))^2);

        W(328+3*i-2,1)=Full(i,2)-Full_comp(i,3);
        W(328+3*i-1,1)=Full(i,3)-Full_comp(i,4);
        W(328+3*i,1-0)=Full(i,4)-Full_comp(i,5);
    end
     
    % For weighted EOP
    A(341:424,1:84)=-eye(84);
    
    for i=1:14
        P(340+6*i-5,340+6*i-5)=(1/(20*(10^-2))^2);
        P(340+6*i-4,340+6*i-4)=(1/(20*(10^-2))^2);
        P(340+6*i-3,340+6*i-3)=(1/(25*(10^-2))^2);
        P(340+6*i-2,340+6*i-2)=(1/((5*pi)/(3600*180)))^2;
        P(340+6*i-1,340+6*i-1)=(1/((5*pi)/(3600*180)))^2;
        P(340+6*i-0,340+6*i-0)=(1/((5*pi)/(3600*180)))^2;

        W(340+6*i-5,1)=EOP(i,3)-EOPA(i,3);
        W(340+6*i-4,1)=EOP(i,4)-EOPA(i,4);
        W(340+6*i-3,1)=EOP(i,5)-EOPA(i,5);
        W(340+6*i-2,1)=EOP(i,6)-EOPA(i,6);
        W(340+6*i-1,1)=EOP(i,7)-EOPA(i,7);
        W(340+6*i,1-0)=EOP(i,8)-EOPA(i,8);
    end
    
    %(3) Compute unknown parameters with Least square method
    dU=-inv(A'*P*A)*A'*P*W;
    normU=norm(dU);
    
    % Adding residulas vector to our approximate exterior orientation parameteres
    for i=1:14
        EOPA(i,3)=EOPA(i,3)+dU(6*i-5);
        EOPA(i,4)=EOPA(i,4)+dU(6*i-4);
        EOPA(i,5)=EOPA(i,5)+dU(6*i-3);
        EOPA(i,6)=EOPA(i,6)+dU(6*i-2);
        EOPA(i,7)=EOPA(i,7)+dU(6*i-1);
        EOPA(i,8)=EOPA(i,8)+dU(6*i-0);
    end
    
    % Adding residulas vector to our approximate Tie points coordinates
    for i=1:50
        TieA(i,2)=TieA(i,2)+dU(84+3*i-2);
        TieA(i,3)=TieA(i,3)+dU(84+3*i-1);
        TieA(i,4)=TieA(i,4)+dU(84+3*i-0);
    end
    
    % Adding residulas vector to our weighted Control points coordinates
    for i=1:4
        Full_comp(i,3)=Full_comp(i,3)+dU(234+3*i-2);
        Full_comp(i,4)=Full_comp(i,4)+dU(234+3*i-1);
        Full_comp(i,5)=Full_comp(i,5)+dU(234+3*i-0);       
    end
    
    c=c+1;
end

Tie_code=unique(DataPoints(DataPoints(:,6)==0,3));
Xg_tie=TieA(:,2);
Yg_tie=TieA(:,3);
Zg_tie=TieA(:,4);
%--------------------------------------------------------------------------

%(4)Plot
scatter3(Xg_tie,Yg_tie,Zg_tie,'* blue');
hold on;
scatter3(Full_comp(:,3),Full_comp(:,4),Full_comp(:,5),'* red');
legend('Tie Points','Full Control Points');
title('Ground coordinates');
xlabel('X');
ylabel('Y');
zlabel('Z');
%--------------------------------------------------------------------------

%(6)Table
Ran_num=DataPoints(:,1);Image_num=DataPoints(:,2);Point_code=DataPoints(:,3);
Xi=DataPoints(:,4);Yi=DataPoints(:,5);Point_type=DataPoints(:,6);Point_num=DataPoints(:,7);
DataPoint=table(Ran_num,Image_num,Point_code,Xi,Yi,Point_type,Point_num);

Point_code=Full_comp(:,2);
Xg=Full_comp(:,3);Yg=Full_comp(:,4);Zg=Full_comp(:,5);
Full_comp=table(Point_code,Xg,Yg,Zg);

Ran_num=EOP(:,1);Image_num=EOP(:,2);
XC=EOP(:,3);YC=EOP(:,4);ZC=EOP(:,5);
Omega=EOP(:,6);Phi=EOP(:,7);Kappa=EOP(:,8);
EOP=table(Ran_num,Image_num,XC,YC,ZC,Omega,Phi,Kappa);

Tie_Ground_coordinates=table(Tie_code,Xg_tie,Yg_tie,Zg_tie);
%--------------------------------------------------------------------------

clear dx dy i image_num IPA j kappa Kappa M Mk Mo Mp omega Omega phi Phi Point_code;
clear Point_num Point_type Ran_num Tie_code X x0 X0 X_col X_comp X_tie XC Xg XG_gcp Xi;
clear Y y0 Y0 Y_col Y_comp Y_tie YC Yg YG_gcp Yi Z Z0 ZC Zg ZG_gcp Zg_tie Yw Zw TieA;
clear Acc  Con_comp DataPoints FAeo FAg Fx Fy Image_num img_num kapa x Xg_tie Xw y Yg_tie;
%End
