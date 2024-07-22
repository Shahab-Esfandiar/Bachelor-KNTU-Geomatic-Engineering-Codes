clc;
clear all;
close all;
format long g;

%ProjectCase3_ApplicationalPhotogrametery
%ShahabEsfandiar_9819373
%HomaGangali_9929953
%--------------------------------------------------------------------------

%(0)Data
DataPoints=xlsread('data_main.xlsx');
TieA=xlsread('tie.xlsx');
GCP=importdata('Control.txt');
EOP=xlsread('EOP.xlsx');
Con_comp=xlsread('con_com.xlsx');

DataPoints(:,4)=DataPoints(:,4)/1000;
DataPoints(:,5)=DataPoints(:,5)/1000;
DataPoints(:,12)=0;
Tie=TieA;
f=0.153692;
%--------------------------------------------------------------------------

%(1)Adding GCPs number to our main data
for i=1:length(DataPoints)
    
    for j=1:length(GCP)
        
        if DataPoints(i,3)== GCP(j,1)
           DataPoints(i,12)=j;
        end
    end
end
%--------------------------------------------------------------------------

%(2)Collinearity Equation
syms X Y Z X0 Y0 Z0 omega phi kapa x y

M=[cos(kapa) sin(kapa) 0;-sin(kapa) cos(kapa) 0;0 0 1]*[cos(phi) 0 -sin(phi);0 1 0;sin(phi) 0 cos(phi)]*[1 0 0;0 cos(omega) sin(omega);0 -sin(omega) cos(omega)];

Fx=x+f*((M(1,1)*(X-X0) + M(1,2)*(Y-Y0) + M(1,3)*(Z-Z0))/(M(3,1)*(X-X0) + M(3,2)*(Y-Y0) + M(3,3)*(Z-Z0))); 
Fy=y+f*((M(2,1)*(X-X0) + M(2,2)*(Y-Y0) + M(2,3)*(Z-Z0))/(M(3,1)*(X-X0) + M(3,2)*(Y-Y0) + M(3,3)*(Z-Z0)));
   
FAeo=[jacobian(Fx,[X0 Y0 Z0 omega phi kapa]);
      jacobian(Fy,[X0 Y0 Z0 omega phi kapa])];
   
FAg=[jacobian(Fx,[X Y Z]);
     jacobian(Fy,[X Y Z])];
%--------------------------------------------------------------------------

%(3)Compute the residuals and weight verctor
normU=1;
A=zeros(405,267);
Acc(1:376,1)=(1/(7*(10^-6))^2);
P=diag(Acc);
c=1;

while normU >10^-10
    
    for i=1:length(DataPoints)
        img_num=DataPoints(i,2);

        X0=EOP(img_num,3);
        Y0=EOP(img_num,4);
        Z0=EOP(img_num,5);
        omega=EOP(img_num,6);
        phi=EOP(img_num,7);
        kapa=EOP(img_num,8);
      
      % Design matrix A 
        if DataPoints(i,6)==0
           X=TieA(DataPoints(i,7),2);
           Y=TieA(DataPoints(i,7),3);
           Z=TieA(DataPoints(i,7),4);
    
           A(2*i-1:2*i,84+3*DataPoints(i,7)-2:84+3*DataPoints(i,7))=[(2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2);
                                                                    -(2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2)];
        else
           X=Con_comp(DataPoints(i,12),3);
           Y=Con_comp(DataPoints(i,12),4);
           Z=Con_comp(DataPoints(i,12),5);      

           A(2*i-1:2*i,234+3*DataPoints(i,12)-2:234+3*DataPoints(i,12))=[(2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2);
                                                                        -(2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) + (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2)];                                            
        end
      
    A(2*i-1:2*i,6*img_num-5:6*img_num)=[(2768668935719305*sin(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*cos(kapa)*cos(phi))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(phi)*sin(omega)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*cos(omega)*cos(phi)*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*(sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(omega)*cos(phi)*(Y - Y0) + cos(phi)*sin(omega)*(Z - Z0))*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*((sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Y - Y0) - (cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Z - Z0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(kapa)*sin(phi)*(X - X0) + cos(kapa)*cos(omega)*cos(phi)*(Z - Z0) - cos(kapa)*cos(phi)*sin(omega)*(Y - Y0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*(cos(phi)*(X - X0) - cos(omega)*sin(phi)*(Z - Z0) + sin(omega)*sin(phi)*(Y - Y0))*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2),  (2768668935719305*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
                                        (2768668935719305*sin(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) + (2768668935719305*cos(phi)*sin(kapa))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), - (2768668935719305*(cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*cos(phi)*sin(omega)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), (2768668935719305*cos(omega)*cos(phi)*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*(cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))), (2768668935719305*(cos(omega)*cos(phi)*(Y - Y0) + cos(phi)*sin(omega)*(Z - Z0))*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2) - (2768668935719305*((cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Y - Y0) - (cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Z - Z0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))),   (2768668935719305*(sin(kapa)*sin(phi)*(X - X0) + cos(omega)*cos(phi)*sin(kapa)*(Z - Z0) - cos(phi)*sin(kapa)*sin(omega)*(Y - Y0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))) - (2768668935719305*(cos(phi)*(X - X0) - cos(omega)*sin(phi)*(Z - Z0) + sin(omega)*sin(phi)*(Y - Y0))*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0))^2), -(2768668935719305*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)))];
                                      
    % Claculate the misclosure verctor(w)                                  
    W(2*i-1,1)=DataPoints(i,4)+(2768668935719305*((cos(omega)*sin(kapa) + cos(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (sin(kapa)*sin(omega) - cos(kapa)*cos(omega)*sin(phi))*(Z - Z0) + cos(kapa)*cos(phi)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
    W(2*i,1)=DataPoints(i,5)+(2768668935719305*((cos(kapa)*cos(omega) - sin(kapa)*sin(omega)*sin(phi))*(Y - Y0) + (cos(kapa)*sin(omega) + cos(omega)*sin(kapa)*sin(phi))*(Z - Z0) - cos(phi)*sin(kapa)*(X - X0)))/(18014398509481984*(sin(phi)*(X - X0) + cos(omega)*cos(phi)*(Z - Z0) - cos(phi)*sin(omega)*(Y - Y0)));
    end
    
    for i=1:11
    
    % Claculate the Covariance matrix(P) 
        if i==1 || i==2 || i==3
           A(376+3*i-2,234+3*i-2)=-1;
           A(376+3*i-1,234+3*i-1)=-1;
           A(376+3*i,234+3*i)=-1;

           P(376+3*i-2,376+3*i-2)=(1/(15*(10^-2))^2);
           P(376+3*i-1,376+3*i-1)=(1/(15*(10^-2))^2);
           P(376+3*i,376+3*i)=(1/(20*(10^-2))^2);

           W(376+3*i-2,1)=GCP(i,2)-Con_comp(i,3);
           W(376+3*i-1,1)=GCP(i,3)-Con_comp(i,4);
           W(376+3*i,1)=GCP(i,4)-Con_comp(i,5);
        

        elseif i==4
               A(376+3*i-4,234+3*i-2)=-1;
               A(376+3*i-3,234+3*i-1)=-1;
        
               P(376+3*i-4,376+3*i-4)=(1/(15*(10^-2))^2);
               P(376+3*i-3,376+3*i-3)=(1/(15*(10^-2))^2);
             
               W(376+3*i-4,1)=GCP(i,2)-Con_comp(i,3);
               W(376+3*i-3,1)=GCP(i,3)-Con_comp(i,4);
            
        elseif i==5 
               A(376+3*i-2,234+3*i)=-1;
               P(376+3*i-2,376+3*i-2)=(1/(20*(10^-2))^2);
               W(376+3*i-2,1)=GCP(i,4)-Con_comp(i,5);
        
        
        elseif i==6 || i==7 || i==8 || i==9
               A(376+3*i-5,234+3*i-2)=-1;
               A(376+3*i-4,234+3*i-1)=-1;
               A(376+3*i-3,234+3*i)=-1;

               P(376+3*i-5,376+3*i-5)=(1/(15*(10^-2))^2);
               P(376+3*i-4,376+3*i-4)=(1/(15*(10^-2))^2);
               P(376+3*i-3,376+3*i-3)=(1/(20*(10^-2))^2);

               W(376+3*i-5,1)=GCP(i,2)-Con_comp(i,3);
               W(376+3*i-4,1)=GCP(i,3)-Con_comp(i,4);
               W(376+3*i-3,1)=GCP(i,4)-Con_comp(i,5);
           
        elseif i==10 
               A(376+3*i-5,234+3*i-2)=-1;
               A(376+3*i-4,234+3*i-1)=-1;
        
               P(376+3*i-5,376+3*i-5)=(1/(15*(10^-2))^2);
               P(376+3*i-4,376+3*i-4)=(1/(15*(10^-2))^2);
             
               W(376+3*i-5,1)=GCP(i,2)-Con_comp(i,3);
               W(376+3*i-4,1)=GCP(i,3)-Con_comp(i,4);
           
        else
            A(376+3*i-6,234+3*i-2)=-1;
            A(376+3*i-5,234+3*i-1)=-1;
            A(376+3*i-4,234+3*i)=-1;

            P(376+3*i-6,376+3*i-6)=(1/(15*(10^-2))^2);
            P(376+3*i-5,376+3*i-5)=(1/(15*(10^-2))^2);
            P(376+3*i-4,376+3*i-4)=(1/(20*(10^-2))^2);

            W(376+3*i-6,1)=GCP(i,2)-Con_comp(i,3);
            W(376+3*i-5,1)=GCP(i,3)-Con_comp(i,4);
            W(376+3*i-4,1)=GCP(i,4)-Con_comp(i,5);
        end
    end
%--------------------------------------------------------------------------
    
    %(4) Compute unknown parameters with Least square method
    dU=-inv(A'*P*A)*A'*P*W;
    normU=norm(dU);
    
    % Adding residulas vector to our approximate exterior orientation parameteres
    for i=1:14
        EOP(i,3)=EOP(i,3)+dU(6*i-5);
        EOP(i,4)=EOP(i,4)+dU(6*i-4);
        EOP(i,5)=EOP(i,5)+dU(6*i-3);
        EOP(i,6)=EOP(i,6)+dU(6*i-2);
        EOP(i,7)=EOP(i,7)+dU(6*i-1);
        EOP(i,8)=EOP(i,8)+dU(6*i);
    end
    
    % Adding residulas vector to our approximate Tie points coordinates
    for i=1:50
        TieA(i,2)=TieA(i,2)+dU(84+3*i-2);
        TieA(i,3)=TieA(i,3)+dU(84+3*i-1);
        TieA(i,4)=TieA(i,4)+dU(84+3*i);
    end
   
    % Adding residulas vector to our weighted Control points coordinates
    for i=1:11
        Con_comp(i,3)=Con_comp(i,3)+dU(234+3*i-2);
        Con_comp(i,4)=Con_comp(i,4)+dU(234+3*i-1);
        Con_comp(i,5)=Con_comp(i,5)+dU(234+3*i);
    
    end
     
    c=c+1;
end

Tie_code=unique(DataPoints(DataPoints(:,6)==0,3));
Xg_tie=TieA(:,2);
Yg_tie=TieA(:,3);
Zg_tie=TieA(:,4);
%--------------------------------------------------------------------------

%(5)Plot
scatter3(Xg_tie,Yg_tie,Zg_tie,'* blue');
hold on;
scatter3(Con_comp(:,3),Con_comp(:,4),Con_comp(:,5),'* red');
legend('Tie Points','Weighted Control Points');
title('Ground coordinates');
xlabel('X');
ylabel('Y');
zlabel('Z');
%--------------------------------------------------------------------------

%(6)Table
Ran_num=DataPoints(:,1);Image_num=DataPoints(:,2);Point_code=DataPoints(:,3);
Xi=DataPoints(:,4);Yi=DataPoints(:,5);Point_type=DataPoints(:,6);Point_num=DataPoints(:,7);
XG_gcp=DataPoints(:,8);YG_gcp=DataPoints(:,9);ZG_gcp=DataPoints(:,10);
DataPoint=table(Ran_num,Image_num,Point_code,Xi,Yi,Point_type,Point_num,XG_gcp,YG_gcp,ZG_gcp);

Point_code=GCP(:,1);
Xg=GCP(:,2);Yg=GCP(:,3);Zg=GCP(:,4);
GCP=table(Point_code,Xg,Yg,Zg);

Ran_num=EOP(:,1);Image_num=EOP(:,2);
XC=EOP(:,3);YC=EOP(:,4);ZC=EOP(:,5);
Omega=EOP(:,6);Phi=EOP(:,7);Kappa=EOP(:,8);
EOP=table(Ran_num,Image_num,XC,YC,ZC,Omega,Phi,Kappa);

Point_code=Con_comp(:,2);
Xw=Con_comp(:,3);Yw=Con_comp(:,4);Zw=Con_comp(:,5);
Cont_comp=table(Point_code,Xw,Yw,Zw);

Tie_Ground_coordinates=table(Tie_code,Xg_tie,Yg_tie,Zg_tie);
%--------------------------------------------------------------------------

clear dx dy i image_num IPA j kappa Kappa M Mk Mo Mp omega Omega phi Phi Point_code;
clear Point_num Point_type Ran_num Tie_code X x0 X0 X_col X_comp X_tie XC Xg XG_gcp Xi;
clear Y y0 Y0 Y_col Y_comp Y_tie YC Yg YG_gcp Yi Z Z0 ZC Zg ZG_gcp Zg_tie Yw Zw TieA;
clear Acc  Con_comp DataPoints FAeo FAg Fx Fy Image_num img_num kapa x Xg_tie Xw y Yg_tie;
%End

