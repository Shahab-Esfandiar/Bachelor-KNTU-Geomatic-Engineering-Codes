clc;
clear all;
close all;
format long g;

%ProjectCase2_ApplicationalPhotogrametery
%ShahabEsfandiar_9819373
%HomaGangali_9929953
%--------------------------------------------------------------------------

%(0)Data
IPA=xlsread('IPA.xlsx');
GCPA=importdata('GCPA.txt');
TieA=xlsread('TieA.xlsx');
EOPA=xlsread('EOPA.xlsx');

DataPoints=IPA;
GCP=GCPA;
EOP=EOPA;

X_tie=TieA(:,2); %From Conformal
Y_tie=TieA(:,3); %From Conformal
Z_tie=TieA(:,4); %From Intersection

DataPoints(:,4)=DataPoints(:,4)/1000;
DataPoints(:,5)=DataPoints(:,5)/1000;
f=153.692*10^-3;
x0=0;y0=0;
%--------------------------------------------------------------------------

%(1)Adding GCPs coordinates to our main data
for i=1:length(DataPoints)
    
    for j=1:length(GCP)
        
        if DataPoints(i,3)==GCP(j,1)
        DataPoints(i,8)=GCP(j,2);
        DataPoints(i,9)=GCP(j,3);
        DataPoints(i,10)=GCP(j,4);
        end
    end 
end
%--------------------------------------------------------------------------

%(2)Collinearity Equation
syms X Y Z X0 Y0 Z0 omega phi kappa 

Mo=[1,           0,            0;
    0,  cos(omega),   sin(omega);
    0, -sin(omega),   cos(omega)];

Mp=[cos(phi), 0, -sin(phi);
    0       , 1,         0;
    sin(phi), 0,  cos(phi)];

Mk=[ cos(kappa),  sin(kappa), 0;
    -sin(kappa),  cos(kappa), 0;
    0          ,  0         , 1];

M=Mk*Mp*Mo;

X_col=f*((M(1,1)*(X-X0)+ M(1,2)*(Y-Y0)+ M(1,3)*(Z-Z0))/(M(3,1)*(X-X0)+ M(3,2)*(Y-Y0)+ M(3,3)*(Z-Z0))); 
Y_col=f*((M(2,1)*(X-X0)+ M(2,2)*(Y-Y0)+ M(2,3)*(Z-Z0))/(M(3,1)*(X-X0)+ M(3,2)*(Y-Y0)+ M(3,3)*(Z-Z0)));
    
dx=jacobian(X_col,[X0 Y0 Z0 omega phi kappa X Y Z]);    
dy=jacobian(Y_col,[X0 Y0 Z0 omega phi kappa X Y Z]);
%--------------------------------------------------------------------------

%(3)Compute the residuals verctor
normU=1;  
c=1;

while normU >10^-6

      for i=1:length(DataPoints)
          j=DataPoints(i,2);
          
          % Design matrix for EOP(Aeo)     
          if DataPoints(i,6)==0
             Aeo(2*i-1,6*j-5)= eval(subs(dx(1,1),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i-1,6*j-4)= eval(subs(dx(1,2),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i-1,6*j-3)= eval(subs(dx(1,3),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i-1,6*j-2)= eval(subs(dx(1,4),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i-1,6*j-1)= eval(subs(dx(1,5),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i-1,6*j-0)= eval(subs(dx(1,6),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
           
             Aeo(2*i,6*j-5)= eval(subs(dy(1,1),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i,6*j-4)= eval(subs(dy(1,2),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i,6*j-3)= eval(subs(dy(1,3),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i,6*j-2)= eval(subs(dy(1,4),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i,6*j-1)= eval(subs(dy(1,5),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Aeo(2*i,6*j-0)= eval(subs(dy(1,6),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
          
          else
             Aeo(2*i-1,6*j-5)= eval(subs(dx(1,1),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i-1,6*j-4)= eval(subs(dx(1,2),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i-1,6*j-3)= eval(subs(dx(1,3),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i-1,6*j-2)= eval(subs(dx(1,4),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i-1,6*j-1)= eval(subs(dx(1,5),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i-1,6*j-0)= eval(subs(dx(1,6),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
 
             Aeo(2*i,6*j-5)= eval(subs(dy(1,1),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i,6*j-4)= eval(subs(dy(1,2),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i,6*j-3)= eval(subs(dy(1,3),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i,6*j-2)= eval(subs(dy(1,4),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i,6*j-1)= eval(subs(dy(1,5),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Aeo(2*i,6*j-0)= eval(subs(dy(1,6),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
          end
      end
      
      % Design matrix for Tie points(Ag)
      for i=1:length(DataPoints)
          j=DataPoints(i,2);
          
          if DataPoints(i,6)==0
             Ag(2*i-1,3*DataPoints(i,7)-2)= eval(subs(dx(1,7),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Ag(2*i-1,3*DataPoints(i,7)-1)= eval(subs(dx(1,8),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Ag(2*i-1,3*DataPoints(i,7)-0)= eval(subs(dx(1,9),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));

             Ag(2*i,3*DataPoints(i,7)-2)= eval(subs(dy(1,7),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Ag(2*i,3*DataPoints(i,7)-1)= eval(subs(dy(1,8),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Ag(2*i,3*DataPoints(i,7)-0)= eval(subs(dy(1,9),[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
          end
      end
      
      A=[Aeo,Ag];  
      
      % Claculate the misclosure verctor(wp)
      for i=1:length(DataPoints)
          j=DataPoints(i,2);
          
          if DataPoints(i,6)==0
             X_comp(2*i-1,1)= eval(subs(X_col,[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             Y_comp(2*i,1)= eval(subs(Y_col,[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) X_tie(DataPoints(i,7),1) Y_tie(DataPoints(i,7),1) Z_tie(DataPoints(i,7),1)]));
             
          else
             X_comp(2*i-1,1)= eval(subs(X_col,[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
             Y_comp(2*i,1)= eval(subs(Y_col,[X0 Y0 Z0 omega phi kappa X Y Z],[EOP(j,3) EOP(j,4) EOP(j,5) EOP(j,6) EOP(j,7) EOP(j,8) DataPoints(i,8) DataPoints(i,9) DataPoints(i,10)]));
          end
      end

      for i=1:length(DataPoints)
          W(2*i-1,1)=DataPoints(i,4)+X_comp(2*i-1,1);
          W(2*i,1)=DataPoints(i,5)+Y_comp(2*i,1);
      end
      %----------------------------------------------------------------------------

      %(4) Compute unknown parameters with Least square method
      dU=-(inv(A'*A))*A'*W;
      normU=norm(dU);

      % Adding residulas vector to our approximate exterior orientation parameteres      
      for i = 1:14
          EOP(i,3)=EOP(i,3)+dU(6*i-5);
          EOP(i,4)=EOP(i,4)+dU(6*i-4);
          EOP(i,5)=EOP(i,5)+dU(6*i-3);
          EOP(i,6)=EOP(i,6)+dU(6*i-2);
          EOP(i,7)=EOP(i,7)+dU(6*i-1);
          EOP(i,8)=EOP(i,8)+dU(6*i-0);
      end
 
      % Adding residulas vector to our approximate Tie points coordinates     
      for i=1:50
          X_tie(i,1)=X_tie(i,1)+dU(84+3*i-2);
          Y_tie(i,1)=Y_tie(i,1)+dU(84+3*i-1);
          Z_tie(i,1)=Z_tie(i,1)+dU(84+3*i-0);
      end
       
      c=c+1;
end

Tie_code=unique(DataPoints(DataPoints(:,6)==0,3));
Xg_tie=X_tie;
Yg_tie=Y_tie;
Zg_tie=Z_tie;
%--------------------------------------------------------------------------

%(5)Plot
scatter3(Xg_tie,Yg_tie,Zg_tie,'* blue');
hold on;
scatter3(GCP(:,2),GCP(:,3),GCP(:,4),'* red');
legend('Tie Points','Control Points');
title('Ground coordinates');
xlabel('X');
ylabel('Y');
zlabel('Z');
%--------------------------------------------------------------------------

%(6)Table
Ran_num=DataPoints(:,1);Image_num=DataPoints(:,2);Point_code=DataPoints(:,3);
Xi=DataPoints(:,4);Yi=DataPoints(:,5);Point_type=DataPoints(:,6);Point_num=DataPoints(:,7);
XG_gcp=DataPoints(:,8);YG_gcp=DataPoints(:,9);ZG_gcp=DataPoints(:,10);
DataPoints=table(Ran_num,Image_num,Point_code,Xi,Yi,Point_type,Point_num,XG_gcp,YG_gcp,ZG_gcp);

Point_code=GCP(:,1);
Xg=GCP(:,2);Yg=GCP(:,3);Zg=GCP(:,4);
GCP=table(Point_code,Xg,Yg,Zg);

Ran_num=EOP(:,1);Image_num=EOP(:,2);
XC=EOP(:,3);YC=EOP(:,4);ZC=EOP(:,5);
Omega=EOP(:,6);Phi=EOP(:,7);Kappa=EOPA(:,8);
EOP=table(Ran_num,Image_num,XC,YC,ZC,Omega,Phi,Kappa);

Tie_Ground_coordinates=table(Tie_code,Xg_tie,Yg_tie,Zg_tie);
%--------------------------------------------------------------------------

clear dx dy GCPA i image_num IPA j kappa Kappa M Mk Mo Mp omega Omega phi Phi Point_code;
clear Point_num Point_type Ran_num Tie_code X x0 X0 X_col X_comp X_tie XC Xg XG_gcp Xi;
clear Y y0 Y0 Y_col Y_comp Y_tie YC Yg YG_gcp Yi Z Z0 ZC Zg ZG_gcp Zg_tie Xg_tie Yg_tie;
%End
