clc;
clear all;
close all;
format long g;

%Project1_CaseApproximate_ApplicationalPhotogrametery
%ShahabEsfandiar_9819373
%HomaGangali_9929953

%--------------------------------------------------------------------------

%(0)Data
Image_Points=xlsread('D:\UNI\Term 7\Karbord Photogrametery\Amaliat\Project1\1\data.xlsx');
Ground_Control_Points=xlsread('D:\UNI\Term 7\Karbord Photogrametery\Amaliat\Project1\1\data_control.xlsx');
EOPA = xlsread('D:\UNI\Term 7\Karbord Photogrametery\Amaliat\Project1\1\EOPA.xlsx');

IP=Image_Points;
GCP=Ground_Control_Points;
f=153.692*10^-3;
x0=0;
y0=0;
%--------------------------------------------------------------------------

%(1)Specify the type of points : TiePoint=0,FullPoint=1,FlatPoint=2,ElevationPoint=3
for i=1:length(IP)
    
    for j=1:length(GCP)
        
        if IP(i,3)==GCP(j,1) & GCP(j,2:4)~=0 
           IP(i,6)=1;
           
        else if IP(i,3)==GCP(j,1) & GCP(j,4)==0
                IP(i,6)=2;    
                
        else if IP(i,3)==GCP(j,1) & GCP(j,2:3)==0
                IP(i,6)=3;  
        end
        end
        end
    end
end

Tie_code=unique(IP(IP(:,6)==0,3));
%--------------------------------------------------------------------------

%(2)Numbering Points
IP(:,7)=0;
control_id_temp=unique(IP(:,3));
len_temp=length(control_id_temp);
T=1;C=101;H=1001;V=10001;

check1_temp=false;
check2_temp=false;
check3_temp=false;
check4_temp=false;

for i=1:len_temp
    
    if check1_temp==true
       T=T+1;
    end
    
    if check2_temp==true
       C=C+1;
    end
    
    if check3_temp==true
       H=H+1;
    end
    
    if check4_temp==true
       V=V+1;
    end
    
    check1_temp=false;
    check2_temp=false;
    check3_temp=false;
    check4_temp=false;
    
    for j=1:length(IP)
        
        if IP(j,6)==0 && control_id_temp(i)==IP(j,3)
           IP(j,7)=T;
           check1_temp=true;
           
        elseif IP(j,6)==1 && control_id_temp(i)==IP(j,3)
               IP(j, 7)=C;
               check2_temp=true;
            
        elseif IP(j, 6)==2 && control_id_temp(i)==IP(j,3)
               IP(j, 7)=H;
               check3_temp=true;
            
        elseif IP(j, 6)==3 && control_id_temp(i)==IP(j,3)
               IP(j,7)=V;
               check4_temp=true;
        end
    end
end

Points=unique(IP(:,3));
nTie=T-1;
nFull=C-100;
nHor=H-1001;
nVer=V-10001;
%--------------------------------------------------------------------------

%(3)Create binary matrix that shows observed points in each image
for i=1:length(IP)
    
    if IP(i,1)==2
        IP(i,2)=IP(i,2)+7;
    end
    
end

nImage=length(unique(IP(:,2)));
biT=zeros(nImage,nTie);
biF=zeros(nImage,nFull);

%For Tie points
for i=1:nImage
    [r,c]=find(IP(:,2)==i);
    
    for j=1:length(r)
        
        for k=1:nTie
            
            if IP(r(j),6)==0 && IP(r(j),7)==k
                biT(i,k)=1;
            end
        end
    end
end

%For Full control points
for i=1:nImage
    [r,c]=find(IP(:,2)==i);
    
    for j=1:length(r)
        
        for k=1:nFull
            
            if IP(r(j),6)==1 && IP(r(j),7)==k+100
                biF(i,k)=1;
            end
        end
    end
end
%--------------------------------------------------------------------------

%(4)Create coefficient and observation matrix and compute initial values
IP(:,4)=IP(:,4)*10^-3;
IP(:,5)=IP(:,5)*10^-3;
L=zeros(2*length(IP),1);

for i=1:length(IP)
    
    j=IP(i,2);
    
    Au(2*i-1,4*j-3)=IP(i,4);
    Au(2*i-1,4*j-2)=IP(i,5);
    Au(2*i-1,4*j-1)=1;
    Au(2*i-1,4*j)=0; 
     
    Au(2*i,4*j-3)=IP(i,5);
    Au(2*i,4*j-2)=-IP(i,4);
    Au(2*i,4*j-1)=0; 
    Au(2*i,4*j)=1;

    if IP(i,6)==0
       k=IP(i,7);
       Ag(2*i-1:2*i,2*k-1:2*k)=-eye(2,2);
       
       
    else if IP(i,6)==3
            Ag(2*i-1:2*i,101:102)=-eye(2,2);
            
            
    else if IP(i,6)==1 || IP(i,6)==2
            
            for j=1:length(GCP)
                
                if IP(i,3)==GCP(j,1)
                   L(2*i-1:2*i,1)=[GCP(j,2);GCP(j,3)];
                end
            end
    end
    end
    end
end

A=[Au,Ag];
x_cap=inv(A'*A)*A'*L; 
%--------------------------------------------------------------------------

%(5)Calculate external orientation parameteres
Z_avg=sum(GCP(:,4))/9;
omega=0;
phi=0;
kappa=0;

for i=1:nImage
    
    cRan=1;
    image_num=i;
    a=x_cap(4*i-3);
    b=x_cap(4*i-2);
    c=x_cap(4*i-1);
    d=x_cap(4*i);
    landa=sqrt(a^2+b^2);
    X0=c;
    Y0=d;
    Z0=landa*f+Z_avg;
    
    if i > 7
       cRan=2;
    end
    
    if b>0 && a>0
       kappa=atan(b/a);
       
    else if b>0 && a<0
         kappa=pi-atan(abs(b/a));
        
    else if b<0 && a<0
         kappa=pi+atan(abs(b/a));
        
    else if b<0 && a>0
         kappa=2*pi-atan(abs(b/a));
    end
    end
    end
    end
    
    EOP(i,:)=[cRan image_num X0 Y0 Z0 omega phi kappa];
end
%--------------------------------------------------------------------------

%(6)Space intersection

%For tie Points
for i=1:nTie
    [r,c]=find(biT(:,i)==1);
    [r1,c1]=find(IP(:,7)==i & IP(:,2)==r(1));
    [r2,c2]=find(IP(:,7)==i & IP(:,2)==r(2));
    
    kappaL=EOP(r(1),8);
    kappaR=EOP(r(2),8);
    
    ML=[cos(kappaL), -sin(kappaL), 0;
        sin(kappaL), cos(kappaL) , 0;
        0          , 0           , 1];
    
    MR=[cos(kappaR), -sin(kappaR), 0;
        sin(kappaR), cos(kappaR) , 0;
        0          , 0           , 1];
    
    UL=(ML(1,1)*IP(r1,4))+(ML(2,1)*IP(r1,5))-(ML(3,1)*f);
    VL=(ML(1,2)*IP(r1,4))+(ML(2,2)*IP(r1,5))-(ML(3,2)*f);
    WL=(ML(1,3)*IP(r1,4))+(ML(2,3)*IP(r1,5))-(ML(3,3)*f);
    
    UR=(MR(1,1)*IP(r2,4))+(MR(2,1)*IP(r2,5))-(MR(3,1)*f);
    VR=(MR(1,2)*IP(r2,4))+(MR(2,2)*IP(r2,5))-(MR(3,2)*f);
    WR=(MR(1,3)*IP(r2,4))+(MR(2,3)*IP(r2,5))-(MR(3,3)*f);
    
    K=(((EOP(r(2),3)-(EOP(r(1),3)))*VL)-((EOP(r(2),4)-(EOP(r(1),4)))*UL))/(VR*UL-VL*UR);
    
    XL_tie(i)=(K*UL)+(EOP(r(1),3));
    YL_tie(i)=(K*VL)+(EOP(r(1),4));
    ZL_tie(i)=(K*WL)+(EOP(r(1),5));
    
    XR_tie(i)=(K*UR)+(EOP(r(2),3));
    YR_tie(i)=(K*VR)+(EOP(r(2),4));
    ZR_tie(i)=(K*WR)+(EOP(r(2),5));
end

Xt=((XL_tie+XR_tie)/2)';
Yt=((YL_tie+YR_tie)/2)';
Zt=((ZL_tie+ZR_tie)/2)';

%For Full control Points
for i=1:nFull
    [r,c]=find(biF(:,i)==1);
    [r1,c1]=find(IP(:,7)==i+100 & IP(:,2)==r(1));
    [r2,c2]=find(IP(:,7)==i+100 & IP(:,2)==r(2));
   
    kappaL=EOP(r(1),8);
    kappaR=EOP(r(2),8);
    
    ML=[cos(kappaL), -sin(kappaL), 0;
        sin(kappaL), cos(kappaL) , 0;
        0          , 0           , 1];
    
    MR=[cos(kappaR), -sin(kappaR), 0;
        sin(kappaR), cos(kappaR) , 0;
        0          , 0           , 1];
    
    UL=(ML(1,1)*IP(r1,4))+(ML(2,1)*IP(r1,5))-(ML(3,1)*f);
    VL=(ML(1,2)*IP(r1,4))+(ML(2,2)*IP(r1,5))-(ML(3,2)*f);
    WL=(ML(1,3)*IP(r1,4))+(ML(2,3)*IP(r1,5))-(ML(3,3)*f);
    
    UR=(MR(1,1)*IP(r2,4))+(MR(2,1)*IP(r2,5))-(MR(3,1)*f);
    VR=(MR(1,2)*IP(r2,4))+(MR(2,2)*IP(r2,5))-(MR(3,2)*f);
    WR=(MR(1,3)*IP(r2,4))+(MR(2,3)*IP(r2,5))-(MR(3,3)*f);
    
    K=(((EOP(r(2),3)-(EOP(r(1),3)))*VL)-((EOP(r(2),4)-(EOP(r(1),4)))*UL))/(VR*UL-VL*UR);
    
    XL_Control(i)=(K*UL)+(EOP(r(1),3));
    YL_Control(i)=(K*VL)+(EOP(r(1),4));
    ZL_Control(i)=(K*WL)+(EOP(r(1),5));
    
    XR_Control(i)=(K*UR)+(EOP(r(2),3));
    YR_Control(i)=(K*VR)+(EOP(r(2),4));
    ZR_Control(i)=(K*WR)+(EOP(r(2),5));
end

Xc=((XL_Control+XR_Control)/2)';
Yc=((YL_Control+YR_Control)/2)';
Zc=((ZL_Control+ZR_Control)/2)';

%For Flat control Points
r=zeros(2,1);

for i=1:nHor
    
    if i+1000==1001
       r(1)=4;
       r(2)=5;
       r1=46;
       r2=58;
    end
       
    if i+1000==1002
       r(1)=11;
       r(2)=12;
       r1=141;
       r2=154;
    end
   
    kappaL=EOP(r(1),8);
    kappaR=EOP(r(2),8);
    
    ML=[cos(kappaL), -sin(kappaL), 0;
        sin(kappaL), cos(kappaL) , 0;
        0          , 0           , 1];
    
    MR=[cos(kappaR), -sin(kappaR), 0;
        sin(kappaR), cos(kappaR) , 0;
        0          , 0           , 1];
    
    UL=(ML(1,1)*IP(r1,4))+(ML(2,1)*IP(r1,5))-(ML(3,1)*f);
    VL=(ML(1,2)*IP(r1,4))+(ML(2,2)*IP(r1,5))-(ML(3,2)*f);
    WL=(ML(1,3)*IP(r1,4))+(ML(2,3)*IP(r1,5))-(ML(3,3)*f);
    
    UR=(MR(1,1)*IP(r2,4))+(MR(2,1)*IP(r2,5))-(MR(3,1)*f);
    VR=(MR(1,2)*IP(r2,4))+(MR(2,2)*IP(r2,5))-(MR(3,2)*f);
    WR=(MR(1,3)*IP(r2,4))+(MR(2,3)*IP(r2,5))-(MR(3,3)*f);
    
    K=(((EOP(r(2),3)-(EOP(r(1),3)))*VL)-((EOP(r(2),4)-(EOP(r(1),4)))*UL))/(VR*UL-VL*UR);
    
    ZL_Hor(i)=(K*WL)+(EOP(r(1),5));
    ZR_Hor(i)=(K*WR)+(EOP(r(2),5));
end

Zh=((ZL_Hor+ZR_Hor)/2)';

%For Elevation control Point
r(1)=3;
r(2)=4;
r1=37;
r2=49;
   
kappaL=EOP(r(1),8);
kappaR=EOP(r(2),8);
    
ML=[cos(kappaL), -sin(kappaL), 0;
    sin(kappaL), cos(kappaL) , 0;
    0          , 0           , 1];
    
MR=[cos(kappaR), -sin(kappaR), 0;
    sin(kappaR), cos(kappaR) , 0;
    0          , 0           , 1];
    
UL=(ML(1,1)*IP(r1,4))+(ML(2,1)*IP(r1,5))-(ML(3,1)*f);
VL=(ML(1,2)*IP(r1,4))+(ML(2,2)*IP(r1,5))-(ML(3,2)*f);
WL=(ML(1,3)*IP(r1,4))+(ML(2,3)*IP(r1,5))-(ML(3,3)*f);
    
UR=(MR(1,1)*IP(r2,4))+(MR(2,1)*IP(r2,5))-(MR(3,1)*f);
VR=(MR(1,2)*IP(r2,4))+(MR(2,2)*IP(r2,5))-(MR(3,2)*f);
WR=(MR(1,3)*IP(r2,4))+(MR(2,3)*IP(r2,5))-(MR(3,3)*f);
    
K=(((EOP(r(2),3)-(EOP(r(1),3)))*VL)-((EOP(r(2),4)-(EOP(r(1),4)))*UL))/(VR*UL-VL*UR);
   
XL_Ver=(K*UL)+(EOP(r(1),3));
YL_Ver=(K*VL)+(EOP(r(1),4));
XR_Ver=(K*UR)+(EOP(r(2),3));
YR_Ver=(K*VR)+(EOP(r(2),4));
    
Xv=((XL_Ver+XR_Ver)/2)';
Yv=((YL_Ver+YR_Ver)/2)';

GCP(4,4)=Zh(1);
GCP(10,4)=Zh(2);
GCP(6,2)=Xv(1);
GCP(6,3)=Yv(1);
%--------------------------------------------------------------------------

%(7)RMSE
for i=1:length(GCP)
    
    if GCP(i,2:4)~=0
       FP(i,:)=GCP(i,:);
    end
end

FP(4,:)=[];
FP(5,:)=[];
FP(8,:)=[];
dX=FP(:,2)-XR_Control';
dY=FP(:,3)-YR_Control';
dZ=FP(:,4)-ZL_Control';

R=sqrt(dX.^2 + dY.^2);
RMSE_2D=sqrt((R'*R)/(length(R)-1))
RMSE_Elev=sqrt((dZ'*dZ)/(length(dZ)-1))
%--------------------------------------------------------------------------

%(8)Plot
plot(Xt,Yt,'* blue');
hold on;
plot(XR_Control,YR_Control,'o red');
plot(GCP(:,2),GCP(:,3),'^ black');
legend('Tie Points','Computed Control Points','Control Points');
title('Ground Coordinates');
xlabel('X');
ylabel('Y');
%--------------------------------------------------------------------------

%(9)Export
result_coordinates=zeros(50,3);
for i=1:50
    result_coordinates(i,1:2)=[x_cap(14*4+2*i-1),x_cap(14*4+2*i)];
end

result_coordinates(:,3)=ZR_tie';

xlswrite('TieCoor_Approx.xls',result_coordinates);
xlswrite('EOP_Approx.xls',EOP);
xlswrite('GCP_Approx.xls',GCP);
xlswrite('IP_Approx.xls',IP);
%--------------------------------------------------------------------------

%(10)Table
Ran_num=IP(:,1);Image_num=IP(:,2);Point_code=IP(:,3);
Xi=IP(:,4);Yi=IP(:,5);Point_type=IP(:,6);Point_num=IP(:,7);
IP=table(Ran_num,Image_num,Point_code,Xi,Yi,Point_type,Point_num);

Point_code=GCP(:,1);
Xg=GCP(:,2);Yg=GCP(:,3);Zg=GCP(:,4);
GCP=table(Point_code,Xg,Yg,Zg);

Ran_num=EOP(:,1);Image_num=EOP(:,2);
XC=EOP(:,3);YC=EOP(:,4);ZC=EOP(:,5);
Omega=EOP(:,6);Phi=EOP(:,7);Kappa=EOPA(:,8);
EOP=table(Ran_num,Image_num,XC,YC,ZC,Omega,Phi,Kappa);

Tie_Ground_coordinates=table(Tie_code,Xt,Yt,Zt);
%--------------------------------------------------------------------------

clear Tie_code cran image_num Kappa Omega Phi Point_code Point_num Point_type Ran_num a b control_id_temp C c c1 c2 T Fu H V d data_control_points EOPA FP image_num i j k K kappa kappaL check1_temp ;
clear Tie_code Xc XC Xg Xi Xt Yc YC Yg Yi Yt Zc ZC Zg Zt check2_temp check3_temp check4_temp len_temp ML MR n omega p phi q r R r1 r2 UL UR VL VR WL WR X X0 x0 x Y y y0 Y0 Z z Z0 Z_avg kappaR landa;
%End
 