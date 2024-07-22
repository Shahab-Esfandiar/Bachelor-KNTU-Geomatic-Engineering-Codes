clc;
clear all;
close all;
format long g;

%project_7: Pointwise_transformation_figure2
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

%WDA_method  

%ECP_DistanceCalc
    r1=[];
    r2=[];
    r3=[];
    r4=[];
    for j=1:length(X)
        dx_w=Xc(1)-X(j);
        dy_w=Yc(1)-Y(j);
        Lw(j)=sqrt(dx_w.^2+dy_w.^2);
        
        if dx_w>0 & dy_w>0;
           r1(j)=Lw(j);
            
        else if dx_w<0 & dy_w>0;
                r2(j)=Lw(j);
                
        else if dx_w<0 & dy_w<0;
                r3(j)=Lw(j);
                
        else dx_w>0 & dy_w<0;
             r4(j)=Lw(j); 
        end
        end
        end
        
    end
    r1=sort(nonzeros(r1));
    r2=sort(nonzeros(r2));
    r3=sort(nonzeros(r3));
    r4=sort(nonzeros(r4));
    
    [m IX]=sort(Lw);
    m1=[r1(1) r2(1) r3(1) r4(1)];
    IX1=[2 1 10 8];
    
    
    w=1./(m1(1:4).^2);
    XI1=dx(IX1(1:4));
    YI1=dy(IX1(1:4));
    dx_ch1(1)=sum(w*XI1)/sum(w);
    dy_ch1(1)=sum(w*YI1)/sum(w);

    r1=[];
    r2=[];
    r3=[];
    r4=[];
    for j=1:length(X)
        dx_w=Xc(2)-X(j);
        dy_w=Yc(2)-Y(j);
        Lw(j)=sqrt(dx_w.^2+dy_w.^2);
        
        if dx_w>0 & dy_w>0;
           r1(j)=Lw(j);
            
        else if dx_w<0 & dy_w>0;
                r2(j)=Lw(j);
                
        else if dx_w<0 & dy_w<0;
                r3(j)=Lw(j);
                
        else dx_w>0 & dy_w<0;
             r4(j)=Lw(j); 
        end
        end
        end
        
    end
    r1=sort(nonzeros(r1));
    r2=sort(nonzeros(r2));
    r3=sort(nonzeros(r3));
    r4=sort(nonzeros(r4));
    
    [m IX]=sort(Lw);
    m2=[r1(1) r2(1) r3(1)];
    IX2=[4 3 5];
    
    
    w=1./(m2(1:3).^2);
    XI2=dx(IX2(1:3));
    YI2=dy(IX2(1:3));
    dx_ch1(2)=sum(w*XI2)/sum(w);
    dy_ch1(2)=sum(w*YI2)/sum(w);
    
    r1=[];
    r2=[];
    r3=[];
    r4=[];
    for j=1:length(X)
        dx_w=Xc(3)-X(j);
        dy_w=Yc(3)-Y(j);
        Lw(j)=sqrt(dx_w.^2+dy_w.^2);
        
        if dx_w>0 & dy_w>0;
           r1(j)=Lw(j);
            
        else if dx_w<0 & dy_w>0;
                r2(j)=Lw(j);
                
        else if dx_w<0 & dy_w<0;
                r3(j)=Lw(j);
                
        else dx_w>0 & dy_w<0;
             r4(j)=Lw(j); 
        end
        end
        end
        
    end
    r1=sort(nonzeros(r1));
    r2=sort(nonzeros(r2));
    r3=sort(nonzeros(r3));
    r4=sort(nonzeros(r4));
    
    [m IX]=sort(Lw);
    m3=[r2(1) r3(1) r4(1)];
    IX3=[1 7 2];
    
    
    w=1./(m3(1:3).^2);
    XI3=dx(IX3(1:3));
    YI3=dy(IX3(1:3));
    dx_ch1(3)=sum(w*XI3)/sum(w);
    dy_ch1(3)=sum(w*YI3)/sum(w);


%FinalCP_Coordinates
for i=1:length(Xc)
    Er(i,1)=Lch(i*2-1);
    Nr(i,1)=Lch(i*2);
    
    Xf1(i)=Er(i)+dx_ch1(i);
    Yf1(i)=Nr(i)+dy_ch1(i);
end

for i=1:3
ex(i)=Ec(i)-Xf1(i);
ey(i)=Nc(i)-Yf1(i);
end

et1=[ex;ey];

%RMSECalc
Vx=Xf1'-Ec;
Vy=Yf1'-Nc;
for i=1:length(Vx)
    dr(i)=sqrt(Vx(i)^2+Vy(i)^2);
end

RMSE1=sqrt(sum(dr.^2)/(length(Vx)-1));

%MA_method
for i=1:length(Xc)
    Lw=[];
    for j=1:length(X)
        dx_w=Xc(i)-X(j);
        dy_w=Yc(i)-Y(j);
        
        Lw(j)=sqrt(dx_w.^2+dy_w.^2);
    end
    [m IX]=sort(Lw);
    XI=dx(IX(1:4));
    YI=dy(IX(1:4));
    xi=E(IX(1:4));
    yi=N(IX(1:4));
    
    k=1;
for z=1:2:length(X)*2
    Am(z,:)=[X(k),Y(k),0,0,1,0];
    Am(z+1,:)=[0,0,X(k),Y(k),0,1];
    Lm(z,1)=dx(k);
    Lm(z+1,1)=dy(k);
    k=k+1;
end
    Xm_cap=inv(Am'*Am)*Am'*Lm;
    dx_ch2(i)=Xm_cap(1)*Er(i)+Xm_cap(2)*Nr(i)+Xm_cap(5);
    dy_ch2(i)=Xm_cap(3)*Er(i)+Xm_cap(4)*Nr(i)+Xm_cap(6);
end

%FinalCP_Coordinates
for i=1:length(Xc)
    
    Xf2(i)=Er(i)+dx_ch2(i);
    Yf2(i)=Nr(i)+dy_ch2(i);
end

for i=1:3
ex(i)=Ec(i)-Xf2(i);
ey(i)=Nc(i)-Yf2(i);
end

et2=[ex;ey];

%RMSECalc
Vx=Xf2'-Ec;
Vy=Yf2'-Nc;
for i=1:length(Vx)
    dr(i)=sqrt(Vx(i)^2+Vy(i)^2);
end

RMSE2=sqrt(sum(dr.^2)/(length(Vx)-1));

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
disp(" WDA_method")
disp("-----------------------------")
disp(" matrix of distance :")
disp("    ")
disp("dxU")
disp("    ")
disp(dx_ch1)
disp("    ")
disp("dyU")
disp("    ")
disp(dy_ch1)
disp("-----------------------------")
disp(" Final calculated coordinates :")
disp("    ")
disp([Xf1;Yf1])
disp("-----------------------------")
disp(" Error of calculated data and input data :")
disp("    ")
disp(et1)
disp("-----------------------------")
disp(" RMSE Measure =")
disp("    ")
disp(RMSE1)
disp("-----------------------------")
disp(" MA_method")
disp("-----------------------------")
disp(" matrix of distance :")
disp("    ")
disp("dxU")
disp("    ")
disp(dx_ch2)
disp("    ")
disp("dyU")
disp("    ")
disp(dy_ch2)
disp("-----------------------------")
disp(" Final calculated coordinates :")
disp("    ")
disp([Xf2;Yf2])
disp("-----------------------------")
disp(" Error of calculated data and input data :")
disp("    ")
disp(et2)
disp("-----------------------------")
disp(" RMSE Measure =")
disp("    ")
disp(RMSE2)

%plots
subplot(2,2,1);
plot(Xp,Yp,'b->');
title('image system');

subplot(2,2,2);
plot(Eg,Ng,'b->');
title('ground system');

%end
