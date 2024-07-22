clc;
clear all;
close all;
format long g;

%Project4_PreProgramCase4_ApplicationalPhotogrametery
%ShahabEsfandiar_9819373
%HomaGangali_9929953
%--------------------------------------------------------------------------

%(0)Data
DataPoints=xlsread('data_main.xlsx');
GCP=importdata('Control.txt');
sFull=importdata('full.txt');
%--------------------------------------------------------------------------

%(1)Find the full control points of model's corners
j=1;

for i=1:length(GCP)
    
    if GCP(i,5)==1
       Full(j,:)=GCP(i,:);
       j=j+1;
    end
end
%--------------------------------------------------------------------------

%(2)Find Tie and control points for our new model
k=1;
q=1;
p=1;
DataPoints(:,8)=[];
DataPoints(:,9)=[];

for i=1:length(DataPoints)
    
    if DataPoints(i,6)==0 || DataPoints(i,6)==1
        DupData(k,:)=DataPoints(i,:);
        k=k+1;
    end
end

% Selected full control Points : 1109,1709,2189,2789
for i=1:length(DupData) 
    
    if DupData(i,3)~=1128 && DupData(i,3)~=1309 && DupData(i,3)~=1648 && DupData(i,3)~=2389
       DataNew(q,:)=DupData(i,:);
       q=q+1;
    end
end

% Numbering full control Points
DataNew(:,8)=0; 

for i=1:length(DataNew)
    
    if (DataNew(i,6)==1 && DataNew(i,8)==0)
          
        for j=1:length(DataNew)

            if (DataNew(i,3)==DataNew(j,3) &&  DataNew(j,8)==0)
                DataNew(j,8)=p;     
            end   
        end
        
    DataNew(i,8)=p;
    p=p+1;
    end
end

FID=Full(:,1);
xFull=Full(:,2);
yFull=Full(:,3);
XsFull=sFull(:,2);
YsFull=sFull(:,3);
%--------------------------------------------------------------------------

%(3)Plot
plot(xFull,yFull,'* blue');
text(xFull,yFull,int2str(FID));
hold on;
plot(XsFull,YsFull,'* red');
legend('Full Control Points','Selected Control Points');
title('Ground coordinates');
xlabel('X');
ylabel('Y');
%--------------------------------------------------------------------------

%(4)Export
xlswrite('PreFull.xls',Full);
xlswrite('sFull.xls',sFull);
xlswrite('DataNew.xls',DataNew);
%--------------------------------------------------------------------------

clear DupData FID i j k p q xFull XsFull yFull YsFull
%end


