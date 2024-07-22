clc;
clear all;
close all;
format long;

%ShahabEsfandiar_9819373

%Data
t=[13.23,14.26,15.21,16.24,17.27,15.29,16.23,14.29,15.21,15.22,13.24,17.25];
p=[1026.47,1026.43,1027.44,1024.41,1026.42,1026.48,1024.49,1026.43,1027.46,1027.48,1027.41,1026.42];
e=[17.063,18.062,19.061,17.066,18.067,16.065,17.067,18.068,16.063,17.069,17.065,19.064];
t_avg=15.20;
p_avg=1026.40;
e_avg=17.060;
A=0.283617;
B=0.0016062;
C=0.00001345;
Q=0.0468;
R=0.00051;
ngs=1.0003;
alpha=0.003661;
L_vs=0.52;
beta=(0.7868 - (0.118*t_avg))*10^-6;

%Solve N for visible spectrum

%p,e=const and t=unknown
Nvs_t=1+((((p_avg*((1+beta*p_avg)./(1+alpha*t)))*(A+(B/L_vs^2)+(C/L_vs^4)))-((e_avg./1+alpha*t)*(Q-(R/L_vs^2))))*10^-6);

%e,t=const and p=unknown
Nvs_p=1+((((p.*((1+beta.*p)/(1+alpha*t_avg)))*(A+(B/L_vs^2)+(C/L_vs^4)))-((e_avg/1+alpha*t_avg)*(Q-(R/L_vs^2))))*10^-6);

%t,p=const and e=unknown
Nvs_e=1+((((p_avg*((1+beta*p_avg)/(1+alpha*t_avg)))*(A+(B/L_vs^2)+(C/L_vs^4)))-((e./1+alpha*t_avg)*(Q-(R/L_vs^2))))*10^-6);


%Solve N for Microwave

%p,e=const and t=unknown
Nmw_t=1+((((p_avg-e_avg).*(77.624./273.16+t))+((e_avg*64.7./273.16+t).*(1+(5748./273.16+t))))*10^-6);

%e,t=const and p=unknown
Nmw_p=1+((((p-e_avg).*(77.624/273.16+t_avg))+((e_avg*64.7/273.16+t_avg)*(1+(5748/273.16+t_avg))))*10^-6);

%t,p=const and e=unknown
Nmw_e=1+((((p_avg-e).*(77.624/273.16+t_avg))+((e.*64.7/273.16+t_avg)*(1+(5748/273.16+t_avg))))*10^-6);


%Solve N for Radiofrequency

%p,e=const and t=unknown
Nrf_t=1+((ngs-1./1+alpha*t)*(p_avg/1013.25))-(0.042*10^-6*e_avg./1+alpha*t);

%e,t=const and p=unknown
Nrf_p=1+((ngs-1/1+alpha*t_avg).*(p/1013.25))-(0.042*10^-6*e_avg/1+alpha*t_avg);

%t,p=const and e=unknown
Nrf_e=1+((ngs-1/1+alpha*t_avg)*(p_avg/1013.25))-(0.042*10^-6*e./1+alpha*t_avg);


%plot

%for visible spectrum
figure
title('visible spectrum')
subplot(1,3,1);plot(t,Nvs_t,'b-')
xlabel('temperature')
ylabel('Refractive index')
subplot(1,3,2);plot(p,Nvs_p,'r-')
xlabel('air pressure')
subplot(1,3,3);plot(e,Nvs_e,'g-')
xlabel('Relative pressure of watervapor')

%for Microwave
figure
title('Microwave')
subplot(1,3,1);plot(t,Nmw_t,'b-')
xlabel('temperature')
ylabel('Refractive index')
subplot(1,3,2);plot(p,Nmw_p,'r-')
xlabel('air pressure')
subplot(1,3,3);plot(e,Nmw_e,'g-')
xlabel('Relative pressure of watervapor')

%for Radiofrequency
figure
subplot(1,3,1);plot(t,Nrf_t,'b-')
xlabel('temperature')
ylabel('Refractive index')
subplot(1,3,2);plot(p,Nrf_p,'r-')
xlabel('air pressure')
subplot(1,3,3);plot(e,Nrf_e,'g-')
xlabel('Relative pressure of watervapor')

%end

