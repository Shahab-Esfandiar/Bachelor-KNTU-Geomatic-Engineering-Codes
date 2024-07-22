clc;
clear all;
close all;
format long;

%Gruop : Esfandiar_Tahmasebi

landa=0.44:10^-6:0.64;
L=landa;
ngs=1+((287.604+(4.8863./(L.^2))+(0.0682./(L.^4)))*10^-6);

plot(L,ngs,'r-')
%end