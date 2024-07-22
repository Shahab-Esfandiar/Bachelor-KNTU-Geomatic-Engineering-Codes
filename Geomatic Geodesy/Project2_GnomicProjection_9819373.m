clc;
clear all;
close all;
format long g;

%Project2_GeometricalGeodesy
%ShahabEsfandiar_9819373
%HomaGangali_9929953
%--------------------------------------------------------------------------

%Observed Data        
l=[81.139;101.675;124.928;80.082;52.005;40.935;
   81.162;94.543;148.659;151.476;92.075;51.484;
   101.677;94.531;60.452;116.592;57.774;68.634;
   124.928;148.647;60.458;95.905;73.299;108.466;
   80.081;151.438;116.585;95.905;66.25;100.999;
   52.003;92.056;57.773;73.297;66.255;42;
   40.936;51.472;68.636;108.464;100.299;42];

a=dms2degrees([84 28 32;23 24 01;354 44 19;304 35 41;52 05 03;
               302 50 04;13 04 23;322 44 46;337 02 10;328 02 46;
               339 03 02;290 21 56;76 35 33;21 24 47;322 25 30;
               53 46 59;20 42 28;93 38 55;50 03 12;35 15 34;
               89 52 30;69 47 16;31 09 40;49 44 07;67 04 56;
               298 41 12;224 21 24;171 0 29;84 21 17;309 44 57;
               282 21 10;159 55 59;57 02 19;26 28 10;331 56 21]);
%--------------------------------------------------------------------------

%Equatios
syms Xcoor1 Xcoor2 Xcoor3 Xcoor4 Xcoor5 Xcoor6 Xcoor7 ;
syms Ycoor1 Ycoor2 Ycoor3 Ycoor4 Ycoor5 Ycoor6 Ycoor7 ;
syms Zcoor1 Zcoor2 Zcoor3 Zcoor4 Zcoor5 Zcoor6 Zcoor7 ;

f(1) = sqrt((Xcoor2-Xcoor1)^2+(Ycoor2-Ycoor1)^2+(Zcoor2-Zcoor1)^2);
f(2) = sqrt((Xcoor3-Xcoor1)^2+(Ycoor3-Ycoor1)^2+(Zcoor3-Zcoor1)^2);
f(3) = sqrt((Xcoor4-Xcoor1)^2+(Ycoor4-Ycoor1)^2+(Zcoor4-Zcoor1)^2);
f(4) = sqrt((Xcoor5-Xcoor1)^2+(Ycoor5-Ycoor1)^2+(Zcoor5-Zcoor1)^2);
f(5) = sqrt((Xcoor6-Xcoor1)^2+(Ycoor6-Ycoor1)^2+(Zcoor6-Zcoor1)^2);
f(6) = sqrt((Xcoor7-Xcoor1)^2+(Ycoor7-Ycoor1)^2+(Zcoor7-Zcoor1)^2);

f(7) = sqrt((Xcoor1-Xcoor2)^2+(Ycoor1-Ycoor2)^2+(Zcoor1-Zcoor2)^2);
f(8) = sqrt((Xcoor3-Xcoor2)^2+(Ycoor3-Ycoor2)^2+(Zcoor3-Zcoor2)^2);
f(9) = sqrt((Xcoor4-Xcoor2)^2+(Ycoor4-Ycoor2)^2+(Zcoor4-Zcoor2)^2);
f(10) = sqrt((Xcoor5-Xcoor2)^2+(Ycoor5-Ycoor2)^2+(Zcoor5-Zcoor2)^2);
f(11) = sqrt((Xcoor6-Xcoor2)^2+(Ycoor6-Ycoor2)^2+(Zcoor6-Zcoor2)^2);
f(12) = sqrt((Xcoor7-Xcoor2)^2+(Ycoor7-Ycoor2)^2+(Zcoor7-Zcoor2)^2);

f(13) = sqrt((Xcoor1-Xcoor3)^2+(Ycoor1-Ycoor3)^2+(Zcoor1-Zcoor3)^2);
f(14) = sqrt((Xcoor2-Xcoor3)^2+(Ycoor2-Ycoor3)^2+(Zcoor2-Zcoor3)^2);
f(15) = sqrt((Xcoor4-Xcoor3)^2+(Ycoor4-Ycoor3)^2+(Zcoor4-Zcoor3)^2);
f(16) = sqrt((Xcoor5-Xcoor3)^2+(Ycoor5-Ycoor3)^2+(Zcoor5-Zcoor3)^2);
f(17) = sqrt((Xcoor6-Xcoor3)^2+(Ycoor6-Ycoor3)^2+(Zcoor6-Zcoor3)^2);
f(18) = sqrt((Xcoor7-Xcoor3)^2+(Ycoor7-Ycoor3)^2+(Zcoor7-Zcoor3)^2);

f(19) = sqrt((Xcoor1-Xcoor4)^2+(Ycoor1-Ycoor4)^2+(Zcoor1-Zcoor4)^2);
f(20) = sqrt((Xcoor2-Xcoor4)^2+(Ycoor2-Ycoor4)^2+(Zcoor2-Zcoor4)^2);
f(21) = sqrt((Xcoor3-Xcoor4)^2+(Ycoor3-Ycoor4)^2+(Zcoor3-Zcoor4)^2);
f(22) = sqrt((Xcoor5-Xcoor4)^2+(Ycoor5-Ycoor4)^2+(Zcoor5-Zcoor4)^2);
f(23) = sqrt((Xcoor6-Xcoor4)^2+(Ycoor6-Ycoor4)^2+(Zcoor6-Zcoor4)^2);
f(24) = sqrt((Xcoor7-Xcoor4)^2+(Ycoor7-Ycoor4)^2+(Zcoor7-Zcoor4)^2);

f(25) = sqrt((Xcoor1-Xcoor5)^2+(Ycoor1-Ycoor5)^2+(Zcoor1-Zcoor5)^2);
f(26) = sqrt((Xcoor2-Xcoor5)^2+(Ycoor2-Ycoor5)^2+(Zcoor2-Zcoor5)^2);
f(27) = sqrt((Xcoor3-Xcoor5)^2+(Ycoor3-Ycoor5)^2+(Zcoor3-Zcoor5)^2);
f(28) = sqrt((Xcoor4-Xcoor5)^2+(Ycoor4-Ycoor5)^2+(Zcoor4-Zcoor5)^2);
f(29) = sqrt((Xcoor6-Xcoor5)^2+(Ycoor6-Ycoor5)^2+(Zcoor6-Zcoor5)^2);
f(30) = sqrt((Xcoor7-Xcoor5)^2+(Ycoor7-Ycoor5)^2+(Zcoor7-Zcoor5)^2);

f(31) = sqrt((Xcoor1-Xcoor6)^2+(Ycoor1-Ycoor6)^2+(Zcoor1-Zcoor6)^2);
f(32) = sqrt((Xcoor2-Xcoor6)^2+(Ycoor2-Ycoor6)^2+(Zcoor2-Zcoor6)^2);
f(33) = sqrt((Xcoor3-Xcoor6)^2+(Ycoor3-Ycoor6)^2+(Zcoor3-Zcoor6)^2);
f(34) = sqrt((Xcoor4-Xcoor6)^2+(Ycoor4-Ycoor6)^2+(Zcoor4-Zcoor6)^2);
f(35) = sqrt((Xcoor5-Xcoor6)^2+(Ycoor5-Ycoor6)^2+(Zcoor5-Zcoor6)^2);
f(36) = sqrt((Xcoor7-Xcoor6)^2+(Ycoor7-Ycoor6)^2+(Zcoor7-Zcoor6)^2);

f(37) = sqrt((Xcoor1-Xcoor7)^2+(Ycoor1-Ycoor7)^2+(Zcoor1-Zcoor7)^2);
f(38) = sqrt((Xcoor2-Xcoor7)^2+(Ycoor2-Ycoor7)^2+(Zcoor2-Zcoor7)^2);
f(39) = sqrt((Xcoor3-Xcoor7)^2+(Ycoor3-Ycoor7)^2+(Zcoor3-Zcoor7)^2);
f(40) = sqrt((Xcoor4-Xcoor7)^2+(Ycoor4-Ycoor7)^2+(Zcoor4-Zcoor7)^2);
f(41) = sqrt((Xcoor5-Xcoor7)^2+(Ycoor5-Ycoor7)^2+(Zcoor5-Zcoor7)^2);
f(42) = sqrt((Xcoor6-Xcoor7)^2+(Ycoor6-Ycoor7)^2+(Zcoor6-Zcoor7)^2);

f(43) = acosd(([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]*[Xcoor2 - Xcoor1, Ycoor2 - Ycoor1]')/ (norm([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]) * norm([Xcoor2 - Xcoor1, Ycoor2 - Ycoor1])));
f(44) = acosd(([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]*[Xcoor3 - Xcoor1, Ycoor3 - Ycoor1]')/ (norm([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]) * norm([Xcoor3 - Xcoor1, Ycoor3 - Ycoor1])));
f(45) = 360-acosd(([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]*[Xcoor4 - Xcoor1, Ycoor4 - Ycoor1]')/ (norm([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]) * norm([Xcoor4 - Xcoor1, Ycoor4 - Ycoor1])));
f(46) = 360-acosd(([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]*[Xcoor5 - Xcoor1, Ycoor5 - Ycoor1]')/ (norm([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]) * norm([Xcoor5 - Xcoor1, Ycoor5 - Ycoor1])));
f(47) = acosd(([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]*[Xcoor7 - Xcoor1, Ycoor7 - Ycoor1]')/ (norm([Xcoor6 - Xcoor1, Ycoor6 - Ycoor1]) * norm([Xcoor7 - Xcoor1, Ycoor7 - Ycoor1])));

f(48) = 360-acosd(([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]*[Xcoor1 - Xcoor2, Ycoor1 - Ycoor2]')/ (norm([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]) * norm([Xcoor1 - Xcoor2, Ycoor1 - Ycoor2])));
f(49) = acosd(([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]*[Xcoor3 - Xcoor2, Ycoor3 - Ycoor2]')/ (norm([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]) * norm([Xcoor3 - Xcoor2, Ycoor3 - Ycoor2])));
f(50) = 360-acosd(([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]*[Xcoor5 - Xcoor2, Ycoor5 - Ycoor2]')/ (norm([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]) * norm([Xcoor5 - Xcoor2, Ycoor5 - Ycoor2])));
f(51) = 360-acosd(([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]*[Xcoor6 - Xcoor2, Ycoor6 - Ycoor2]')/ (norm([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]) * norm([Xcoor6 - Xcoor2, Ycoor6 - Ycoor2])));
f(52) = 360-acosd(([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]*[Xcoor7 - Xcoor2, Ycoor7 - Ycoor2]')/ (norm([Xcoor4 - Xcoor2, Ycoor4 - Ycoor2]) * norm([Xcoor7 - Xcoor2, Ycoor7 - Ycoor2])));

f(53) = 360-acosd(([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]*[Xcoor1 - Xcoor3, Ycoor1 - Ycoor3]')/ (norm([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]) * norm([Xcoor1 - Xcoor3, Ycoor1 - Ycoor3])));
f(54) = 360-acosd(([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]*[Xcoor2 - Xcoor3, Ycoor2 - Ycoor3]')/ (norm([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]) * norm([Xcoor2 - Xcoor3, Ycoor2 - Ycoor3])));
f(55) = acosd(([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]*[Xcoor4 - Xcoor3, Ycoor4 - Ycoor3]')/ (norm([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]) * norm([Xcoor4 - Xcoor3, Ycoor4 - Ycoor3])));
f(56) = acosd(([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]*[Xcoor5 - Xcoor3, Ycoor5 - Ycoor3]')/ (norm([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]) * norm([Xcoor5 - Xcoor3, Ycoor5 - Ycoor3])));
f(57) = 360-acosd(([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]*[Xcoor7 - Xcoor3, Ycoor7 - Ycoor3]')/ (norm([Xcoor6 - Xcoor3, Ycoor6 - Ycoor3]) * norm([Xcoor7 - Xcoor3, Ycoor7 - Ycoor3])));

f(58) = acosd(([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]*[Xcoor1 - Xcoor4, Ycoor1 - Ycoor4]')/ (norm([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]) * norm([Xcoor1 - Xcoor4, Ycoor1 - Ycoor4])));
f(59) = acosd(([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]*[Xcoor2 - Xcoor4, Ycoor2 - Ycoor4]')/ (norm([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]) * norm([Xcoor2 - Xcoor4, Ycoor2 - Ycoor4])));
f(60) = acosd(([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]*[Xcoor5 - Xcoor4, Ycoor5 - Ycoor4]')/ (norm([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]) * norm([Xcoor5 - Xcoor4, Ycoor5 - Ycoor4])));
f(61) = acosd(([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]*[Xcoor6 - Xcoor4, Ycoor6 - Ycoor4]')/ (norm([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]) * norm([Xcoor6 - Xcoor4, Ycoor6 - Ycoor4])));
f(62) = acosd(([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]*[Xcoor7 - Xcoor4, Ycoor7 - Ycoor4]')/ (norm([Xcoor3 - Xcoor4, Ycoor3 - Ycoor4]) * norm([Xcoor7 - Xcoor4, Ycoor7 - Ycoor4])));

f(63) = acosd(([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]*[Xcoor1 - Xcoor5, Ycoor1 - Ycoor5]')/ (norm([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]) * norm([Xcoor1 - Xcoor5, Ycoor1 - Ycoor5])));
f(64) = acosd(([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]*[Xcoor2 - Xcoor5, Ycoor2 - Ycoor5]')/ (norm([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]) * norm([Xcoor2 - Xcoor5, Ycoor2 - Ycoor5])));
f(65) = acosd(([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]*[Xcoor3 - Xcoor5, Ycoor3 - Ycoor5]')/ (norm([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]) * norm([Xcoor3 - Xcoor5, Ycoor3 - Ycoor5])));
f(66) = acosd(([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]*[Xcoor6 - Xcoor5, Ycoor6 - Ycoor5]')/ (norm([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]) * norm([Xcoor6 - Xcoor5, Ycoor6 - Ycoor5])));
f(67) = acosd(([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]*[Xcoor7 - Xcoor5, Ycoor7 - Ycoor5]')/ (norm([Xcoor4 - Xcoor5, Ycoor4 - Ycoor5]) * norm([Xcoor7 - Xcoor5, Ycoor7 - Ycoor5])));

f(68) = 360-acosd(([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]*[Xcoor2 - Xcoor6, Ycoor2 - Ycoor6]')/ (norm([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]) * norm([Xcoor2 - Xcoor6, Ycoor2 - Ycoor6])));
f(69) = 360-acosd(([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]*[Xcoor3 - Xcoor6, Ycoor3 - Ycoor6]')/ (norm([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]) * norm([Xcoor3 - Xcoor6, Ycoor3 - Ycoor6])));
f(70) = acosd(([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]*[Xcoor4 - Xcoor6, Ycoor4 - Ycoor6]')/ (norm([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]) * norm([Xcoor4 - Xcoor6, Ycoor4 - Ycoor6])));
f(71) = acosd(([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]*[Xcoor5 - Xcoor6, Ycoor5 - Ycoor6]')/ (norm([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]) * norm([Xcoor5 - Xcoor6, Ycoor5 - Ycoor6])));
f(72) = 360-acosd(([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]*[Xcoor7 - Xcoor6, Ycoor7 - Ycoor6]')/ (norm([Xcoor1 - Xcoor6, Ycoor1 - Ycoor6]) * norm([Xcoor7 - Xcoor6, Ycoor7 - Ycoor6])));

f(73) = 360-acosd(([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]*[Xcoor1 - Xcoor7, Ycoor1 - Ycoor7]')/ (norm([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]) * norm([Xcoor1 - Xcoor7, Ycoor1 - Ycoor7])));
f(74) = acosd(([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]*[Xcoor2 - Xcoor7, Ycoor2 - Ycoor7]')/ (norm([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]) * norm([Xcoor2 - Xcoor7, Ycoor2 - Ycoor7])));
f(75) = acosd(([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]*[Xcoor3 - Xcoor7, Ycoor3 - Ycoor7]')/ (norm([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]) * norm([Xcoor3 - Xcoor7, Ycoor3 - Ycoor7])));
f(76) = acosd(([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]*[Xcoor4 - Xcoor7, Ycoor4 - Ycoor7]')/ (norm([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]) * norm([Xcoor4 - Xcoor7, Ycoor4 - Ycoor7])));
f(77) = 360-acosd(([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]*[Xcoor5 - Xcoor7, Ycoor5 - Ycoor7]')/ (norm([Xcoor6 - Xcoor7, Ycoor6 - Ycoor7]) * norm([Xcoor5 - Xcoor7, Ycoor5 - Ycoor7])));

A=jacobian(f,[Xcoor1 Xcoor2 Xcoor3 Xcoor4 Xcoor5 Xcoor6 Xcoor7 Ycoor1 Ycoor2 Ycoor3 Ycoor4 Ycoor5 Ycoor6 Ycoor7 Zcoor1 Zcoor2 Zcoor3 Zcoor4 Zcoor5 Zcoor6 Zcoor7]);
L=[l;a];
%--------------------------------------------------------------------------

%Inital coordinates
Xcoor1=561741;
Ycoor1=3607310;
Zcoor1=1724;

Xcoor2=561751.9;
Ycoor2=3607230;
Zcoor2=1719.04;

Xcoor3=561835.8;
Ycoor3=3607273;
Zcoor3=1722.917;

Xcoor4=561864.9;
Ycoor4=3607326;
Zcoor4=1727.232;

Xcoor5=561784;
Ycoor5=3607378;
Zcoor5=1729.57;

Xcoor6=561793;
Ycoor6=3607312;
Zcoor6=1724.959;

Xcoor7=561767.3;
Ycoor7=3607279;
Zcoor7=1722.358;      
%--------------------------------------------------------------------------

%Least Square adjustment
for i=1:100
    
    F=eval(f);
    dl=L-F';
    An=eval(A);
    dx=pinv(An'*An)*An'*dl;
    
    Xcoor1 = Xcoor1 + dx(1);
    Xcoor2 = Xcoor2 + dx(2);
    Xcoor3 = Xcoor3 + dx(3);
    Xcoor4 = Xcoor4 + dx(4);
    Xcoor5 = Xcoor5 + dx(5);
    Xcoor6 = Xcoor6 + dx(6);
    Xcoor7 = Xcoor7 + dx(7);
    Ycoor1 = Ycoor1 + dx(8);
    Ycoor2 = Ycoor2 + dx(9);
    Ycoor3 = Ycoor3 + dx(10);
    Ycoor4 = Ycoor4 + dx(11);
    Ycoor5 = Ycoor5 + dx(12);
    Ycoor6 = Ycoor6 + dx(13);
    Ycoor7 = Ycoor7 + dx(14);
    Zcoor1 = Zcoor1 + dx(15);
    Zcoor2 = Zcoor2 + dx(16);
    Zcoor3 = Zcoor3 + dx(17);
    Zcoor4 = Zcoor4 + dx(18);
    Zcoor5 = Zcoor5 + dx(19);
    Zcoor6 = Zcoor6 + dx(20);
    Zcoor7 = Zcoor7 + dx(21);
     
    if abs(dx)<10^-3
        break
    end
end

xcap=[Xcoor1;Xcoor2;Xcoor3;Xcoor4;Xcoor5;Xcoor6;Xcoor7;Ycoor1;Ycoor2;Ycoor3;Ycoor4;Ycoor5;Ycoor6;Ycoor7;Zcoor1;Zcoor2;Zcoor3;Zcoor4;Zcoor5;Zcoor6;Zcoor7];
x=xcap(1:7);
y=xcap(8:14);
z=xcap(15:21);
%--------------------------------------------------------------------------

%Convert to Geodetic coordinates
f=1/298.257223563;
a=6378137;
b=-((a*f)-a);
e=(sqrt((a^2)-(b^2)))/a;
e2=e^2/(1-(e^2));
e1=(1-sqrt(1-(e^2)))/(1+sqrt(1-(e^2)));
k0=0.9996;
zone=39;

if zone>30
    la0=((zone-30)*6)-3;
else
    la0=(-(zone-30)*6)-3;
end

M=(x./k0);
mu=M./(a*(1-((e^2)/4)-((3*(e^4))/64)-((5*(e^6))/256)));
phi1=mu+(((3*e1)/2)-((27*(e1^3))/32))*sin(2*mu)+(((21*(e1^2))/16)-((55*(e1^4))/32))*sin(4*mu)+(((151*(e1^3))/96))*sin(6*mu);
T1=tan(phi1);
C1=e2*(cos(phi1)).^2; 
M1=sqrt((a*(1-(e^2)))./(1-(e^2)*(sin(phi1)).^2).^3);
N1=a./sqrt(1-(e^2)*(sin(phi1)).^2);

phi=phi1+((y.^2).*T1./((k0^2)*2.*M1.*N1)-((y.^4).*T1./((k0^4)*24.*M1.*(N1.^3))).*(5+(3*(T1.^2))+(6*C1)-(6*C1.*(T1.^2))-(3*(C1.^4))-(9*(T1.^2).*(C1.^2)))+(((y.^6).*T1)./((k0^6)*720.*M1.*(N1.^5))).*(61+90*(T1.^2)+45*(T1.^4)+107*C1-162*(T1.^2).*C1)-45*(T1.^4).*C1); 
landa=la0+(y./(k0*N1.*cos(phi1)))-((y.^3)./((k0^3)*6*(N1.^3).*cos(phi1))).*(1+2*(T1.^2)+C1)+((y.^5)./((k0^5)*120*(N1.^5).*cos(phi1))).*(5+28*(T1.^2)+24*(T1.^4)+6*C1+8*(T1.^2).*C1);
%--------------------------------------------------------------------------

%Gnomonic projection
x_gnomonic=(cosd(phi).*sind(landa-51.6660))./((sind(32.6539)-sind(landa))+(cosd(32.6539).*cosd(phi).*cosd(landa-51.6660)));
y_gnomonic=((cosd(32.6539).*sind(phi))-(cosd(32.6539).*cosd(phi).*cosd(landa-51.6660)))./((sind(32.6539)-sind(landa))+(cosd(32.6539).*cosd(phi).*cosd(landa-51.6660))); 
%--------------------------------------------------------------------------

%Table&Plot
n=[1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 4 4 4 5 5 6];
s=[2 3 4 5 6 7 3 4 5 6 7 4 5 6 7 5 6 7 6 7 7];
g=graph(n,s);

connect=[1 2;1 3;1 4;1 5;1 6;1 7;2 3;2 4;2 5;2 6;2 7;3 4;3 5;3 6;3 7;4 5;4 6;4 7;5 6;5 7;6 7];
ID=[1,2,3,4,5,6,7];
id=ID;
temp=[x_gnomonic,y_gnomonic,ID'];
xg=temp(:,1);
yg=temp(:,2);
FID=temp(:,3);

figure;
plot(g,'xdata',x,'ydata',y);
title('UTM Geodetic Network');
xlabel('X');
ylabel('Y');

figure;
plot(x_gnomonic,y_gnomonic,'ro')
text(xg,yg,int2str(FID));
title('Gnomonic Geodetic Network');
xlabel('X');
ylabel('Y');

hold on; 
for i = 1:size(connect, 1)
    plot([x_gnomonic(connect(i, 1)), x_gnomonic(connect(i, 2))], [y_gnomonic(connect(i, 1)), y_gnomonic(connect(i, 2))], 'b-'); 
end
hold off; 

figure;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
axesm ('gnomonic', 'Frame', 'on', 'Grid', 'on');
geoshow(landareas,'FaceColor',[1 1 .5],'EdgeColor',[.6 .6 .6]);
tissot;

UTM=table(id',x,y,z);
Geo=table(id',landa,phi);
Gnom=table(id',xg,yg);
%--------------------------------------------------------------------------

%Display
disp(' The exact UTM coordinates of our Points is : ')
fprintf('\n')
disp(UTM)
disp('--------------------------------------------- ')
fprintf('\n')

disp(' The Geodetic coordinates of our Points is : ')
fprintf('\n')
disp(Geo)
disp('--------------------------------------------- ')
fprintf('\n')

disp(' The Gnomonic coordinates of our Points is : ')
fprintf('\n')
disp(Gnom)

%End