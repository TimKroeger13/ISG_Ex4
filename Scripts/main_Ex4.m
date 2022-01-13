%--------------------------------------------------------------------------
%   
%          Introduction to Space Geodesy - Main
%   Assignment 4: Reference systems and transformations
% 
%   Author         : Group B
%   Version        : January 13, 2021
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format longG

%% Task 1







%yyyy=1990:2030
yyyy=2000;
mm=1;
dd=1;
ut1=12;
minute=0;
second=0;

[jd,~] = gre2jd(yyyy,mm,dd,ut1,minute,second); %Julian Day
d = jd -2451545.0; %Based on year 2000
T = d/36525; %2000 based Julian centurie

zeta_v1 = De_Mi_Se_To_DecDeg(Expand_DecDeg(2306.2181));
zeta_v2 = De_Mi_Se_To_DecDeg(Expand_DecDeg(0.30188));
UpperZeta_v1 =  De_Mi_Se_To_DecDeg(Expand_DecDeg(2306.2181));
UpperZeta_v2 =  De_Mi_Se_To_DecDeg(Expand_DecDeg(1.09468));
teta_v1 =  De_Mi_Se_To_DecDeg(Expand_DecDeg(2004.3109));
teta_v2 =  De_Mi_Se_To_DecDeg(Expand_DecDeg(0.42665));


zeta= zeta_v1 * T + zeta_v2 * T^2;
UpperZeta= UpperZeta_v1 * T + UpperZeta_v2 * T^2;
teta= teta_v1 * T + teta_v2 * T^2;



alpha = 0;
delta = 0;

[x,y,z] = SphToCart(1,delta,alpha);

rot3d(-z)















function [r,lamda , phi] = CartToSphe(x,y,z)
format long
r= sqrt(x^2 +y^2 +z^2);
lamda_rad=atan2(y,x); %rad
lamda=lamda_rad*12/pi; %rad to hour angle
%we are using arctan2 because
%it can evaluate the domain between -pi to pi
phi_rad=atan(z/(sqrt(x^2+y^2))); %rad
phi= rad2deg(phi_rad) ; %degree
end








%this function transfor spherical coorsinates to cartesian
% everything should be in radian
function [x,y,z] = SphToCart(r,phi,lambda)
phi=0
lambda=0
r=1
x=r*cos(lambda)*cos(phi);
y= r*sin(lambda)*cos(phi);
z=r*sin(phi);
end






function [P] = perc(JD)

d=JD - 2451545;
T= d/36525;

a1=[0 38 26.2181]; a1degree=dms2degrees(a1); %2306".2181
a2=[0 0 0.30188]; a2degree=dms2degrees(a2);
a3=[0 0 1.09468]; a3degree=dms2degrees(a3);
a4=[0 33 24.3109]; a4degree=dms2degrees(a4); %2004.3109"
a5=[0 0 0.42665]; a5degree=dms2degrees(a5);

zeta= a1degree* T + a2degree *T^2;
z= a1degree* T + a3degree *T^2;
teta = a4degree * T + a5degree *T^2;







R2=[cosd(teta) 0 -sind(teta)
    0 1 0
    sind(teta) 0 cosd(teta)];

R3_z =[cosd(-z) sind(-z) 0
    -sind(-z) cosd(-z) 0
    0 0 1];

R3_zeta =[cosd(-zeta) sind(-zeta) 0
    -sind(-zeta) cosd(-zeta) 0
    0 0 1];
P=R3_z * R2 * R3_zeta;

end









