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
%yyyy=2000;
yyyy=1990;
mm=1;
dd=1;
ut1=12;
minute=0;
second=0;

[jd,~] = gre2jd(yyyy,mm,dd,ut1,minute,second); %Julian Day
P = precession(jd); %precession

alpha = 0;
delta = 0;


[x,y,z] = SphToCart(1,delta,alpha);

precCart = P * [x,y,z]';

[r,alpha , delta] = CartToSphe(precCart(1),precCart(2),precCart(3));

[De_Mi_Se_Matrix] = DecDeg_To_De_Mi_Se(alpha);

alpha = sum(De_Mi_Se_Matrix(1)+De_Mi_Se_Matrix(2)/60+De_Mi_Se_Matrix(3)/3600);

%% Exercise1 Loop

yyyy=1990:2030;
mm=1;
dd=1;
ut1=12;
minute=0;
second=0;

% coordinatesofa fictitious celestial object at the vernal equinox (0; 0) 
J2000_alpha = 0; 
J2000_delta = 0;

%[jd_2000,~] = gre2jd(2000,1,1,12,0,0); %J2000

[x,y,z] = SphToCart(1,J2000_delta,J2000_alpha);

counter = 0;

Calculated_delta = zeros(1,length(yyyy));
Calculated_alpha = zeros(1,length(yyyy));

for i = yyyy

counter = counter+1;

[jd,~] = gre2jd(i,mm,dd,ut1,minute,second); %Running Julian Day

P = precession(jd); %precession

precCart = P * [x,y,z]';

[r,alpha , delta] = CartToSphe(precCart(1),precCart(2),precCart(3));

Calculated_delta(counter) = delta;
Calculated_alpha(counter) = alpha;

end

%This ia bad, but i don't know right now the fix the offet or if
%the calculations are even right here
a = Calculated_alpha/180*12
a(a < -6) = a(a < -6)+12;


plotyy(yyyy,Calculated_delta,yyyy,a)






