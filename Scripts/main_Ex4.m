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

%% Task 1 - Precession

yyyy=1990:2030;
mm=1;
dd=1;
ut1=12;
minute=0;
second=0;

% coordinatesofa fictitious celestial object at the vernal equinox (0; 0) 
J2000_alpha = 0; 
J2000_delta = 0;

[x,y,z] = SphToCart(1,J2000_delta,J2000_alpha); 

counter = 0;

Calculated_delta = zeros(1,length(yyyy));
Calculated_alpha = zeros(1,length(yyyy));

%Loop for all Year
for i = yyyy

    counter = counter+1;

    [jd,~] = gre2jd(i,mm,dd,ut1,minute,second); %Running Julian Day

    P = precession(jd); %precession

    precCart = P * [x,y,z]';

    [r,alpha , delta] = CartToSphe(precCart(1),precCart(2),precCart(3));

    Calculated_delta(counter) = delta;
    Calculated_alpha(counter) = alpha;

end

alpha_ha = Calculated_alpha*12/pi; %To Hour angles
delta_deg = Calculated_delta/pi*180; %To degree

% Plot Task 1
figure % new figure
[hAx,~,~] = plotyy(yyyy,alpha_ha,yyyy,delta_deg);

title('Precession of a fictutious celestial obaject at (\alpha=0,\delta=0)')
xlabel('time (year)')

ylabel(hAx(1),'\alpha (h)') % left y-axis 
ylabel(hAx(2),'\delta (deg)') % right y-axis


%% Task 2 - Nutation

yyyy=1990:2030;
mm=1:12;
dd=1;
ut1=12;
minute=0;
second=0;

% coordinatesofa fictitious celestial object at the vernal equinox (0; 0) 
J2000_alpha = 0; 
J2000_delta = 0;

[x,y,z] = SphToCart(1,J2000_delta,J2000_alpha); 

counter = 0;

Calculated_delta = zeros(1,length(yyyy)*length(mm));
Calculated_alpha = zeros(1,length(yyyy)*length(mm));

%Loop for all Year
for i = yyyy
    for k = mm

        counter = counter+1;

        [jd,~] = gre2jd(i,k,dd,ut1,minute,second); %Running Julian Day

        N = nutation(jd); %precession

        NutationCart = N * [x,y,z]';

        [r,alpha , delta] = CartToSphe(NutationCart(1),NutationCart(2),NutationCart(3));

        Calculated_delta(counter) = delta;
        Calculated_alpha(counter) = alpha;

    end
end

alpha_ha = Calculated_alpha*12/pi; %To Hour angles
delta_deg = Calculated_delta/pi*180; %To degree

% Plot Task 2
YYYYMM = 1:(length(yyyy)*length(mm));

figure % new figure
[hAx,hLine1,hLine2] = plotyy(YYYYMM,alpha_ha,YYYYMM,delta_deg);

set(gca, 'XTick', 1:60:(length(yyyy)*length(mm))); % Change x-axis ticks
set(gca, 'XTickLabel', yyyy(1:5:length(yyyy))); % Change x-axis ticks labels.

title('Nutation of a fictutious celestial obaject at (\alpha=0,\delta=0)')
xlabel('time (year)')

ylabel(hAx(1),'\alpha (h)') % left y-axis 
ylabel(hAx(2),'\delta (deg)') % right y-axis

%% Task 3 - xxx

[jd,mjd] = gre2jd(2021,6,11,0,0,0); %Running Julian Day

yyyy=2021;
mm=3;
dd=6;
hour=2;
minute=13;
second=17;

polarMotion2021(yyyy,mm,dd,hour,minute,second)

