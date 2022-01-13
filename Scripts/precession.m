%% precession
% Discription:
% Calculates the precession based on Newcomb parameters.
% usage:
% [p] = precession(jd)
% input:
% jd = Julian date
% output:
% p = precession Matrix
% external calls:
% rot3d
% Expand_DecDeg
% De_Mi_Se_To_DecDeg
% Author: Tim Kröger

function [p] = precession(jd)

% Value checks:
if (~isnumeric(jd))
    error("jd is not numeric")
end

% computations

d = jd -2451545.0; %Based on year 2000
T = d/36525; %2000 based Julian centurie

%Zeta values
zeta_v1 = De_Mi_Se_To_DecDeg(Expand_DecDeg(2306.2181));
zeta_v2 = De_Mi_Se_To_DecDeg(Expand_DecDeg(0.30188));

%Upper Zeta values
UpperZeta_v1 = De_Mi_Se_To_DecDeg(Expand_DecDeg(2306.2181));
UpperZeta_v2 = De_Mi_Se_To_DecDeg(Expand_DecDeg(1.09468));

%teta Zeta values
teta_v1 = De_Mi_Se_To_DecDeg(Expand_DecDeg(2004.3109));
teta_v2 = De_Mi_Se_To_DecDeg(Expand_DecDeg(0.42665));

zeta= zeta_v1 * T + zeta_v2 * T^2;
UpperZeta= UpperZeta_v1 * T + UpperZeta_v2 * T^2;
teta= teta_v1 * T + teta_v2 * T^2;

p = rot3d(-UpperZeta,3) * rot3d(teta,2) * rot3d(-zeta,3);

end


