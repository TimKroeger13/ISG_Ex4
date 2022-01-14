%% nutation
% Discription:
% Calculates the nutation based on J2000.0.
% usage:
% [N] = nutation(jd)
% input:
% jd = Julian date
% output:
% N = nutation Matrix
% external calls:
% rot3d
% Expand_DecDeg
% De_Mi_Se_To_DecDeg
% Author: Tim Kröger

function [N] = nutation(jd)

% Value checks:
if (~isnumeric(jd))
    error("jd is not numeric")
end

% computations

d = jd -2451545.0; %Based on year 2000
T = d/36525; %2000 based Julian centurie


epsilon_A_v1 = De_Mi_Se_To_DecDeg(Expand_DecDeg(84381.448));
epsilon_A_v2 = De_Mi_Se_To_DecDeg(Expand_DecDeg(46.8150));

delta_Psi_v1 = 0.0048;
delta_Psi_v2 = 0.0004;

delta_epsilon_v1 = 0.0026;
delta_epsilon_v2 = 0.0002;

f1_v1 = 125;
f1_v2 = 0.05295;

f2_v1 = 200.9;
f2_v2 = 1.97129;

f1 = f1_v1 - f1_v2 * d;
f2 = f2_v1 + f2_v2 * d;


epsilon_A = epsilon_A_v1 - epsilon_A_v2 * T;

delta_Psi = -delta_Psi_v1 * sind(f1) - delta_Psi_v2 * sind(f2);

delta_epsilon = delta_epsilon_v1 * cosd(f1) - delta_epsilon_v2 * cosd(f2);

N = rot3d(-epsilon_A-delta_epsilon,1) *  rot3d(-delta_Psi,3) *  rot3d(epsilon_A,1);

end