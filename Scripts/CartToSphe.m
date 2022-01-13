%% CartToSphe
% Discription:
% Converts Catesian coordiantes to spherical coordinates
% usage:
% [r,lamda , phi] = CartToSphe(x,y,z)
% input:
% x = The x coordinate in the cartesian sythem
% y = The x coordinate in the cartesian sythem
% z = The x coordinate in the cartesian sythem
% output:
% r <= Radus of the sphere.
% phi <= degree phi.
% lambda <= degree lambda.
% external calls:
% none
% Author: Delaram + Tim Kröger

function [r,lamda , phi] = CartToSphe(x,y,z)
r = sqrt(x^2 + y^2 + z^2);
lamda = atan2d(y,x);
phi = atan2d(sqrt(x^2+y^2),z);
end




