%% SphToCart
% Discription:
% Converts Catesian coordiantes to spherical coordinates
% everything should be in radian
% usage:
% [x,y,z] = SphToCart(r,phi,lambda)
% input:
% r <= Radus of the sphere.
% phi <= degree phi.
% lambda <= degree lambda.
% output:
% x = The x coordinate in the cartesian sythem
% y = The x coordinate in the cartesian sythem
% z = The x coordinate in the cartesian sythem
% external calls:
% none
% Author: Delaram's Group (Tempname)

function [x,y,z] = SphToCart(r,phi,lambda)
x=r*cos(lambda)*cos(phi);
y= r*sin(lambda)*cos(phi);
z=r*sin(phi);
end

