%% SphToCart
% Discription:
% Converts spherical coordiantes to Catesian coordinates
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
% Author: Delaram + Tim Kröger

function [x,y,z] = SphToCart(r,phi,lambda)
x = r * sind(phi) * cosd(lambda);
y = r * sind(phi) * sind(lambda);
z = r * cosd(phi);
end

