%% ref3d
% Discription:
% Creates a [3x3] reflection matrix
% usage:
% R_ref = ref3d(axis)
% input:
% axis  <= roation axis (1,2 or 3)
% output:
% R_i   <= [3x3] reflection matrix
% external calls:
% none

%% Function

function [R_ref] = ref3d(axis)

% Value checks:

if (axis ~= 1 && axis ~= 2 && axis ~= 3)
    error("2nd input argument 'axis' must be 1, 2 or 3.")
end

% computations

R_ref = eye(3);

switch axis
    case 1
        R_ref(1,1) = -1;
    case 2 
        R_ref(2,2) = -1;
    case 3
        R_ref(3,3) = -1;
end
end






