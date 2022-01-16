%% F
% Discription:
% Get the rest from the fix function.
% usage:
% Fx=F(x)
% input:
% x = number (real)
% output:
% Fx = number (real) 0>= Fx <= 1
% external calls:
% none

function Fx=F(x)

if (~isnumeric(x))
    error("x is not numeric")
end

Fx= x -fix(x);
end