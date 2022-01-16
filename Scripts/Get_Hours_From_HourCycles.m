%% Get_Hours_From_HourCycles
% Discription:
% Get the Hour  from multiple cycles of Hours.
% usage:
% [hour] = Get_Hours_From_HourCycles(hour)
% input:
% hour = hour
% output:
% hour = The hour where 0>= hour <= 24
% external calls:
% none

function [hour] = Get_Hours_From_HourCycles(hour)

% Value checks:

if (~isnumeric(hour))
    error("hour is not numeric")
end

%calculations

while hour>=24

    hour = hour - 24;

end
while hour<0

    hour = hour + 24;

end
end


