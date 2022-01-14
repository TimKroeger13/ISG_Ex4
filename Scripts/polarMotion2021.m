%% polarMotion2021
% Discription:
% Calculates the Polar motion in the year 2021.
% usage:
% [W] = polarMotion2021(yyyy,mm,dd,hour,minute,second)
% input:
% yyyy = Gregorian calender years
% mm = Gregorian calender month
% dd = Gregorian calender days
% hour = Gregorian calender hours
% minute = Gregorian calender minutes
% second = Gregorian calenderseconds
% output:
% W = polar motion matrix
% external calls:
% rot3d
% Author: Tim Kröger

function [W] = polarMotion2021(yyyy,mm,dd,hour,minute,second)

% Value checks:

if (~isnumeric(yyyy))
    error("yyyy is not numeric")
end
if (~isnumeric(mm))
    error("mm is not numeric")
end
if (~isnumeric(dd))
    error("dd is not numeric")
end
if (~isnumeric(hour))
    error("hour is not numeric")
end
if (~isnumeric(minute))
    error("minute is not numeric")
end
if (~isnumeric(second))
    error("second is not numeric")
end


if (mm<1 || mm>12)
    error("month must be a number between 1 and 12")
end
if (dd<1 || dd>31)
    error("Day must be a number between 1 and 31")
end
if (hour<0 || hour>23)
    error("Hour must be a number between 1 and 24")
end
if (minute<0 || minute>59)
    error("Minute must be a number between 1 and 60")
end
if (second<0 || second>59)
    error("Second must be a number between 1 and 60")
end
if(yyyy ~= 2021)
    error("currently only operating on year 2021")
end

if(mm==12 && dd > 12 && (minute>0 || seconds>0))
    error("Sorry but the Database is currently only based on data " + ...
        "till the 12.12.2021 at midnight")
end

% computations

IAU1980 = importdata('EOP_14_C04_IAU1980_yearly_files_2021_data.txt');
IAU1980_Names = ["Year","Month","Day","MJD","x","y","UT1-UTC",...
    "LOD","dPsi","dEps","x Err","y Err","UT1-UTC Err",...
    "LOD Err","dPsi Err","dEpsilon Err"];

%Get adjacent diurnal vaules

year_fit = IAU1980(:,1) == yyyy;
month_fit = IAU1980(:,2) == mm;
day_fit = IAU1980(:,3) == dd;

fit = year_fit & month_fit & day_fit;

Pre_date = IAU1980(fit,:);
Post_date = IAU1980(find(fit)+1,:);

%find the X coordiante
colum_x = IAU1980_Names == "x";
Pre_x = Pre_date(colum_x);
Post_x = Post_date(colum_x);

%find the Y coordiante
colum_y = IAU1980_Names == "y";
Pre_y = Pre_date(colum_y);
Post_y = Post_date(colum_y);

%Interploate

%Value to interploate to
% (when the Pre_date is 0 and the Post_date 1)

TimePercent = hour/24 + minute/24/60 + second/24/60/60;

x = interp1([0,Pre_x;1,Post_x],1+TimePercent,'linear'); %arcsec
y = interp1([0,Pre_y;1,Post_y],1+TimePercent,'linear'); %arcsec

%arcseconds to degree

x = x(2)/3600;
y = y(2)/3600;

W = rot3d(-x,2) * rot3d(-y,1);

end