function [UTC_sec_Jan_2000]=convert_UTC_sec_Jan_2000(year,month,day,HH,MM,SS)
sec_in_day=24*60*60;
days_months=[31,28,31,30,31,30,31,31,30,31,30,31];

%compute the leap days
years_from_2000 = (2000:1:year);
leap_days = 0;
for i_year=1:length(years_from_2000)
    if mod(years_from_2000(i_year),4)~=0
    elseif mod(years_from_2000(i_year),100)~=0
        leap_days=leap_days+1;
    elseif mod(years_from_2000(i_year),400)~=0
    else
        leap_days=leap_days+1;
    end
end

% total number of seconds since 1st January 2000
if month > 1
    UTC_sec_Jan_2000 = (year-2000).*365*sec_in_day+(sum(days_months(1:month-1)))*sec_in_day...
        +(day-1)*sec_in_day+HH*60*60+MM*60+SS;
else
    UTC_sec_Jan_2000 = (year-2000).*365*sec_in_day+...
        +(day-1)*sec_in_day+HH*60*60+MM*60+SS;
end

%add leap days
if month>2
    UTC_sec_Jan_2000 = UTC_sec_Jan_2000 + leap_days*sec_in_day;
else
    UTC_sec_Jan_2000 = UTC_sec_Jan_2000 + (leap_days-1)*sec_in_day;
end

end