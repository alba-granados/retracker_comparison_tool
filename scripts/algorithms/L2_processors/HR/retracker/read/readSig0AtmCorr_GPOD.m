function data = readSig0AtmCorr_GPOD(Sig0AtmCorr_path, filename_L1B, data)

% Finding the correct input file
[~,name,ext]    = fileparts(filename_L1B);
%initial time of acquisition
year            = name(17:20);
month           = name(21:22);
day             = name(23:24);
hour_aux        = str2double(name(26:27));

%list all the files and from there take the ones with the day of initial
%acquisition: help in taking the next file in case hour of interest is
%larger than 18h (should take next day--> change of month or even year)
files=dir(strcat(Sig0AtmCorr_path,filesep,'*.nc'));
aux=struct2cell(files); aux=aux(1,:); %Keep the
true_false_idx=(~cellfun(@isempty,strfind(aux,strcat(year,month,day))));
idx_int =find(true_false_idx);
idx_int = [idx_int,idx_int(end)+1];

%find the two closest maps in time (where it lies in between)
%let's define by the closest time 
hours_sampled = [0,6,12,18,24];
d=sort(abs(hours_sampled-hour_aux));
if (d(1) == d(2))
      vals = find(abs(hours_sampled-hour_aux)==d(1));
      lowest = vals(1);
      second_lowest = vals(2);
else
      lowest=find(abs(hours_sampled-hour_aux)==d(1));
      second_lowest=find(abs(hours_sampled-hour_aux)==d(2));
end
idx_selec = sort([lowest,second_lowest]);

%% Interpolation with lat/lon for the two closest maps
% -------------------- Map left time -----------------------------------
filename_Sig0AtmCorr_Map_left_time    = [Sig0AtmCorr_path,char(aux(idx_int(idx_selec(1))))];

Sig0AtmCorr     = double(ncread(filename_Sig0AtmCorr_Map_left_time,'sig0_atmos_corr'));
auxlat          = double(ncread(filename_Sig0AtmCorr_Map_left_time,'lat')).';
auxlon          = double(ncread(filename_Sig0AtmCorr_Map_left_time,'lon')).';

[lat,lon]       = meshgrid(auxlat,auxlon);

F = scatteredInterpolant(lon(:),lat(:),Sig0AtmCorr(:));
Sig0AtmCorr_left_time    = F(wrapTo360(data.GEO.LON),data.GEO.LAT);


% -------------------- Map right time -----------------------------------
filename_Sig0AtmCorr_Map_right_time    = [Sig0AtmCorr_path,char(aux(idx_int(idx_selec(2))))];

Sig0AtmCorr     = double(ncread(filename_Sig0AtmCorr_Map_right_time,'sig0_atmos_corr'));
auxlat          = double(ncread(filename_Sig0AtmCorr_Map_right_time,'lat')).';
auxlon          = double(ncread(filename_Sig0AtmCorr_Map_right_time,'lon')).';

[lat,lon]       = meshgrid(auxlat,auxlon);

F = scatteredInterpolant(lon(:),lat(:),Sig0AtmCorr(:));
Sig0AtmCorr_right_time    = F(wrapTo360(data.GEO.LON),data.GEO.LAT);


%% Interpolation in time between the two maps
% Convert the file time of the two maps to UTC seconds since 2000-01-01 00:00:00
% From L1B we have UTC per surface
% dumm = char(aux(idx_int(idx_selec(1))));
%use the UTC time in data.GEO.TAI.total (tot i dir TAI en ppi es UTC
%seconds, check) to perform an interpolation between the data from the two
%associated times
% UTC_sec_left_time = convert_UTC_sec_Jan_2000(str2num(dumm(14:17)),str2num(dumm(18:19)),...
%                                              str2num(dumm(20:21)),str2num(dumm(23:24)),0,0);
% dumm=char(aux(idx_int(idx_selec(2))));
%UTC_sec_right_time = convert_UTC_sec_Jan_2000(str2num(dumm(14:17)),str2num(dumm(18:19)),...
%                                              str2num(dumm(20:21)),str2num(dumm(23:24)),0,0);

% Define number of seconds since UTC standard zero time to UTC2000
t0 = posixtime(datetime(2000, 1, 1, 0, 0, 0));

dumm_l = char(aux(idx_int(idx_selec(1))));
year_l  = dumm_l(14:17);
month_l = dumm_l(18:19);
day_l   = dumm_l(20:21);
hour_l  = dumm_l(23:24);
UTC2000_sec_left_time_2 = posixtime(datetime(str2double(year_l), str2double(month_l), str2double(day_l), str2double(hour_l), 0, 0)) - t0;

dumm_r=char(aux(idx_int(idx_selec(2))));
year_r  = dumm_r(14:17);
month_r = dumm_r(18:19);
day_r   = dumm_r(20:21);
hour_r  = dumm_r(23:24);                                           
UTC2000_sec_right_time_2 = posixtime(datetime(str2double(year_r), str2double(month_r), str2double(day_r), str2double(hour_r), 0, 0)) - t0;
                      
% Interpolate
% (actually, solve a triangle "regla de tres" segons:
%   x_i = x_A + ((t_i - t_A)/(t_B - t_B))*(x_B - x_A)
int_sig0val = Sig0AtmCorr_left_time + ...
    ((data.GEO.TAI.total - UTC2000_sec_left_time_2)./(UTC2000_sec_right_time_2 - UTC2000_sec_left_time_2)).*(Sig0AtmCorr_right_time - Sig0AtmCorr_left_time);

% Store in output object along with input data structure content:
data.COR_sig0.Sig0AtmCorr = int_sig0val;

end




