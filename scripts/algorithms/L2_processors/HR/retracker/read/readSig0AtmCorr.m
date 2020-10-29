function data = readSig0AtmCorr(Sig0AtmCorr_path, filename_L1B, data)

% Finding the correct input file
[~,name,ext]    = fileparts(filename_L1B);
year            = name(17:20);
month           = name(21:22);
day             = name(23:24);
hour_aux        = str2double(name(26:27));

%let's define by the closest time 
hours_sampled = [0,6,12,18,24];
hours_sampled_string = {'00','06','12','18','24'};
[~,idx]=min(abs(hours_sampled-hour_aux));

%select the closest two maps in time
if hours_sampled(idx)==24
    %need to select next day
    hour_map_1=hours_sampled(idx-1);
    month
    day_map_1=day;
else
    
end


if hours_sampled(idx(1))==24
    hour = '00';
    day  = num2str(str2num(day)+1);
else
    hour = char(hours_sampled_string(idx(1)));
end


% if hour_aux < 6
%     hour = '00';
% elseif hour_aux < 12
%     hour = '06';
% elseif hour_aux < 18
%     hour = '12';
% else
%     hour = '18';
% end

filename_Sig0AtmCorr    = [Sig0AtmCorr_path,year,'\gfs_atten_ku_',year,month,day,'_',hour,'.nc'];


Sig0AtmCorr     = double(ncread(filename_Sig0AtmCorr,'sig0_atmos_corr'));
auxlat          = double(ncread(filename_Sig0AtmCorr,'lat')).';
auxlon          = double(ncread(filename_Sig0AtmCorr,'lon')).';

latsize         = length(auxlat);
lonsize         = length(auxlon);

[lat,lon]       = meshgrid(auxlat,auxlon);

% lat             = repmat(auxlat,[1,lonsize]);
% lon             = repmat(auxlon,[1,latsize]).';

F = scatteredInterpolant(lon(:),lat(:),Sig0AtmCorr(:));
data.COR.Sig0AtmCorr    = F(wrapTo360(data.GEO.LON),data.GEO.LAT);
%wrapping to 360 is required as definition of the longitude is - West and
%positive E from ISD and GPOD

% % Get the lat/lon/elev values inside our track
% ind_track       = inpolygon(lon,lat,[min(data.GEO.LON),max(data.GEO.LON)],[min(data.GEO.LAT),max(data.GEO.LAT)]);
% if max(max(ind_track)) || min(min(ind_track))
%     ind_dem     = inpolygon(data.GEO.LON,data.GEO.LAT,[min(lon(ind_track)),max(lon(ind_track))],[min(lat(ind_track)),max(lat(ind_track))]);
%     
%     
%     % Interpolatation
%     F = scatteredInterpolant(lon(ind_track),lat(ind_track),Sig0AtmCorr(ind_track));
%     data.COR.Sig0AtmCorr    = F(data.GEO.LON,data.GEO.LAT);
%         
% else %area is too small --> then we get the closest value for the whole track
%     [lon_dif_min,lon_ind_min] = min(abs(round(min(data.GEO.LON))-auxlon));
%     [lon_dif_max,lon_ind_max] = min(abs(round(max(data.GEO.LON))-auxlon));
%     if lon_ind_min ~= lon_ind_max
%         if lon_dif_min < lon_dif_max
%             lon_ind = lon_ind_min;
%         else
%             lon_ind = lon_ind_max;
%         end
%     else
%         lon_ind = lon_ind_min;
%     end
%     
%     [lat_dif_min,lat_ind_min] = min(abs(round(min(data.GEO.LAT))-auxlat));
%     [lat_dif_max,lat_ind_max] = min(abs(round(max(data.GEO.LAT))-auxlat));
%     if lat_ind_min ~= lat_ind_max
%         if lat_dif_min < lat_dif_max
%             lat_ind = lat_ind_min;
%         else
%             lat_ind = lat_ind_max;
%         end
%     else
%         lat_ind = lat_ind_min;
%     end
%     
%     data.COR.Sig0AtmCorr = Sig0AtmCorr(lat_ind,lon_ind) * ones(1,length(data.GEO.LAT));
%     
% end

end