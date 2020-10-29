%% Comparison retrievals SWH DeDop and L2 IPF

file_DeDop_bs1 = 'C:\Users\eduard.makhoul\isardSAT\projects\CCI_SEA_STATE\data\dedop-processed\L2\S3\data\L2_measurement_l1a_cropped_s3-baselineWAT.nc';
file_DeDop_bs2 = 'C:\Users\eduard.makhoul\isardSAT\projects\CCI_SEA_STATE\data\dedop-processed\L2\S3_zp2_Hamm\data\L2_measurement_l1a_cropped_hammingWAT.nc';
file_L2_IPF = 'C:\Users\eduard.makhoul\isardSAT\projects\CCI_SEA_STATE\data\S3A_SR_2_LAN____20171005T210218_20171005T215242_20180520T144733_3024_023_071______LR1_R_NT_003.SEN3\standard_measurement_cropped.nc';

filter_outlier = 0;
filter_latitude = 1;

SWH_DeDop_bs1=ncread(file_DeDop_bs1,'swh_analytical_SWH_MSSfixed_20_ku');
lat_DeDop_bs1=ncread(file_DeDop_bs1,'lat_20_ku');
lon_DeDop_bs1=ncread(file_DeDop_bs1,'lon_20_ku');
SWH_DeDop_bs2=ncread(file_DeDop_bs2,'swh_analytical_SWH_MSSfixed_20_ku');
lat_DeDop_bs2=ncread(file_DeDop_bs2,'lat_20_ku');
SWH_L2_IPF = ncread(file_L2_IPF,'swh_ocean_20_ku');
lat_L2_IPF = ncread(file_L2_IPF,'lat_20_ku');
lon_L2_IPF = ncread(file_L2_IPF,'lon_20_ku');

lla2kmlWaveforms_noimage('track_DeDop_CCI_SS','DeDop',lat_DeDop_bs1,lon_DeDop_bs1,zeros(1,length(lat_DeDop_bs1)), '.');
lla2kmlWaveforms_noimage('track_IPF_CCI_SS','IPF',lat_L2_IPF,lon_L2_IPF,zeros(1,length(lat_L2_IPF)), '.');


%% Filtering by latitude
if filter_latitude==1
lat_min=48.10;
lat_max=51.70;

idx_lat =find(lat_DeDop_bs1<lat_min | lat_DeDop_bs1>lat_max);
SWH_DeDop_bs1(idx_lat)=NaN;
clear idx_lat;

idx_lat =find(lat_DeDop_bs2<lat_min | lat_DeDop_bs2>lat_max);
SWH_DeDop_bs2(idx_lat)=NaN;
clear idx_lat;

idx_lat =find(lat_L2_IPF<lat_min | lat_L2_IPF>lat_max);
SWH_L2_IPF(idx_lat)=NaN;
clear idx_lat;
end

%% Outliers filtering
if filter_outlier==1
hampel_wind=3;% size of window half size
hampel_sigma=3; %number of std deviations to which a sample differ from local median

[~,idx_outliers] = hampel(lat_DeDop_bs1,SWH_DeDop_bs1,hampel_wind,hampel_sigma);
SWH_DeDop_bs1(idx_outliers)=NaN;
clear idx_outliers;

[~,idx_outliers] = hampel(lat_DeDop_bs2,SWH_DeDop_bs2,hampel_wind,hampel_sigma);
SWH_DeDop_bs2(idx_outliers)=NaN;
clear idx_outliers;

[~,idx_outliers] = hampel(lat_L2_IPF,SWH_L2_IPF,hampel_wind,hampel_sigma);
SWH_L2_IPF(idx_outliers)=NaN;
clear idx_outliers;
end

%% Statistics
win_size_stat = 20;
[~,~,std_SWH_DeDop_bs1,~]=block_statistics(SWH_DeDop_bs1,win_size_stat);
[~,~,std_SWH_DeDop_bs2,~]=block_statistics(SWH_DeDop_bs2,win_size_stat);
[~,~,std_SWH_L2_IPF,~]=block_statistics(SWH_L2_IPF,win_size_stat);

%% Ploting

figure; 
plot(lat_L2_IPF,SWH_L2_IPF,'*m');
hold on;
plot(lat_DeDop_bs1,SWH_DeDop_bs1,'+b');
plot(lat_DeDop_bs2,SWH_DeDop_bs2,'sg');
ylabel('SWH [m]');
xlabel('Latitude [deg.]');
title(strcat({'std [cm]: ESA IPF = '},num2str(std_SWH_L2_IPF.*100),{', DeDop S-3 = '},num2str(std_SWH_DeDop_bs1*100),{', DeDop zp2 Hamm = '},num2str(std_SWH_DeDop_bs2*100)))
legend('ESA IPF','DeDop S-3','DeDop zp2 Hamm');