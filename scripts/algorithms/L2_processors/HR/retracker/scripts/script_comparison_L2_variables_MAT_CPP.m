close all;
common_path = '\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Integration_test\L2\';
filename_L2_CPP = strcat(common_path,'inputs\S6A_OS21_P4__HR__GR_20170305T070423_20170305T070628_0001.nc');
filename_L2_MAT = strcat(common_path,'inputs\S6A_OS21_P4__HR__L2_20170305T070423_20170305T070628_0001_isd.nc');
resultPath = strcat(common_path,'results\');
mkdir(resultPath);
TRP=0;

global c_cst sec_in_day_cst
c_cst=299792458;
sec_in_day_cst=86400;

%configuration
cnf_p.retracker_name='ANALYTICAL';
cnf_p.looks_index_method='Look_angle';
cnf_p.look_ang_method='approximate';
cnf_p.mode = 'RAW';

set_default_plot;
set(0,'defaultFigureVisible','on');

figure_format='jpg';
res_fig='-r300';

switch lower(figure_format)
    case 'pdf'
        file_ext='.pdf';
        print_file='-dpdf';
    case 'eps'
        file_ext='.eps';
        print_file='-depsc';
    case 'png'
        file_ext='.png';
        print_file='-dpng';
    case 'jpg'
        file_ext='.jpg';
        print_file='-djpeg';        
end


%% ---------------------- read L2 MAT -------------------------------------
data_MAT.GEO.TAI        = double(ncread(filename_L2_MAT,'time_20_ku').');
data_MAT.GEO.LAT        = double(ncread(filename_L2_MAT,'lat_20_ku').');
data_MAT.GEO.LON        = double(ncread(filename_L2_MAT,'lon_20_ku').');
data_MAT.GEO.H          = double(ncread(filename_L2_MAT,'alt_20_ku').');
data_MAT.MEA.win_delay  = double(ncread(filename_L2_MAT,'retracked_range_analytical_SWH_MSSfixed_20_ku').').*2/c_cst;
data_MAT.MEA.ssh        = double(ncread(filename_L2_MAT,'ssh_analytical_SWH_MSSfixed_20_ku').');
data_MAT.MEA.swh        = double(ncread(filename_L2_MAT,'swh_analytical_SWH_MSSfixed_20_ku').');
data_MAT.MEA.sigma0     = double(ncread(filename_L2_MAT,'sig0_analytical_SWH_MSSfixed_20_ku').');
data_MAT.MEA.corr       = double(ncread(filename_L2_MAT,'Pearson_corr_analytical_SWH_MSSfixed_20_ku').');


%% ---------------------- read L2 CPP -------------------------------------
data_CPP.GEO.TAI        = double(ncread(filename_L2_CPP,'data_20/ku/time_tai').');
data_CPP.GEO.LAT        = double(ncread(filename_L2_CPP,'data_20/ku/latitude').');
data_CPP.GEO.LON        = double(wrapTo180(ncread(filename_L2_CPP,'data_20/ku/longitude').'));
data_CPP.GEO.H          = double(ncread(filename_L2_CPP,'data_20/ku/com_altitude').');
data_CPP.MEA.win_delay  = double(ncread(filename_L2_CPP,'data_20/ku/range').').*2/c_cst;
data_CPP.MEA.ssh        = double(ncread(filename_L2_CPP,'data_20/ku/ssh').');
data_CPP.MEA.swh        = double(ncread(filename_L2_CPP,'data_20/ku/hs').');
data_CPP.MEA.sigma0     = double(ncread(filename_L2_CPP,'data_20/ku/sigma0').');
data_CPP.MEA.corr       = double(ncread(filename_L2_CPP,'data_20/ku/fitting_cor_coef').')*100.0;

N_surfaces=length(data_CPP.GEO.TAI);



if TRP
    idx_int=[28,122,216,309];
else
    idx_int=1:N_surfaces;
end

set_default_plot;
set(0,'defaultFigureVisible','on');

%% ------------ GEO: geographical information -----------------------------

%time_tai_20_ku
figure('Name','time_tai_20_ku');
subplot(2,1,1)
plot(1:N_surfaces,data_MAT.GEO.TAI(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.GEO.TAI,'-r'); 
legend('MAT','CPP')
title('time_tai_20_ku','Interpreter','none');
ylabel('[s]'); xlabel('Surface'); 
subplot(2,1,2)
[hAX]=plotyy(1:N_surfaces,(data_CPP.GEO.TAI-data_MAT.GEO.TAI(idx_int)),...
       1:N_surfaces,abs(data_CPP.GEO.TAI-data_MAT.GEO.TAI(idx_int))./(data_MAT.GEO.TAI(idx_int))*100);
title('Error time_tai_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[s]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'time_tai_20_ku_comparison',file_ext)); 

%latitude_20_ku
figure('Name','latitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.LAT(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.GEO.LAT,'-r'); 
legend('MAT','CPP')
title('latitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.GEO.LAT-data_MAT.GEO.LAT(idx_int)),...
       1:N_surfaces,abs((data_CPP.GEO.LAT-data_MAT.GEO.LAT(idx_int))./(data_MAT.GEO.LAT(idx_int)))*100);
title('Error latitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'latitude_20_ku_comparison',file_ext)); 


%longitude_20_ku
figure('Name','longitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.LON(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.GEO.LON,'-r'); 
legend('MAT','CPP')
title('longitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.GEO.LON-data_MAT.GEO.LON(idx_int)),...
       1:N_surfaces,abs((data_CPP.GEO.LON-data_MAT.GEO.LON(idx_int))./(data_MAT.GEO.LON(idx_int)))*100);
title('Error longitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'longitude_20_ku_comparison',file_ext)); 

%com_altitude_20_ku
figure('Name','alt_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.H(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.GEO.H,'-r'); 
legend('MAT','CPP')
title('alt_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.GEO.H-data_MAT.GEO.H(idx_int)),...
       1:N_surfaces,abs((data_CPP.GEO.H-data_MAT.GEO.H(idx_int))./(data_MAT.GEO.H(idx_int)))*100);
title('Error alt_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'com_altitude_20_ku_comparison',file_ext)); 

%close all;

%% ------------ MEA: measurements -----------------------------------------
%tracker_range_calibrated_20_ku
figure('Name','range_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.win_delay(idx_int).*c_cst/2,'-b'); 
hold on; plot(1:N_surfaces,data_CPP.MEA.win_delay.*c_cst/2,'-r'); 
legend('MAT','CPP')
title('range_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.MEA.win_delay-data_MAT.MEA.win_delay(idx_int)).*c_cst/2,...
       1:N_surfaces,abs((data_CPP.MEA.win_delay-data_MAT.MEA.win_delay(idx_int))./(data_MAT.MEA.win_delay(idx_int)))*100);
title('Error range_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'range_20_ku_comparison',file_ext));

%% --------------- GEOPHYSICAL RETRIEVALS ---------------------------------
% --------------- SSH -----------------------------------------------------
figure('Name','ssh_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.ssh(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.MEA.ssh,'-r'); 
legend('MAT','CPP')
title('ssh_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.MEA.ssh-data_MAT.MEA.ssh(idx_int)),...
       1:N_surfaces,abs((data_CPP.MEA.ssh-data_MAT.MEA.ssh(idx_int))./(data_MAT.MEA.ssh(idx_int)))*100);
title('Error ssh_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'ssh_20_ku_comparison',file_ext));

% --------------- SWH -----------------------------------------------------
figure('Name','swh_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.swh(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.MEA.swh,'-r'); 
legend('MAT','CPP')
title('swh_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.MEA.swh-data_MAT.MEA.swh(idx_int)),...
       1:N_surfaces,abs((data_CPP.MEA.swh-data_MAT.MEA.swh(idx_int))./(data_MAT.MEA.swh(idx_int)))*100);
title('Error swh_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'swh_20_ku_comparison',file_ext));

% --------------- sigma0 --------------------------------------------------
figure('Name','sigma0_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.sigma0(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.MEA.sigma0,'-r'); 
legend('MAT','CPP')
title('sigma0_20_ku','Interpreter','none');
ylabel('[dB]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.MEA.sigma0-data_MAT.MEA.sigma0(idx_int)),...
       1:N_surfaces,abs((data_CPP.MEA.sigma0-data_MAT.MEA.sigma0(idx_int))./(data_MAT.MEA.sigma0(idx_int)))*100);
title('Error sigma0_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[dB]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'sigma0_20_ku_comparison',file_ext));

% --------------- Pearson correlation coefficient -------------------------
figure('Name','fitting_cor_coef_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.corr(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_CPP.MEA.corr,'-r'); 
legend('MAT','CPP')
title('fitting_cor_coef_20_ku','Interpreter','none');
ylabel('[%]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_CPP.MEA.corr-data_MAT.MEA.corr(idx_int)),...
       1:N_surfaces,abs((data_CPP.MEA.corr-data_MAT.MEA.corr(idx_int))./(data_MAT.MEA.corr(idx_int)))*100);
title('Error fitting_cor_coef_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[%]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'fitting_cor_coef_20_ku_comparison',file_ext));

% max_image=10*log10(max([max(stack_real(:)),max(xxall(:))]));
% xxall(stack_real==0)=0;
% min_image=max_image-40.0;
% 
% figure;
% subplot(2,2,1)
% imagesc(10*log10(stack_real)); colormap('jet'); c=colorbar;
% ylabel(c,'[dB]'); caxis([min_image,max_image])
% xlabel('Samples'); ylabel('Beams'); 
% title('GR (MAT)');
% 
% subplot(2,2,2)
% imagesc(10*log10(xxall)); colormap('jet'); c=colorbar;
% ylabel(c,'[dB]'); caxis([min_image,max_image])
% title('GR (C++)');
% 
% subplot(2,2,3)
% imagesc(abs(stack_real-xxall)); colormap('jet'); c=colorbar;
% xlabel('Samples'); ylabel('Beams'); 
% title('Absolute error (MAT - C++)');
% 
% subplot(2,2,4)
% imagesc(abs(stack_real-xxall)./(stack_real)); colormap('jet'); c=colorbar;
% ylabel(c,'[%]')
% caxis([0,100]);
% xlabel('Samples'); ylabel('Beams'); 
% title('Relative error (MAT - C++)');
