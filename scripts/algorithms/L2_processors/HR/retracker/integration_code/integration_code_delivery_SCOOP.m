close all;

filename_L2_original = 'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\Phase_2\test_portion_code_SCOOP\L2\original_code_SCOOP\CR2_SR_2_WAT____20130102T134943_20130102T135354_20181013T171909_isd.nc';
%filename_L2_delivery = 'C:\Users\eduard.makhoul\Desktop\test_L2\L2\data\CR2_SR_2_WAT____20130102T134943_20130102T135354_20190130T110227_isd.nc';
filename_L2_delivery ='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\Phase_2\test_portion_code_SCOOP\L2\only_portion_code_SCOOP\from_code_DL380\CR2_SR_2_WAT____20130102T134943_20130102T135354_20190130T153336_isd.nc';
resultPath = 'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\Phase_2\test_portion_code_SCOOP\L2\integration_comparison\';
mkdir(resultPath);
TRP=0;

nccmp(filename_L2_original,filename_L2_delivery);

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


%% ---------------------- read L2 delivery -------------------------------------
data_delivery.GEO.TAI        = double(ncread(filename_L2_delivery,'time_20_ku').');
data_delivery.GEO.LAT        = double(ncread(filename_L2_delivery,'lat_20_ku').');
data_delivery.GEO.LON        = double(ncread(filename_L2_delivery,'lon_20_ku').');
data_delivery.GEO.H          = double(ncread(filename_L2_delivery,'alt_20_ku').');
data_delivery.MEA.win_delay  = double(ncread(filename_L2_delivery,'retracked_range_analytical_SWH_MSSfixed_20_ku').').*2/c_cst;
data_delivery.MEA.ssh        = double(ncread(filename_L2_delivery,'ssh_analytical_SWH_MSSfixed_20_ku').');
data_delivery.MEA.swh        = double(ncread(filename_L2_delivery,'swh_analytical_SWH_MSSfixed_20_ku').');
data_delivery.MEA.sigma0     = double(ncread(filename_L2_delivery,'sig0_analytical_SWH_MSSfixed_20_ku').');
data_delivery.MEA.corr       = double(ncread(filename_L2_delivery,'Pearson_corr_analytical_SWH_MSSfixed_20_ku').');
data_delivery.MEA.Flag_validity       = double(ncread(filename_L2_delivery,'Flag_validity_L1B_wvfm_20_ku').');


%% ---------------------- read L2 original -------------------------------------
data_original.GEO.TAI        = double(ncread(filename_L2_original,'time_20_ku').');
data_original.GEO.LAT        = double(ncread(filename_L2_original,'lat_20_ku').');
data_original.GEO.LON        = double(ncread(filename_L2_original,'lon_20_ku').');
data_original.GEO.H          = double(ncread(filename_L2_original,'alt_20_ku').');
data_original.MEA.win_delay  = double(ncread(filename_L2_original,'retracked_range_analytical_SWH_MSSfixed_20_ku').').*2/c_cst;
data_original.MEA.ssh        = double(ncread(filename_L2_original,'ssh_analytical_SWH_MSSfixed_20_ku').');
data_original.MEA.swh        = double(ncread(filename_L2_original,'swh_analytical_SWH_MSSfixed_20_ku').');
data_original.MEA.sigma0     = double(ncread(filename_L2_original,'sig0_analytical_SWH_MSSfixed_20_ku').');
data_original.MEA.corr       = double(ncread(filename_L2_original,'Pearson_corr_analytical_SWH_MSSfixed_20_ku').');
data_original.MEA.Flag_validity       = double(ncread(filename_L2_original,'Flag_validity_L1B_wvfm_20_ku').');

N_surfaces=length(data_original.GEO.TAI);



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
plot(1:N_surfaces,data_delivery.GEO.TAI(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.GEO.TAI,'-r'); 
legend('delivery','original')
title('time_tai_20_ku','Interpreter','none');
ylabel('[s]'); xlabel('Surface'); 
subplot(2,1,2)
[hAX]=plotyy(1:N_surfaces,(data_original.GEO.TAI-data_delivery.GEO.TAI(idx_int)),...
       1:N_surfaces,abs(data_original.GEO.TAI-data_delivery.GEO.TAI(idx_int))./(data_delivery.GEO.TAI(idx_int))*100);
title('Error time_tai_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[s]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'time_tai_20_ku_comparison',file_ext)); 

%latitude_20_ku
figure('Name','latitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.GEO.LAT(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.GEO.LAT,'-r'); 
legend('delivery','original')
title('latitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.GEO.LAT-data_delivery.GEO.LAT(idx_int)),...
       1:N_surfaces,abs((data_original.GEO.LAT-data_delivery.GEO.LAT(idx_int))./(data_delivery.GEO.LAT(idx_int)))*100);
title('Error latitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'latitude_20_ku_comparison',file_ext)); 


%longitude_20_ku
figure('Name','longitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.GEO.LON(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.GEO.LON,'-r'); 
legend('delivery','original')
title('longitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.GEO.LON-data_delivery.GEO.LON(idx_int)),...
       1:N_surfaces,abs((data_original.GEO.LON-data_delivery.GEO.LON(idx_int))./(data_delivery.GEO.LON(idx_int)))*100);
title('Error longitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'longitude_20_ku_comparison',file_ext)); 

%com_altitude_20_ku
figure('Name','alt_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.GEO.H(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.GEO.H,'-r'); 
legend('delivery','original')
title('alt_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.GEO.H-data_delivery.GEO.H(idx_int)),...
       1:N_surfaces,abs((data_original.GEO.H-data_delivery.GEO.H(idx_int))./(data_delivery.GEO.H(idx_int)))*100);
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
plot(1:N_surfaces,data_delivery.MEA.win_delay(idx_int).*c_cst/2,'-b'); 
hold on; plot(1:N_surfaces,data_original.MEA.win_delay.*c_cst/2,'-r'); 
legend('delivery','original')
title('range_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.MEA.win_delay-data_delivery.MEA.win_delay(idx_int)).*c_cst/2,...
       1:N_surfaces,abs((data_original.MEA.win_delay-data_delivery.MEA.win_delay(idx_int))./(data_delivery.MEA.win_delay(idx_int)))*100);
title('Error range_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'range_20_ku_comparison',file_ext));

%% --------------- GEOPHYSICAL RETRIEVALS ---------------------------------
% --------------- SSH -----------------------------------------------------
figure('Name','ssh_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.MEA.ssh(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.MEA.ssh,'-r'); 
legend('delivery','original')
title('ssh_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.MEA.ssh-data_delivery.MEA.ssh(idx_int)),...
       1:N_surfaces,abs((data_original.MEA.ssh-data_delivery.MEA.ssh(idx_int))./(data_delivery.MEA.ssh(idx_int)))*100);
title('Error ssh_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'ssh_20_ku_comparison',file_ext));

% --------------- SWH -----------------------------------------------------
figure('Name','swh_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.MEA.swh(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.MEA.swh,'-r'); 
legend('delivery','original')
title('swh_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.MEA.swh-data_delivery.MEA.swh(idx_int)),...
       1:N_surfaces,abs((data_original.MEA.swh-data_delivery.MEA.swh(idx_int))./(data_delivery.MEA.swh(idx_int)))*100);
title('Error swh_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'swh_20_ku_comparison',file_ext));

% --------------- sigma0 --------------------------------------------------
figure('Name','sigma0_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.MEA.sigma0(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.MEA.sigma0,'-r'); 
legend('delivery','original')
title('sigma0_20_ku','Interpreter','none');
ylabel('[dB]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.MEA.sigma0-data_delivery.MEA.sigma0(idx_int)),...
       1:N_surfaces,abs((data_original.MEA.sigma0-data_delivery.MEA.sigma0(idx_int))./(data_delivery.MEA.sigma0(idx_int)))*100);
title('Error sigma0_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[dB]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'sigma0_20_ku_comparison',file_ext));

% --------------- Pearson correlation coefficient -------------------------
figure('Name','fitting_cor_coef_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_delivery.MEA.corr(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_original.MEA.corr,'-r'); 
legend('delivery','original')
title('fitting_cor_coef_20_ku','Interpreter','none');
ylabel('[%]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_original.MEA.corr-data_delivery.MEA.corr(idx_int)),...
       1:N_surfaces,abs((data_original.MEA.corr-data_delivery.MEA.corr(idx_int))./(data_delivery.MEA.corr(idx_int)))*100);
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
% title('GR (delivery)');
% 
% subplot(2,2,2)
% imagesc(10*log10(xxall)); colormap('jet'); c=colorbar;
% ylabel(c,'[dB]'); caxis([min_image,max_image])
% title('GR (C++)');
% 
% subplot(2,2,3)
% imagesc(abs(stack_real-xxall)); colormap('jet'); c=colorbar;
% xlabel('Samples'); ylabel('Beams'); 
% title('Absolute error (delivery - C++)');
% 
% subplot(2,2,4)
% imagesc(abs(stack_real-xxall)./(stack_real)); colormap('jet'); c=colorbar;
% ylabel(c,'[%]')
% caxis([0,100]);
% xlabel('Samples'); ylabel('Beams'); 
% title('Relative error (delivery - C++)');
