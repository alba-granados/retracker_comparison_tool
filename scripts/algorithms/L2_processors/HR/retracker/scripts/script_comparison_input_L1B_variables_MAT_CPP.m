close all;
common_path = '\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Integration_test\L1B\';
filename_L1B_MAT = strcat(common_path,'inputs\S6A_OS20_P4__RAW_1B_00000000T000000_99999999T999999_0001.mat');
filename_L1B_GPP = strcat(common_path,'inputs\S6A_OS20_P4__HR__1B_20170305T065123_20170305T065233_0001.nc');
resultPath       = strcat(common_path,'results\');
cnf_p.mode = 'RAW';

mkdir(resultPath);

%load('.\plotting\colormapATDD.mat');
colormap_set = 'hot';%colormapATDD;

TRP=0;

global c_cst sec_in_day_cst flat_coeff_cst semi_major_axis_cst
flat_coeff_cst = 0.003352810664747;
c_cst=299792458;
sec_in_day_cst=86400;
semi_major_axis_cst = 6378137.0;

%configuration
cnf_p.retracker_name='ANALYTICAL';
cnf_p.looks_index_method='Look_angle';
cnf_p.look_ang_method='approximate';


%cst
cst_p.sec_in_day_cst=86400;
cst_p.c_cst=299792458.0; 
cst_p.flat_coeff_cst=0.003352810664747;
cst_p.semi_major_axis_cst = 6378137.0;
%read netcdf file
data_GPP = readL1B_S6_ISD(filename_L1B_GPP,cnf_p,cst_p);
N_surfaces=length(data_GPP.GEO.TAI.total);
%read mat file
data_MAT = readL1B_S6_ISD(filename_L1B_MAT,cnf_p,cst_p);



if TRP
    idx_int=[28,122,216,309];
else
    idx_int=1:N_surfaces;
end

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

%% ------------ GEO: geographical information -----------------------------

%time_tai_20_ku
figure('Name','time_tai_20_ku');
subplot(2,1,1)
plot(1:N_surfaces,data_MAT.GEO.TAI.total(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.TAI.total,'-r'); 
legend('MAT','GPP')
title('time_tai_20_ku','Interpreter','none');
ylabel('[s]'); xlabel('Surface'); 
subplot(2,1,2)
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.TAI.total-data_MAT.GEO.TAI.total(idx_int)),...
       1:N_surfaces,abs(data_GPP.GEO.TAI.total-data_MAT.GEO.TAI.total(idx_int))./(data_MAT.GEO.TAI.total(idx_int))*100);
title('Error time_tai_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),'[s]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'time_tai_20_ku_comparison',file_ext)); 


%latitude_20_ku
figure('Name','latitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.LAT(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.LAT,'-r'); 
legend('MAT','GPP')
title('latitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.LAT-data_MAT.GEO.LAT(idx_int)),...
       1:N_surfaces,abs((data_GPP.GEO.LAT-data_MAT.GEO.LAT(idx_int))./(data_MAT.GEO.LAT(idx_int)))*100);
title('Error latitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'latitude_20_ku_comparison',file_ext)); 


%longitude_20_ku
figure('Name','longitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.LON(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.LON,'-r'); 
legend('MAT','GPP')
title('longitude_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.LON-data_MAT.GEO.LON(idx_int)),...
       1:N_surfaces,abs((data_GPP.GEO.LON-data_MAT.GEO.LON(idx_int))./(data_MAT.GEO.LON(idx_int)))*100);
title('Error longitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'longitude_20_ku_comparison',file_ext)); 


%com_altitude_rate_20_ku
figure('Name','com_altitude_rate_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.H_rate(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.H_rate,'-r'); 
legend('MAT','GPP')
title('com_altitude_rate_20_ku','Interpreter','none');
ylabel('[m/s]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.H_rate-data_MAT.GEO.H_rate(idx_int)),...
       1:N_surfaces,abs((data_GPP.GEO.H_rate-data_MAT.GEO.H_rate(idx_int))./(data_MAT.GEO.H_rate(idx_int)))*100);
title('Error com_altitude_rate_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [m/s]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'com_altitude_rate_20_ku_comparison',file_ext)); 


%norm com_velocity_vector_20_ku
figure('Name','norm com_velocity_vector_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.V(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.V,'-r'); 
legend('MAT','GPP')
title('com_velocity_vector_20_ku','Interpreter','none');
ylabel('[m/s]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.V-data_MAT.GEO.V(idx_int)),...
       1:N_surfaces,abs((data_GPP.GEO.V-data_MAT.GEO.V(idx_int))./(data_MAT.GEO.V(idx_int)))*100);
title('Error norm com_velocity_vector_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [m/s]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'com_velocity_vector_20_ku_comparison',file_ext)); 


%com_altitude_20_ku
figure('Name','com_altitude_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.H(idx_int),'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.H,'-r'); 
legend('MAT','GPP')
title('com_altitude_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.H-data_MAT.GEO.H(idx_int)),...
       1:N_surfaces,abs((data_GPP.GEO.H-data_MAT.GEO.H(idx_int))./(data_MAT.GEO.H(idx_int)))*100);
title('Error com_altitude_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'com_altitude_20_ku_comparison',file_ext)); 


%off_nadir_pitch_angle_pf_20_ku
figure('Name','off_nadir_pitch_angle_pf_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.pitch(idx_int).*180/pi,'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.pitch.*180/pi,'-r'); 
legend('MAT','GPP')
title('off_nadir_pitch_angle_pf_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.pitch-data_MAT.GEO.pitch(idx_int))*180/pi,...
       1:N_surfaces,abs((data_GPP.GEO.pitch-data_MAT.GEO.pitch(idx_int))./(data_MAT.GEO.pitch(idx_int)))*100);
title('Error off_nadir_pitch_angle_pf_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'off_nadir_pitch_angle_pf_20_ku_comparison',file_ext)); 


%off_nadir_roll_angle_pf_20_ku
figure('Name','off_nadir_roll_angle_pf_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.roll(idx_int).*180/pi,'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.roll.*180/pi,'-r'); 
legend('MAT','GPP')
title('off_nadir_roll_angle_pf_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.roll-data_MAT.GEO.roll(idx_int))*180/pi,...
       1:N_surfaces,abs((data_GPP.GEO.roll-data_MAT.GEO.roll(idx_int))./(data_MAT.GEO.roll(idx_int)))*100);
title('Error off_nadir_roll_angle_pf_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'off_nadir_roll_angle_pf_20_ku_comparison',file_ext)); 


%off_nadir_yaw_angle_pf_20_ku
figure('Name','off_nadir_yaw_angle_pf_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.GEO.yaw(idx_int).*180/pi,'-b'); 
hold on; plot(1:N_surfaces,data_GPP.GEO.yaw.*180/pi,'-r'); 
legend('MAT','GPP')
title('off_nadir_yaw_angle_pf_20_ku','Interpreter','none');
ylabel('[deg]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.GEO.yaw-data_MAT.GEO.yaw(idx_int))*180/pi,...
       1:N_surfaces,abs((data_GPP.GEO.yaw-data_MAT.GEO.yaw(idx_int))./(data_MAT.GEO.yaw(idx_int)))*100);
title('Error off_nadir_yaw_angle_pf_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [deg]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'off_nadir_yaw_angle_pf_20_ku_comparison',file_ext)); 


%close all;

%% ------------ MEA: measurements -----------------------------------------
%tracker_range_calibrated_20_ku
figure('Name','tracker_range_calibrated_20_ku');
subplot(2,1,1);
plot(1:N_surfaces,data_MAT.MEA.win_delay(idx_int).*c_cst/2,'-b'); 
hold on; plot(1:N_surfaces,data_GPP.MEA.win_delay.*c_cst/2,'-r'); 
legend('MAT','GPP')
title('tracker_range_calibrated_20_ku','Interpreter','none');
ylabel('[m]'); xlabel('Surface'); 
subplot(2,1,2);
[hAX]=plotyy(1:N_surfaces,(data_GPP.MEA.win_delay-data_MAT.MEA.win_delay(idx_int)).*c_cst/2,...
       1:N_surfaces,((data_GPP.MEA.win_delay-data_MAT.MEA.win_delay(idx_int))./(data_MAT.MEA.win_delay(idx_int)))*100);
title('Error tracker_range_calibrated_20_ku','Interpreter','none');
xlabel('Surface'); 
ylabel(hAX(1),' [m]');
ylabel(hAX(2),'Relative [%]');
print(print_file,res_fig,strcat(resultPath,'tracker_range_calibrated_20_ku_comparison',file_ext)); 


%power waveforms hr_power_waveform_20_ku
max_image=10*log10(max([max(max(data_MAT.HRM.power_wav(:,idx_int))),...
                        max(max(data_GPP.HRM.power_wav))]));
min_image=max_image-40.0;                    
figure('Name','Power Waveforms hr_power_waveform_20_ku'); 
subplot(1,2,1)
imagesc(10*log10(data_MAT.HRM.power_wav(:,idx_int))); 
colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
caxis([min_image,max_image]);
title('MAT','Interpreter','none')
xlabel('Surfaces'); ylabel('Samples');
subplot(1,2,2)
imagesc(10*log10(data_GPP.HRM.power_wav)); 
colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
caxis([min_image,max_image]);
title('GPP','Interpreter','none')
xlabel('Surfaces'); ylabel('Samples');
[axT,hT]=suplabel('Power waveforms comparison hr_power_waveform_20_ku','t');
print(print_file,res_fig,strcat(resultPath,'Power_waveforms_comparison',file_ext)); 


%individual comparative wvfms analysis
if TRP
    for i_trp=1:length(idx_int)
        figure('Name',strcat('Wvfm TRP #',num2str(i_trp)));
        plot(10*log10(data_MAT.HRM.power_wav(:,idx_int(i_trp))),'-b');
        hold on;
        plot(10*log10(data_GPP.HRM.power_wav(:,i_trp)),'-r');
        ylim([min_image,max_image]);
        title(strcat('Wvfm TRP #',num2str(i_trp)));
        ylabel('[dBW]'); xlabel('Range bin');
        legend('MAT','GPP');
    end
end

figure('Name','Error on waveforms hr_power_waveform_20_ku'); 
subplot(1,2,1)
imagesc(abs(data_GPP.HRM.power_wav-data_MAT.HRM.power_wav(:,idx_int))); 
colormap(colormap_set); c=colorbar; ylabel(c,'[Watts]');
title('Absolute','Interpreter','none')
xlabel('Surfaces'); ylabel('Samples');
subplot(1,2,2)
imagesc(abs((data_GPP.HRM.power_wav-data_MAT.HRM.power_wav(:,idx_int))./(data_MAT.HRM.power_wav(:,idx_int))).*100); 
colormap(colormap_set); c=colorbar; ylabel(c,'[%]');
caxis([0,100]);
title('Relative','Interpreter','none')
xlabel('Surfaces'); ylabel('Samples');
[axT,hT]=suplabel('Error on power waveforms hr_power_waveform_20_ku','t');
print(print_file,res_fig,strcat(resultPath,'Error_Power_waveforms_comparison',file_ext)); 


%power waveforms stack_mask_start_stop_20_ku
idx_rows=any(data_GPP.HRM.Doppler_mask-1,2).';
figure('Name','Error on stack_mask_start_stop_20_ku'); 
subplot(1,2,1)
imagesc(abs(data_GPP.HRM.Doppler_mask(idx_rows,:)-data_MAT.HRM.Doppler_mask(:,idx_int))); 
colormap(colormap_set); c=colorbar; ylabel(c,'[samples]');
caxis([0,2]);
title('Absolute','Interpreter','none')
ylabel('Beams'); xlabel('Surfaces');

subplot(1,2,2)
imagesc(abs((data_GPP.HRM.Doppler_mask(idx_rows,:)-data_MAT.HRM.Doppler_mask(:,idx_int))./(data_MAT.HRM.Doppler_mask(:,idx_int))).*100); 
colormap(colormap_set); c=colorbar; ylabel(c,'[%]');
caxis([0,100]);
title('Relative','Interpreter','none')
ylabel('Beams'); xlabel('Surfaces');
[axT,hT]=suplabel('Error on power waveforms stack_mask_start_stop_20_ku','t');
print(print_file,res_fig,strcat(resultPath,'stack_mask_start_stop_20_ku_comparison',file_ext)); 

if strcmp(cnf_p.mode,'RAW') || strcmp(cnf_p.mode,'RMC')
    %look_angle_start_20_ku
    figure('Name','look_angle_start_20_ku');
    subplot(2,1,1);
    plot(1:N_surfaces,data_MAT.HRM.look_ang_start_surf(idx_int).*180/pi,'-b');
    hold on; plot(1:N_surfaces,data_GPP.HRM.look_ang_start_surf.*180/pi,'-r');
    legend('MAT','GPP')
    title('look_angle_start_20_ku','Interpreter','none');
    ylabel('[deg]'); xlabel('Surface');
    subplot(2,1,2);
    [hAX]=plotyy(1:N_surfaces,(data_GPP.HRM.look_ang_start_surf-data_MAT.HRM.look_ang_start_surf(idx_int)).*180/pi,...
        1:N_surfaces,abs((data_GPP.HRM.look_ang_start_surf-data_MAT.HRM.look_ang_start_surf(idx_int))./(data_MAT.HRM.look_ang_start_surf(idx_int)))*100);
    title('Error look_angle_start_20_ku','Interpreter','none');
    xlabel('Surface');
    ylabel(hAX(1),' [deg]');
    ylabel(hAX(2),'Relative [%]');
    print(print_file,res_fig,strcat(resultPath,'look_angle_start_20_ku_comparison',file_ext));
    
    
    %look_angle_stop_20_ku
    figure('Name','look_angle_stop_20_ku');
    subplot(2,1,1);
    plot(1:N_surfaces,data_MAT.HRM.look_ang_stop_surf(idx_int).*180/pi,'-b');
    hold on; plot(1:N_surfaces,data_GPP.HRM.look_ang_stop_surf.*180/pi,'-r');
    legend('MAT','GPP')
    title('look_angle_stop_20_ku','Interpreter','none');
    ylabel('[deg]'); xlabel('Surface');
    subplot(2,1,2);
    [hAX]=plotyy(1:N_surfaces,(data_GPP.HRM.look_ang_stop_surf-data_MAT.HRM.look_ang_stop_surf(idx_int)).*180/pi,...
        1:N_surfaces,abs((data_GPP.HRM.look_ang_stop_surf-data_MAT.HRM.look_ang_stop_surf(idx_int))./(data_MAT.HRM.look_ang_stop_surf(idx_int)))*100);
    title('Error look_angle_stop_20_ku','Interpreter','none');
    xlabel('Surface');
    ylabel(hAX(1),' [deg]');
    ylabel(hAX(2),'Relative [%]');
    print(print_file,res_fig,strcat(resultPath,'look_angle_stop_20_ku_comparison',file_ext));
    
end


% %% -------------- Comparison RAW and RMC data -----------------------------
% filename_L1B_GPP = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\l1b_hr_product_lat00_CFI.nc';
% filename_L1B_MAT = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\l1b_hr_product_lat00.mat';
% 
% filename_L1B_RMC_GPP = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\l1b_hr_rmc_product_lat00.nc';
% filename_L1B_RMC_MAT = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\l1b_hr_rmc_product_lat0.mat';
% 
% path_results = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\comparison_results\';
% mkdir(path_results);
% 
% global c_cst sec_in_day_cst
% c_cst=299792458;
% sec_in_day_cst=86400;
% 
% %configuration
% cnf_p.retracker_name='ANALYTICAL';
% cnf_p.looks_index_method='Look_angle';
% cnf_p.look_ang_method='approximate';
% 
% %--- raw data ---------- read netcdf file
% data_GPP = readL1B_S6_ISD (filename_L1B_GPP,cnf_p);
% %read mat file
% data_MAT = readL1B_S6_ISD (filename_L1B_MAT,cnf_p);
% 
% 
% % -- rmc data ---------- read netcdf file
% data_RMC_GPP = readL1B_S6_ISD (filename_L1B_RMC_GPP,cnf_p);
% %read mat file
% data_RMC_MAT = readL1B_S6_ISD (filename_L1B_RMC_MAT,cnf_p);
% 
% %data_RMC_MAT.HRM.power_wav = data_RMC_MAT.HRM.power_wav.*(10^(124.39/10)/(10^(88.27/10)));
% 
% set_default_plot;
% set(0,'defaultFigureVisible','on');
% 
% 
% figure_format='jpg';
% res_fig='-r300';
% 
% switch lower(figure_format)
%     case 'pdf'
%         file_ext='.pdf';
%         print_file='-dpdf';
%     case 'eps'
%         file_ext='.eps';
%         print_file='-depsc';
%     case 'png'
%         file_ext='.png';
%         print_file='-dpng';
%     case 'jpg'
%         file_ext='.jpg';
%         print_file='-djpeg';        
% end
% 
% 
% 
% %POWER WAVEFORM COMPARISON
% max_image=10*log10(max([max(data_MAT.HRM.power_wav(:)),...
%                    max(data_GPP.HRM.power_wav(:)),...
%                    max(data_RMC_MAT.HRM.power_wav(:)),...
%                    max(data_RMC_GPP.HRM.power_wav(:))]));
% min_image=max_image-30.0;
% figure('Name','Power Waveforms'); 
% subplot(2,2,1)
% imagesc(10*log10(data_MAT.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('RAW-MAT','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,2)
% imagesc(10*log10(data_GPP.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('RAW-GPP','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,3)
% imagesc(10*log10(data_RMC_MAT.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('RMC-MAT','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,4)
% imagesc(10*log10(data_RMC_GPP.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('RMC-GPP','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% [axT,hT]=suplabel('Power waveforms comparison','t');
% 
% set(0,'defaultFigureVisible','off');
% 
% %waveform a waveform comparison
% num_pools=4;
% parpool(num_pools);
% 
% parfor i_surf=1:data_GPP.N_records
%     
%     f1=figure;
%     plot(10*log10(data_MAT.HRM.power_wav(:,i_surf)),'-b');
%     hold on;
%     plot(10*log10(data_GPP.HRM.power_wav(:,i_surf)),'-r');
%     plot(10*log10(data_RMC_MAT.HRM.power_wav(:,i_surf)),'-g');
%     plot(10*log10(data_RMC_GPP.HRM.power_wav(:,i_surf)),'-m');
%     ylim([min_image,max_image]);
%     legend('RAW-MAT','RAW-GPP','RMC-MAT','RMC-GPP');
%     xlabel('Samples'); ylabel('[dBW]');
%     title(strcat('L1B wvfm comparison surf #',num2str(i_surf),...
%         '(LAT=',num2str(data_GPP.GEO.LAT(i_surf)),'[deg.])'))
%     print(print_file,res_fig,strcat(path_results,'stack_comparison_surf_',...
%         num2str(i_surf,'%04.0f'),file_ext));
%     close(f1)
% end
% 
% poolobj = gcp('nocreate');
% delete(poolobj);


% %% -------------- Comparison RAW zeros/no zeros ML data -----------------------------
% filename_L1B_GPP = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\RAW\Z_ML\l1b_hr_product_lat00_CFI.nc';
% filename_L1B_MAT = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\RAW\Z_ML\l1b_hr_product_lat00.mat';
% 
% filename_L1B_NZ_GPP = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\RAW\NZ_ML\l1b_hr_raw_product_lat00_zeros_ml.nc';
% filename_L1B_NZ_MAT = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\RAW\NZ_ML\l1b_hr_product_lat0_NZ_ML.mat';
% 
% path_results = 'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\inputs\L1B_GPP\lat0\RAW\comparison_results\';
% mkdir(path_results);
% 
% global c_cst sec_in_day_cst
% c_cst=299792458;
% sec_in_day_cst=86400;
% 
% %configuration
% cnf_p.retracker_name='ANALYTICAL';
% cnf_p.looks_index_method='Look_angle';
% cnf_p.look_ang_method='approximate';
% cnf_p.mode = 'RAW';
% 
% %--- raw data ---------- read netcdf file
% data_GPP = readL1B_S6_ISD (filename_L1B_GPP,cnf_p);
% %read mat file
% data_MAT = readL1B_S6_ISD (filename_L1B_MAT,cnf_p);
% 
% 
% % -- rmc data ---------- read netcdf file
% data_NZ_GPP = readL1B_S6_ISD (filename_L1B_NZ_GPP,cnf_p);
% %read mat file
% data_NZ_MAT = readL1B_S6_ISD (filename_L1B_NZ_MAT,cnf_p);
% 
% %data_RMC_MAT.HRM.power_wav = data_RMC_MAT.HRM.power_wav.*(10^(124.39/10)/(10^(88.27/10)));
% 
% set_default_plot;
% set(0,'defaultFigureVisible','on');
% 
% 
% figure_format='jpg';
% res_fig='-r300';
% 
% switch lower(figure_format)
%     case 'pdf'
%         file_ext='.pdf';
%         print_file='-dpdf';
%     case 'eps'
%         file_ext='.eps';
%         print_file='-depsc';
%     case 'png'
%         file_ext='.png';
%         print_file='-dpng';
%     case 'jpg'
%         file_ext='.jpg';
%         print_file='-djpeg';        
% end
% 
% 
% 
% %POWER WAVEFORM COMPARISON
% max_image=10*log10(max([max(data_MAT.HRM.power_wav(:)),...
%                    max(data_GPP.HRM.power_wav(:)),...
%                    max(data_NZ_MAT.HRM.power_wav(:)),...
%                    max(data_NZ_GPP.HRM.power_wav(:))]));
% min_image=max_image-30.0;
% figure('Name','Power Waveforms'); 
% subplot(2,2,1)
% imagesc(10*log10(data_MAT.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('MAT (zeros ML)','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,2)
% imagesc(10*log10(data_GPP.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('GPP (zeros ML)','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,3)
% imagesc(10*log10(data_NZ_MAT.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('MAT (no-zeros ML)','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% subplot(2,2,4)
% imagesc(10*log10(data_NZ_GPP.HRM.power_wav)); 
% colormap(colormap_set); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% title('GPP (no-zeros ML)','Interpreter','none')
% xlabel('Surfaces'); ylabel('Samples');
% [axT,hT]=suplabel('Power waveforms comparison','t');
% 
% set(0,'defaultFigureVisible','off');
% 
% %waveform a waveform comparison
% num_pools=2;
% parpool(num_pools);
% 
% parfor i_surf=1:data_GPP.N_records
%     
%     f1=figure;
%     plot(10*log10(data_MAT.HRM.power_wav(:,i_surf)),'-b');
%     hold on;
%     plot(10*log10(data_GPP.HRM.power_wav(:,i_surf)),'-r');
%     plot(10*log10(data_NZ_MAT.HRM.power_wav(:,i_surf)),'-g');
%     plot(10*log10(data_NZ_GPP.HRM.power_wav(:,i_surf)),'-m');
%     ylim([min_image,max_image]);
%     legend('MAT (zeros ML)','GPP (zeros ML)','MAT (no-zeros ML)','GPP (no-zeros ML)');
%     xlabel('Samples'); ylabel('[dBW]');
%     title(strcat('L1B wvfm comparison surf #',num2str(i_surf),...
%         '(LAT=',num2str(data_GPP.GEO.LAT(i_surf)),'[deg.])'))
%     print(print_file,res_fig,strcat(path_results,'wvfm_comparison_cpp_mat_NZ_ML_',...
%         num2str(i_surf,'%04.0f'),file_ext));
%     close(f1)
% end
% 
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% 
% %rename the waveform fitted plots
% input_path='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\L2\lat40\GPP\RMC\plots\fitted_waveforms\';
% file_ext='.png';
% token_break_string='wvfm_';
% files=dir([input_path '*' file_ext]);
% 
% for i_file=1:1%length(files)
%     filename_string = strrep(char(files(i_file).name),file_ext,'');
%     number = strsplit(filename_string,token_break_string);
%     filename_string_nonumber = char(number(1)); 
%     number = num2str(str2num(char(number(end))),'%04.0f');
%     movefile([input_path char(files(i_file).name)],...
%         [input_path strcat(filename_string_nonumber,'wvfm_',number,file_ext)]);
% end
