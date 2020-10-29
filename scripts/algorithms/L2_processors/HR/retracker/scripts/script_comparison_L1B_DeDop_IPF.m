%% -------------- Files defintion for comparison --------------------------
file_DeDop ='C:\Users\eduard.makhoul\isardSAT\projects\DeDop\data\L1B_issues_data_DeDop\data\s3_ocean\L1B_SR_1_measurement_l1a_default.nc';

file_IPF ='C:\Users\eduard.makhoul\isardSAT\projects\DeDop\data\L1B_issues_data_DeDop\data\s3_ocean\S3A_SR_1_SRA____20171105T054940_20171105T064009_20171130T221240_3029_024_119______LN3_O_NT_002.SEN3\measurement.nc';

resultPath ='C:\Users\eduard.makhoul\isardSAT\projects\DeDop\data\L1B_issues_data_DeDop\data\s3_ocean\plots_comparison_DeDop_IPF\'
mkdir(resultPath);

file_CST='C:\Users\eduard.makhoul\isardSAT\projects\DeDop\data\L1B_issues_data_DeDop\configuration_L2\cst_file.json';
file_CHD='C:\Users\eduard.makhoul\isardSAT\projects\DeDop\data\L1B_issues_data_DeDop\configuration_L2\chd_file_S3.json'

%configuration
cnf_p.retracker_name='ANALYTICAL';
cnf_p.looks_index_method='Look_angle';
cnf_p.look_ang_method='approximate';
cnf_p.mode = 'SAR';
cnf_p.input_type =1;
cnf_p.window_type_r ='Boxcar';
cnf_p.window_type_a ='Boxcar';
cnf_p.Doppler_mask_cons_option ='Internal';

%cst: constants
cst_p=read_CST_json(file_CST);
%characterization: characterization
chd_p=read_CHD_json(file_CHD,cnf_p,cst_p);

figure_format='jpg';
res_fig='-r150';



%% -------------- Reading variables of interest ---------------------------


cnf_p.L1proc = 'DeDop';
data_DeDop = readL1B_S3_DeDop (file_DeDop,cnf_p,cst_p,chd_p);
cnf_p.L1proc = 'ESA';
data_IPF   = readL1B_S3_ESA (file_IPF,cnf_p,cst_p,chd_p);


%% ------------------- Plots ----------------------------------------------
set_default_plot;
set(0,'defaultFigureVisible','on');

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
%% ------------------- Geolocation of surfaces ----------------------------
%Compute the positions of the surfaces in cartesian coordinates
%x,y,z
%assuming surface height is zero
IPF_xyz_surf = lla2ecef([data_IPF.GEO.LAT.',data_IPF.GEO.LON.',zeros(data_IPF.N_records,1)],...
                        cst_p.flat_coeff_cst,cst_p.semi_major_axis_cst);
DeDop_xyz_surf = lla2ecef([data_DeDop.GEO.LAT.',data_DeDop.GEO.LON.',zeros(data_DeDop.N_records,1)],...
                        cst_p.flat_coeff_cst,cst_p.semi_major_axis_cst);                    

%consider the N first DeDop surfaces for IPF
idx_DeDop = (zeros(1,data_DeDop.N_records));                    
idx_IPF = (zeros(1,data_IPF.N_records));                    

if data_DeDop.N_records<=data_IPF.N_records
    idx_IPF(1:data_DeDop.N_records-1)=1;
    idx_DeDop(2:data_DeDop.N_records)=1;
else
    idx_IPF(1:data_IPF.N_records)=1;
    idx_DeDop(2:data_IPF.N_records+1)=1;
end
idx_DeDop =logical(idx_DeDop);
idx_IPF =logical(idx_IPF);
N_surfaces_common=length(find(idx_DeDop));
    
%***************** differences geolocation in meters **********************
GEO_loc_errors=sqrt(sum((IPF_xyz_surf(idx_IPF,:)-DeDop_xyz_surf(idx_DeDop,:)).^2,2));

f1=figure; plot(data_DeDop.GEO.LAT(idx_DeDop),GEO_loc_errors,'*b');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Geolocation errors: IPF & DeDop');
xlabel('Latitude DeDop [deg.]'); ylabel('Distance [m]');
print(print_file,res_fig,strcat(resultPath,'Geolocation_error_vs_LAT',file_ext)); 
close(f1);

f1=figure; scatter3(data_DeDop.GEO.LON(idx_DeDop),data_DeDop.GEO.LAT(idx_DeDop),...
    GEO_loc_errors/1e3,5,GEO_loc_errors/1e3); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[Km]'); caxis([0,5])
title(strcat('Geolocation errors: IPF & DeDop [Km]')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
print(print_file,res_fig,strcat(resultPath,'Geolocation_error_vs_LAT_LON',file_ext)); 
close(f1)

%********************  Track over MAP *************************************
f1=figure; geoshow('landareas.shp','FaceColor', rgb('green'));
geoshow(data_DeDop.GEO.LAT,data_DeDop.GEO.LON,'DisplayType','Point','Marker','.','MarkerEdgeColor','blue');
geoshow(data_IPF.GEO.LAT,data_IPF.GEO.LON,'DisplayType','Point','Marker','o','MarkerEdgeColor','red');
title('Surface track')
xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
legend('DeDop','IPF');
print(print_file,res_fig,strcat(resultPath,'Geolocation_map',file_ext)); 
close(f1)

% *********************** separation surfaces *****************************
DeDop_surf_dis=sqrt(sum((diff(DeDop_xyz_surf(idx_DeDop,:),1,1)).^2,2)); % separation between consecutive surfaces
IPF_surf_dis=sqrt(sum((diff(IPF_xyz_surf(idx_IPF,:),1,1)).^2,2));

f1=figure; plot(data_DeDop.GEO.LAT(idx_DeDop(1:end-1)),DeDop_surf_dis,'*r');
hold on; plot(data_DeDop.GEO.LAT(idx_DeDop(1:end-1)),IPF_surf_dis,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Surface separation: IPF & DeDop [m]');
xlabel('Latitude DeDop [deg.]'); ylabel('Distance [m]');
ylim([320,350]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Separation_surfaces_vs_LAT',file_ext)); 
close(f1);

f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON(idx_DeDop(1:end-1)),data_DeDop.GEO.LAT(idx_DeDop(1:end-1)),...
    DeDop_surf_dis,5,DeDop_surf_dis); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m]'); caxis([320,350])
title(strcat('Surface separation: DeDop [m]')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_DeDop.GEO.LON(idx_DeDop(1:end-1)),data_DeDop.GEO.LAT(idx_DeDop(1:end-1)),...
    IPF_surf_dis,5,IPF_surf_dis); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m]'); caxis([320,350])
title(strcat('Surface separation: IPF [m]')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
print(print_file,res_fig,strcat(resultPath,'Separation_surfaces_vs_LAT_LON',file_ext)); 
close(f1)

%% ----------------------------- Orbital Height ---------------------------
max_image =max([max(data_DeDop.GEO.H(:)),max(data_IPF.GEO.H(:))])/1e3;
min_image =min([min(data_DeDop.GEO.H(:)),min(data_IPF.GEO.H(:))])/1e3;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.GEO.H/1e3,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.GEO.H/1e3,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Satellite altitude: IPF & DeDop [Km]');
xlabel('Latitude [deg.]'); ylabel('[Km]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Satellite_altitude_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.GEO.H/1e3,5,data_DeDop.GEO.H/1e3); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[Km]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.GEO.H/1e3,5,data_IPF.GEO.H/1e3); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[Km]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Satellite altitude','t');
print(print_file,res_fig,strcat(resultPath,'Satellite_altitude_vs_LAT_LON',file_ext)); 
close(f1)

%% ----------------------------- Tracker range ---------------------------
max_image =max([max(data_DeDop.MEA.win_delay(:)),max(data_IPF.MEA.win_delay(:))])*cst_p.c_cst/2/1e3;
min_image =min([min(data_DeDop.MEA.win_delay(:)),min(data_IPF.MEA.win_delay(:))])*cst_p.c_cst/2/1e3;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.MEA.win_delay*cst_p.c_cst/2/1e3,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.MEA.win_delay*cst_p.c_cst/2/1e3,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Tracker range: IPF & DeDop [Km]');
xlabel('Latitude [deg.]'); ylabel('[Km]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Tracker_range_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.MEA.win_delay*cst_p.c_cst/2/1e3,5,data_DeDop.MEA.win_delay*cst_p.c_cst/2/1e3); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[Km]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.MEA.win_delay*cst_p.c_cst/2/1e3,5,data_IPF.MEA.win_delay*cst_p.c_cst/2/1e3); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[Km]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Tracker range','t');
print(print_file,res_fig,strcat(resultPath,'Tracker_range_vs_LAT_LON',file_ext)); 
close(f1)

%% ----------------------------- Surface Height ---------------------------
max_image =max([max(data_DeDop.GEO.H(:)-data_DeDop.MEA.win_delay(:)*cst_p.c_cst/2),max(data_IPF.GEO.H(:)-data_IPF.MEA.win_delay(:)*cst_p.c_cst/2)]);
min_image =min([min(data_DeDop.GEO.H(:)-data_DeDop.MEA.win_delay(:)*cst_p.c_cst/2),min(data_IPF.GEO.H(:)-data_IPF.MEA.win_delay(:)*cst_p.c_cst/2)]);

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.GEO.H(:)-data_DeDop.MEA.win_delay(:)*cst_p.c_cst/2,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.GEO.H(:)-data_IPF.MEA.win_delay(:)*cst_p.c_cst/2,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Orbital height - Tracker range: IPF & DeDop [m]');
xlabel('Latitude [deg.]'); ylabel('[m]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Surface_height_vs_LAT',file_ext)); 
close(f1);

max_image=50.0;
min_image=-50.0;
f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.GEO.H(:)-data_DeDop.MEA.win_delay(:)*cst_p.c_cst/2,5,...
    data_DeDop.GEO.H(:)-data_DeDop.MEA.win_delay(:)*cst_p.c_cst/2); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.GEO.H(:)-data_IPF.MEA.win_delay(:)*cst_p.c_cst/2,5,...
    data_IPF.GEO.H(:)-data_IPF.MEA.win_delay(:)*cst_p.c_cst/2); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Orbital height - Tracker range','t');
print(print_file,res_fig,strcat(resultPath,'Surface_height_vs_LAT_LON',file_ext)); 
close(f1)

%% ----------------------------- Height rate ---------------------------
max_image =max([max(data_DeDop.GEO.H_rate(:)),max(data_IPF.GEO.H_rate(:))]);
min_image =min([min(data_DeDop.GEO.H_rate(:)),min(data_IPF.GEO.H_rate(:))]);

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.GEO.H_rate,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.GEO.H_rate,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Height ratge: IPF & DeDop [m/s]');
xlabel('Latitude [deg.]'); ylabel('[m/s]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Height_rate_vs_LAT',file_ext)); 
close(f1);

f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.GEO.H_rate,5,...
    data_DeDop.GEO.H_rate); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m/s]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.GEO.H_rate,5,...
    data_IPF.GEO.H_rate); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[m/s]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Height rate','t');
print(print_file,res_fig,strcat(resultPath,'Height_rate_vs_LAT_LON',file_ext)); 
close(f1)

%% ------------------------- Number of contributing beams -----------------
max_image =max([max(data_DeDop.HRM.Neff),max(data_IPF.HRM.Neff)]);
min_image =min([min(data_DeDop.HRM.Neff),min(data_IPF.HRM.Neff)]);

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.Neff,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.Neff,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Number of contributing beams: IPF & DeDop [m]');
xlabel('Latitude [deg.]'); ylabel('# Beams');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Contributing_beams_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.Neff,5,...
    data_DeDop.HRM.Neff); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'# Beams'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.Neff,5,...
    data_IPF.HRM.Neff); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'# Beams'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Number of contributing beams','t');
print(print_file,res_fig,strcat(resultPath,'Contributing_beams_vs_LAT_LON',file_ext)); 
close(f1)

%% ------------------------- Look angles -----------------
%Start
max_image =max([max(data_DeDop.HRM.look_ang_start_surf),max(data_IPF.HRM.look_ang_start_surf)])*180/cst_p.pi_cst;
min_image =min([min(data_DeDop.HRM.look_ang_start_surf),min(data_IPF.HRM.look_ang_start_surf)])*180/cst_p.pi_cst;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.look_ang_start_surf*180/cst_p.pi_cst,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Look angle start (beam_ang_l1b_echo_sar_ku): IPF & DeDop [deg.]', 'Interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[deg.]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Look_angle_start_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.look_ang_start_surf*180/cst_p.pi_cst,5,...
    data_DeDop.HRM.look_ang_start_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst,5,...
    data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Look angle start (beam_ang_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Look_angle_start_vs_LAT_LON',file_ext)); 
close(f1)

%Stop
max_image =max([max(data_DeDop.HRM.look_ang_stop_surf),max(data_IPF.HRM.look_ang_stop_surf)])*180/cst_p.pi_cst;
min_image =min([min(data_DeDop.HRM.look_ang_stop_surf),min(data_IPF.HRM.look_ang_stop_surf)])*180/cst_p.pi_cst;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.look_ang_stop_surf*180/cst_p.pi_cst,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Look angle stop (beam_ang_l1b_echo_sar_ku): IPF & DeDop [deg.]', 'Interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[deg.]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Look_angle_stop_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.look_ang_stop_surf*180/cst_p.pi_cst,5,...
    data_DeDop.HRM.look_ang_stop_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst,5,...
    data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Look angle start (beam_ang_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Look_angle_stop_vs_LAT_LON',file_ext)); 
close(f1)

% ********** computing the look angles for DeDop from height rate *********
Dopp_angle = asin(data_DeDop.GEO.H_rate./data_DeDop.GEO.V);
data_DeDop.HRM.look_ang_start_surf_recomp=(cst_p.pi_cst/2-(data_DeDop.HRM.look_ang_start_surf-Dopp_angle));
data_DeDop.HRM.look_ang_stop_surf_recomp=(cst_p.pi_cst/2-(data_DeDop.HRM.look_ang_stop_surf-Dopp_angle));

%Start
max_image =max([max(data_DeDop.HRM.look_ang_start_surf_recomp),max(data_IPF.HRM.look_ang_start_surf)])*180/cst_p.pi_cst;
min_image =min([min(data_DeDop.HRM.look_ang_start_surf_recomp),min(data_IPF.HRM.look_ang_start_surf)])*180/cst_p.pi_cst;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.look_ang_start_surf_recomp*180/cst_p.pi_cst,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Look angle start (beam_ang_l1b_echo_sar_ku): IPF & DeDop [deg.]', 'Interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[deg.]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Look_angle_start_corr_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.look_ang_start_surf_recomp*180/cst_p.pi_cst,5,...
    data_DeDop.HRM.look_ang_start_surf_recomp*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst,5,...
    data_IPF.HRM.look_ang_start_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Look angle start (beam_ang_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Look_angle_start_corr_vs_LAT_LON',file_ext)); 
close(f1)

%Stop
max_image =max([max(data_DeDop.HRM.look_ang_stop_surf_recomp),max(data_IPF.HRM.look_ang_stop_surf)])*180/cst_p.pi_cst;
min_image =min([min(data_DeDop.HRM.look_ang_stop_surf_recomp),min(data_IPF.HRM.look_ang_stop_surf)])*180/cst_p.pi_cst;

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.look_ang_stop_surf_recomp*180/cst_p.pi_cst,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Look angle start (beam_ang_l1b_echo_sar_ku): IPF & DeDop [deg.]', 'Interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[deg.]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Look_angle_stop_corr_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.look_ang_stop_surf_recomp*180/cst_p.pi_cst,5,...
    data_DeDop.HRM.look_ang_stop_surf_recomp*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst,5,...
    data_IPF.HRM.look_ang_stop_surf*180/cst_p.pi_cst); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[deg.]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Look angle start (beam_ang_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Look_angle_stop_corr_vs_LAT_LON',file_ext)); 
close(f1)


%% ---------------------- Waveforms ---------------------------------------

% ************** Peak of the waveforms ************************************
max_image=10*log10(max([max(data_DeDop.HRM.power_wav(:)),max(data_IPF.HRM.power_wav(:))]));
min_image=max_image-80.0;

f1=figure; 
plot(data_DeDop.GEO.LAT,10*log10(max(data_DeDop.HRM.power_wav,[],1)),'*r');
hold on; plot(data_IPF.GEO.LAT,10*log10(max(data_IPF.HRM.power_wav,[],1)),'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('Peak waveforms across-track (i2q2_meas_ku_l1b_echo_sar_ku): IPF & DeDop','interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[dBW]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Peak_waveform_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    10*log10(max(data_DeDop.HRM.power_wav,[],1)),5,...
    10*log10(max(data_DeDop.HRM.power_wav,[],1))); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    10*log10(max(data_IPF.HRM.power_wav,[],1)),5,...
    10*log10(max(data_IPF.HRM.power_wav,[],1))); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('Peak waveforms across-track (i2q2_meas_ku_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Peak_waveform_vs_LAT_LON',file_ext)); 
close(f1)

% ---------------------- AGC  -----------------------------------
max_image=max([max(data_DeDop.HRM.AGC(:)),max(data_IPF.HRM.AGC(:))]);
min_image=min([min(data_DeDop.HRM.AGC(:)),min(data_IPF.HRM.AGC(:))]);

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.AGC,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.AGC,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('AGC correction (agc_ku_l1b_echo_sar_ku): IPF & DeDop','interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[dB]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'AGC_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.AGC,5,...
    data_DeDop.HRM.AGC); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dB]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.AGC,5,...
    data_IPF.HRM.AGC); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dB]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('AGC (agc_ku_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'AGC_vs_LAT_LON',file_ext)); 
close(f1)


%----------------------- Waveforms plot -----------------------------------
max_image=10*log10(max([max(data_DeDop.HRM.power_wav(:)),max(data_IPF.HRM.power_wav(:))]));
min_image=max_image-80.0;

f1=figure; 
subplot(1,2,1);
imagesc(1:data_DeDop.N_samples,data_DeDop.GEO.LAT,...
    10*log10(data_DeDop.HRM.power_wav).'); 
colormap('jet'); c=colorbar; view(2); xlabel('Samples'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
subplot(1,2,2);
imagesc(1:data_IPF.N_samples,data_IPF.GEO.LAT,...
    10*log10(data_IPF.HRM.power_wav).'); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('IPF')); 
[axT,hT]=suplabel('Waveforms (i2q2_meas_ku_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Wvfm_vs_LAT_LON',file_ext)); 
close(f1)

max_image=10*log10(max([max(data_DeDop.HRM.power_wav(:)).*128,max(data_IPF.HRM.power_wav(:))]));
min_image=max_image-40.0;

max_image=10*log10(max([mean(data_DeDop.HRM.power_wav(:)).*128,mean(data_IPF.HRM.power_wav(:))]))+25.0;
min_image=max_image-25.0;

f1=figure; 
subplot(1,2,1);
imagesc(1:data_DeDop.N_samples,data_DeDop.GEO.LAT,...
    10*log10(data_DeDop.HRM.power_wav.*128).'); 
colormap('jet'); c=colorbar; view(2); xlabel('Samples'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
subplot(1,2,2);
imagesc(1:data_IPF.N_samples,data_IPF.GEO.LAT,...
    10*log10(data_IPF.HRM.power_wav).'); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dBW]'); caxis([min_image,max_image])
title(strcat('IPF')); 
[axT,hT]=suplabel('Waveforms (i2q2_meas_ku_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Wvfm_128_corr_DeDop_vs_LAT_LON',file_ext)); 
close(f1)

idx_loc_DeDop=find(idx_DeDop);
idx_loc_IPF=find(idx_IPF);
for i_surf=1:N_surfaces_common
    if data_DeDop.GEO.LAT(idx_loc_DeDop(i_surf))>=-20.0 && data_DeDop.GEO.LAT(idx_loc_DeDop(i_surf))<=20.0
        f1=figure;
        plot(1:data_DeDop.N_samples,data_DeDop.HRM.power_wav(:,(idx_loc_DeDop(i_surf)))/...
            (max(data_DeDop.HRM.power_wav(:,(idx_loc_DeDop(i_surf))))),'-r');
        hold on; plot(1:data_IPF.N_samples,data_IPF.HRM.power_wav(:,(idx_loc_IPF(i_surf)))/...
            (max(data_IPF.HRM.power_wav(:,(idx_loc_IPF(i_surf))))),'-b');
        %hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
        title('Power waveforms: DeDop & IPF','interpreter','none');
        xlabel('Samples'); ylabel('Norm. power (max.)');
        ylim([0,1]);
        legend('DeDop','IPF')
        print(print_file,res_fig,strcat(resultPath,'Wvfm_common_num_',num2str(i_surf,'%04d'),file_ext));
        close(f1)
    end
end


%% ---------------------- Scaling factor ----------------------------------

max_image=(max([max(data_DeDop.HRM.s0_sf(:)),max(data_IPF.HRM.s0_sf(:))]));
min_image=(min([min(data_DeDop.HRM.s0_sf(:)),min(data_IPF.HRM.s0_sf(:))]));

f1=figure; 
plot(data_DeDop.GEO.LAT,data_DeDop.HRM.s0_sf,'*r');
hold on; plot(data_IPF.GEO.LAT,data_IPF.HRM.s0_sf,'ob');
%hold on; plot(1:ISD_num_surfaces_filtered,GEO_loc_errors_wdcorr,'or');
title('sigma0 scaling factor (scale_factor_ku_l1b_echo_sar_ku): IPF & DeDop','interpreter','none');
xlabel('Latitude [deg.]'); ylabel('[dBW]');
ylim([min_image,max_image]);
legend('DeDop','IPF')
print(print_file,res_fig,strcat(resultPath,'Sigma0_scale_vs_LAT',file_ext)); 
close(f1);


f1=figure; 
subplot(1,2,1);
scatter3(data_DeDop.GEO.LON,data_DeDop.GEO.LAT,...
    data_DeDop.HRM.s0_sf,5,...
    data_DeDop.HRM.s0_sf); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dB]'); caxis([min_image,max_image])
title(strcat('DeDop')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
subplot(1,2,2);
scatter3(data_IPF.GEO.LON,data_IPF.GEO.LAT,...
    data_IPF.HRM.s0_sf,5,...
    data_IPF.HRM.s0_sf); 
colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]');
ylabel(c,'[dB]'); caxis([min_image,max_image])
title(strcat('IPF')); 
hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
[axT,hT]=suplabel('sigma0 scaling factor (scale_factor_ku_l1b_echo_sar_ku)','t');
print(print_file,res_fig,strcat(resultPath,'Sigma0_scale_vs_LAT_LON',file_ext)); 
close(f1)


disp('End')


% figure; scatter3(data.GEO.LON,data.GEO.LAT,data.GEO.H/1e3,5,data.GEO.H/1e3); 
% colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]'); 
% title(strcat('Satellite altitude (Km)')); 
% hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
% 
% figure; scatter3(data.GEO.LON,data.GEO.LAT,data.MEA.win_delay.*cst_p.c_cst/2/1e3,5,data.MEA.win_delay.*cst_p.c_cst/2/1e3); 
% colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]'); 
% title(strcat('Window delay -Range (Km)')); 
% hold on; geoshow('landareas.shp','FaceColor', rgb('green'));
% 
% figure; scatter3(data.GEO.LON,data.GEO.LAT,data.GEO.H_rate,5,data.GEO.H_rate); 
% colormap('jet'); c=colorbar; view(2); xlabel('Longitude [deg]'); ylabel('Latitude [deg.]'); 
% title(strcat('Height Rate (m/s)')); 
% hold on; geoshow('landareas.shp','FaceColor', rgb('green'));

