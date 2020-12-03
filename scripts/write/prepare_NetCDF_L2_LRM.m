function prepare_NetCDF_L2_LRM(filename_L2,out,cnf_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code implements CODING  algorithm for L2 LR products into
% netCDF format
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados - Eduard Makhoul / isardSAT
%
% Reviewer:         ----- / isardSAT
% 
% Last revision:    Alba Granados / isardSAT V1 07/09/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filename_L2    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       
% RESTRICTIONS/COMMENTS:
% Addapted from /retracker/write/prepare_NetCDF_L2.m by Eduard Makhoul
% 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:



ncid = netcdf.create(filename_L2,'NETCDF4');

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
netcdf_v4_format = 'netcdf4';


ku_rec_dimension = netcdf.defDim(ncid,'time_20_ku',out.N_records);
%nl_dimension = netcdf.defDim(ncid,'max_multi_stack_ind',N_max_beams_stack_chd);
%ns_dimension = netcdf.defDim(ncid,'echo_sample_ind',N_samples*zp_fact_range_cnf);

day_units = 'day';
seconds_units = 'seconds';
degrees_units = 'degrees';
meters_units = 'meters';
dB_units = 'dB';
percent_units='percent';

int16_type = 'NC_SHORT';
int32_type = 'NC_INT';
double_type= 'NC_DOUBLE';




%% --------------------- CODING L2 ----------------------------------------
switch cnf_p.mission    
    case {'S6'}
        mean_offset_height_range = 1300000.0;
end
%--------------------------------------------------------------------------
% ----------------------- TIME/POSITION -----------------------------------
%--------------------------------------------------------------------------

% -------------------------------- TIME -----------------------------------
time_20_ku_name = 'time_20_ku';
id_aux      = netcdf.defVar(ncid,time_20_ku_name,double_type,ku_rec_dimension);
            netcdf.putAtt(ncid,id_aux,std_name_att,'time');
            netcdf.putAtt(ncid,id_aux,long_name_att,'UTC Seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'time at surface of the SAR measurement(multilooked waveform).');
time_20_ku=double(out.TAI.total);            

UTC_day_l1b_20_ku_name = 'UTC_day_20_ku';
id_aux      = netcdf.defVar(ncid,UTC_day_l1b_20_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Days since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,units_att,day_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'days elapsed since 2000-01-01. To be used to link with L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');
UTC_day_20_ku=int16(out.TAI.days);     


UTC_sec_20_ku_name = 'UTC_sec_20_ku';
id_aux = netcdf.defVar(ncid,UTC_sec_20_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,1.844674407370960e+19);
netcdf.putAtt(ncid,id_aux,long_name_att,'Seconds in the day UTC, with microsecond resolution (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
netcdf.putAtt(ncid,id_aux,comment_att,'seconds in the day. To be used to link L1 and L2 records (time_l1b provides the number of seconds since 2000-01-01).');
UTC_sec_20_ku=double(out.TAI.secs+out.TAI.microsecs.*1e-6); 

%-------------------------- POSITION; -------------------------------------
alt_20_ku_name = 'alt_20_ku';
id_aux = netcdf.defVar(ncid,alt_20_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');
alt_20_ku = int32((out.H_orb-mean_offset_height_range).*1e4); 

lat_20_ku_name = 'lat_20_ku';
id_aux = netcdf.defVar(ncid,lat_20_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude (positive N, negative S) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');
lat_20_ku = int32(out.lat.*1e6);

lon_20_ku_name = 'lon_20_ku';
id_aux = netcdf.defVar(ncid,lon_20_ku_name,int32_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-6);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,comment_att,'longitude of measurement [-180, +180]: Positive at East, Negative at West');
lon_20_ku = int32(out.lon.*1e6); 


%--------------------------------------------------------------------------
%------------------------- MEASUREMENTS -----------------------------------
%--------------------------------------------------------------------------
%---------- Altimeter range and Corrections -------------------------------
range_20_ku_name = 'range_20_ku';
id_aux = netcdf.defVar(ncid,range_20_ku_name,int32_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,2147483647);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected measured range for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');
range_20_ku=int32((out.range-mean_offset_height_range).*1e4);

%--------------------------------------------------------------------------
%---------------------------- SCALINGS ------------------------------------
%--------------------------------------------------------------------------
s0_scale_factor_20_ku_name = 's0_scale_factor_20_ku';
id_aux = netcdf.defVar(ncid,s0_scale_factor_20_ku_name,int16_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,32767);
netcdf.putAtt(ncid,id_aux,long_name_att,'Scaling factor for sigma0 evaluation');
netcdf.putAtt(ncid,id_aux,units_att,dB_units);
netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
netcdf.putAtt(ncid,id_aux,comment_att,'This is a scaling factor in order to retrieve sigma-0 from Pu derived by retracker. It includes antenna gains and geometry satellite - surface.');
s0_scale_factor_20_ku=int16(out.s0_sf.*1e2);
% 


%--------------------------------------------------------------------------
%---------------------- RETRACKERS RESULTS --------------------------------
%--------------------------------------------------------------------------
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'LR'}
            epoch_LR_20_ku_name = 'epoch_LR_20_ku';
            id_aux = netcdf.defVar(ncid,epoch_LR_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for Ku band (LR retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the LR retracker. This corresponds to zero-padded sample value.');
%             epoch_LR_20_ku=int32(out.RETRACKER.LR.retracking_cor.*2/cst_p.c_cst.*1e15);
            epoch_LR_20_ku=int32(out.RETRACKER.LR.Epoch.*1e15);
            
            cnf_p.nc_name_surface_height = 'ssh';
            ssh_LR_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_LR_20_ku');
            id_aux = netcdf.defVar(ncid,ssh_LR_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Sea surface height for Ku band (LR ice retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Sea surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
            ssh_LR_20_ku=int32(out.RETRACKER.LR.SSH.*1e4);

            swh_LR_20_ku_name = 'swh_LR_20_ku';
            id_aux = netcdf.defVar(ncid,swh_LR_20_ku_name,int16_type, ku_rec_dimension);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Significant wave height for Ku band (LR retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
            swh_LR_20_ku=int16(out.RETRACKER.LR.Hs.*1e3);
            
            sig0_LR_20_ku_name = 'sig0_LR_20_ku';
            id_aux = netcdf.defVar(ncid,sig0_LR_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (LR retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
            sig0_LR_20_ku=int16(out.RETRACKER.LR.sigma0.*1e2);
            
            Pu_LR_20_ku_name = 'Pu_LR_20_ku';
            id_aux = netcdf.defVar(ncid,Pu_LR_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (LR retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the LR retracker.');
            Pu_LR_20_ku=int16(out.RETRACKER.LR.Pu.*1e2);

            Pearson_corr_analytical_20_ku_name = 'Pearson_corr_Pu_LR_20_ku';
            id_aux = netcdf.defVar(ncid,Pearson_corr_analytical_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Pearson coefficient for Ku band (LR retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,percent_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Pearson correlation coefficient as percentage indicating the goodness of fitting between the real waveform and the fitted one.');
            Pearson_corr_Pu_LR_20_ku=int16(out.RETRACKER.LR.COR.*1e2);
            
            
    end
end


% ----------------------- GLOBAL ATTRIBUTES -------------------------------
%----------  Global Attributes definition ---------------------------------
%---- attributes inherited from Sentinel-3 product description-------------
id_aux = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,id_aux,'creation_time',char(out.date_creation));
netcdf.putAtt(ncid,id_aux,'Conventions',netcdf_v4_format);
netcdf.putAtt(ncid,id_aux,'mission_name',cnf_p.mission);
netcdf.putAtt(ncid,id_aux,'operation_mode',cnf_p.mode);
if isfield(out,'GLOBAL_ATT')
    switch cnf_p.mission
%         case {'S3','S3A','S3B','S6'}
        case {'S6'}
            netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name);
            netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name);
            netcdf.putAtt(ncid,id_aux,'doris_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name);
            netcdf.putAtt(ncid,id_aux,'acq_station_name',out.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name);
            %netcdf.putAtt(ncid,id_aux,'doris_sensor_name',acq_station_name);
            netcdf.putAtt(ncid,id_aux,'first_meas_time',out.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time);
            netcdf.putAtt(ncid,id_aux,'last_meas_time',out.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_level0',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_level1b',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_orbit',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit);
            netcdf.putAtt(ncid,id_aux,'xref_doris_USO',out.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_sar_cal1',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_ku_cal2',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_c_cal2',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2);
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_characterisation',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation);
            netcdf.putAtt(ncid,id_aux,'semi_major_ellipsoid_axis',out.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis);
            netcdf.putAtt(ncid,id_aux,'ellipsoid_flattening',out.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening);
 
            netcdf.putAtt(ncid,id_aux,'pass_number',out.GLOBAL_ATT.DATA_FILE_INFO.pass_number);
            netcdf.putAtt(ncid,id_aux,'cycle_number',out.GLOBAL_ATT.DATA_FILE_INFO.cycle_number);
        otherwise
            fprintf('Mission not contemplated to prepare L2 output netCDF file.\n');
    end
 
end
% -------------- Processing configuration attributes ----------------------
netcdf.putAtt(ncid,id_aux,'L1B_processor',cnf_p.L1proc);
netcdf.putAtt(ncid,id_aux,'Retrackers',char(strjoin(cnf_p.retracker_name,', ')));
if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
    netcdf.putAtt(ncid,id_aux,'Analytical_Retrackers_fitting',char(strjoin(cnf_p.analytical_type_of_fitting,', ')));
end
%netcdf.putAtt(ncid,id_aux,'Seed_information_exploitation',num2str(cnf_p.seed));
netcdf.putAtt(ncid,id_aux,'Geographical_masking_active',num2str(cnf_p.mask_ROI_flag));
if cnf_p.mask_ROI_flag
    netcdf.putAtt(ncid,id_aux,'Geographical_mask_file_KML',out.GLOBAL_ATT.DATA_FILE_INFO.geographical_mask_kml)
end
netcdf.putAtt(ncid,id_aux,'Filtering land surfaces',num2str(cnf_p.filter_land));
if cnf_p.filter_land
    netcdf.putAtt(ncid,id_aux,'Filtering land surfaces type',(cnf_p.filter_land_type));
end
netcdf.putAtt(ncid,id_aux,'Looks_masking_flag',num2str(cnf_p.mask_looks_flag));
if cnf_p.mask_looks_flag
    netcdf.putAtt(ncid,id_aux,'Looks_masking_Nlooks',num2str(cnf_p.Neff_thres));
end

netcdf.putAtt(ncid,id_aux,'Discard_waveform_samples',num2str(cnf_p.wvfm_discard_samples));
if cnf_p.wvfm_discard_samples
    netcdf.putAtt(ncid,id_aux,'Number_of_samples_discarded_from_begining',num2str(cnf_p.wvfm_discard_samples_begin));
    netcdf.putAtt(ncid,id_aux,'Number_of_samples_discarded_from_end',num2str(cnf_p.wvfm_discard_samples_end));
end

netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_flag',num2str(cnf_p.wvfm_portion_selec));
if cnf_p.wvfm_portion_selec
    netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_type',cnf_p.wvfm_portion_selec_type);
    switch lower(cnf_p.wvfm_portion_selec_type)
        case 'ref_height'
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_DEM',cnf_p.wvfm_portion_selec_DEM_ref);            
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_left_samples_ref_pos',num2str(cnf_p.wvfm_portion_selec_l_samples));
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_right_samples_ref_pos',num2str(cnf_p.wvfm_portion_selec_r_samples));
        case 'peak_win'            
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_left_samples_ref_pos',num2str(cnf_p.wvfm_portion_selec_l_samples));
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_right_samples_ref_pos',num2str(cnf_p.wvfm_portion_selec_r_samples));
        case 'peak_thresh'            
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_left_thresh_ref_pos',num2str(cnf_p.wvfm_portion_selec_l_thres));
            netcdf.putAtt(ncid,id_aux,'Wvfm_portion_selection_right_thresh_ref_pos',num2str(cnf_p.wvfm_portion_selec_r_thres));
    end
end

if any(strcmp(cnf_p.retracker_name,'THRESHOLD'))
   netcdf.putAtt(ncid,id_aux,'THRESHOLD_retracker_percentage_peak',num2str(cnf_p.th_retracker.percentage_peak));
end

netcdf.putAtt(ncid,id_aux,'geo_corr_active',num2str(cnf_p.geo_corr_application_flag));
netcdf.putAtt(ncid,id_aux,'geo_corr_sametype_all_records',num2str(cnf_p.force_geocorr_surf_type));
if cnf_p.force_geocorr_surf_type
    netcdf.putAtt(ncid,id_aux,'geo_corr_type_surf',cnf_p.product_type_surface);
end

netcdf.endDef(ncid);
%% --------------------- PAKCING L2 ---------------------------------------
% ----------------------- TIME/POSITION -----------------------------------
% TIME
var_id=netcdf.inqVarID(ncid,'time_20_ku');
netcdf.putVar(ncid,var_id,time_20_ku);

var_id=netcdf.inqVarID(ncid,'UTC_day_20_ku');
netcdf.putVar(ncid,var_id,UTC_day_20_ku);

var_id=netcdf.inqVarID(ncid,'UTC_sec_20_ku');
netcdf.putVar(ncid,var_id,UTC_sec_20_ku);

% POSITION
var_id=netcdf.inqVarID(ncid,'alt_20_ku');
netcdf.putVar(ncid,var_id,alt_20_ku);

var_id=netcdf.inqVarID(ncid,'lat_20_ku');
netcdf.putVar(ncid,var_id,lat_20_ku);

var_id=netcdf.inqVarID(ncid,'lon_20_ku');
netcdf.putVar(ncid,var_id,lon_20_ku);



%--------------------------------------------------------------------------
%------------------------- MEASUREMENTS -----------------------------------
%--------------------------------------------------------------------------
% range
var_id=netcdf.inqVarID(ncid,'range_20_ku');
netcdf.putVar(ncid,var_id,range_20_ku);


%--------------------------------------------------------------------------
%---------------------------- SCALINGS ------------------------------------
%--------------------------------------------------------------------------
var_id=netcdf.inqVarID(ncid,'s0_scale_factor_20_ku');
netcdf.putVar(ncid,var_id,s0_scale_factor_20_ku);


% %--------------------------------------------------------------------------
% %----------------------------- FLAGS --------------------------------------
% %--------------------------------------------------------------------------
% var_id=netcdf.inqVarID(ncid,'Flag_validity_L1B_wvfm_20_ku');
% netcdf.putVar(ncid,var_id,Flag_validity_L1B_wvfm_20_ku);

%--------------------------------------------------------------------------
%---------------------- RETRACKERS RESULTS --------------------------------
%--------------------------------------------------------------------------
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))            
        case 'LR'
            var_id=netcdf.inqVarID(ncid,'epoch_LR_20_ku');
            netcdf.putVar(ncid,var_id,epoch_LR_20_ku);
            
            var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_LR_20_ku'));
            netcdf.putVar(ncid,var_id,ssh_LR_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'swh_LR_20_ku');
            netcdf.putVar(ncid,var_id,swh_LR_20_ku);
                        
            var_id=netcdf.inqVarID(ncid,'sig0_LR_20_ku');
            netcdf.putVar(ncid,var_id,sig0_LR_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'Pu_LR_20_ku');
            netcdf.putVar(ncid,var_id,Pu_LR_20_ku);

            var_id=netcdf.inqVarID(ncid,'Pearson_corr_Pu_LR_20_ku');
            netcdf.putVar(ncid,var_id,Pearson_corr_Pu_LR_20_ku);
            
    end
end


netcdf.close(ncid);



end

