function prepare_NetCDF_L2_SS_CCI(filename_L2,out,cnf_p,cst_p,chd_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code implements CODING & PACKING algorithm for L2 products into
% netCDF format for the Sea State CCI (Input data is assumed from Sentinel-3)
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filename_L2    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -out = structure of the output data for the L2 product
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       
% RESTRICTIONS: 
% Only valid for Sentinel-3 data and only when the analytical retracker
% with SWH fitting is considered 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------



ncid = netcdf.create(filename_L2,'NETCDF4');

long_name_att = 'long_name';
std_name_att = 'standard_name';
calendar_name_att='calendar';
comment_att = 'comment';
units_att = 'units';
scale_factor_att = 'scale_factor';
add_offset_att = 'add_offset';
flag_values_att='flag_values';
flag_desc_att='flag_meanings';
dimensions_key = 'Dimensions';
format_key = 'Format';
data_type_key = 'DataType';
fill_value_key='FillValue';

netcdf_v4_format = 'netcdf4';


ku_rec_dimension = netcdf.defDim(ncid,'time_20_ku',out.N_records);
%nl_dimension = netcdf.defDim(ncid,'max_multi_stack_ind',N_max_beams_stack_chd);
%ns_dimension = netcdf.defDim(ncid,'echo_sample_ind',N_samples*zp_fact_range_cnf);
space_3D_dimension = netcdf.defDim(ncid,'space_3D',3);

day_units = 'day';
seconds_units = 'seconds';
seconds_3dot125d64d1e9_units='3.125/64*1e-9 seconds';
seconds_3dot125d1024d1e9_units='3.125/1024*1e-9 seconds';
number_units = 'count';
degrees_units = 'degrees';
meters_units = 'meters';
meters_per_second_units = 'm/s';
rate_per_second_units='1/s';
dB_units = 'dB';
fft_pow_units='FFT power unit';
Hz_units = 'Hz';
T0d64_units = 'T0/64';
T0d16d64_units = 'T0/16/64';
W_per_count_units = 'Watt/#';
rad_units = 'rad';
percent_units='percent';


int8_type = 'NC_BYTE';
uint8_type = netcdf.getConstant('ubyte');
int16_type = 'NC_SHORT';
uint16_type = netcdf.getConstant('ushort');
int32_type = 'NC_INT';
uint32_type = netcdf.getConstant('uint');
int64_type = netcdf.getConstant('int64');
float_type = 'NC_FLOAT';
double_type= 'NC_DOUBLE';




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


%-------------------------- POSITION; -------------------------------------
alt_20_ku_name = 'alt_20_ku';
id_aux = netcdf.defVar(ncid,alt_20_ku_name,double_type,ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,214748364);
netcdf.putAtt(ncid,id_aux,long_name_att,'altitude of satellite at 20 Hz');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Altitude of the satellite Centre of Mass');
alt_20_ku = double((out.H_orb)); 


lat_20_ku_name = 'lat_20_ku';
id_aux = netcdf.defVar(ncid,lat_20_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,add_offset_att,0.);
netcdf.putAtt(ncid,id_aux,std_name_att,'latitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'latitude at 20 Hz (positive N, negative S) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Latitude of measurement [-90, +90]: Positive at Nord, Negative at South');
lat_20_ku = double(out.lat);

lon_20_ku_name = 'lon_20_ku';
id_aux = netcdf.defVar(ncid,lon_20_ku_name,double_type,ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,std_name_att,'longitude');
netcdf.putAtt(ncid,id_aux,long_name_att,'longitude at 20 Hz (positive E, negative W) (Ku-band)');
netcdf.putAtt(ncid,id_aux,units_att,degrees_units);
netcdf.putAtt(ncid,id_aux,comment_att,'longitude of measurement [-180, +180]: Positive at East, Negative at West');
lon_20_ku = double(out.lon); 


%--------------------------------------------------------------------------
%------------------------- MEASUREMENTS -----------------------------------
%--------------------------------------------------------------------------
%---------- Altimeter range and Corrections -------------------------------
range_20_ku_name = 'range_20_ku';
id_aux = netcdf.defVar(ncid,range_20_ku_name,double_type, ku_rec_dimension);
netcdf.putAtt(ncid,id_aux,long_name_att,'Corrected measured range at 20 Hz for Ku band');
netcdf.putAtt(ncid,id_aux,units_att,meters_units);
netcdf.putAtt(ncid,id_aux,comment_att,'Reference range corrected for USO frequency drift and internal path correction');
range_20_ku=double(out.range);



% %--------------------------------------------------------------------------
% %---------------------------- SCALINGS ------------------------------------
% %--------------------------------------------------------------------------
% s0_scale_factor_20_ku_name = 's0_scale_factor_20_ku';
% id_aux = netcdf.defVar(ncid,s0_scale_factor_20_ku_name,double_type, ku_rec_dimension);
% netcdf.putAtt(ncid,id_aux,long_name_att,'Scaling factor for sigma0 evaluation at 20 Hz');
% netcdf.putAtt(ncid,id_aux,units_att,dB_units);
% netcdf.putAtt(ncid,id_aux,comment_att,'This is a scaling factor in order to retrieve sigma-0 from Pu derived by retracker. It includes antenna gains and geometry satellite - surface.');
% s0_scale_factor_20_ku=double(out.s0_sf);






%--------------------------------------------------------------------------
%---------------------- RETRACKERS RESULTS --------------------------------
%--------------------------------------------------------------------------
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            %---------- ANALYTICAL RETRACKER  -----------------------------------------
            switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                    case 'SWH'
                        
                        %--------------------------------------------------------------------------
                        %----------------------------- FLAGS --------------------------------------
                        %--------------------------------------------------------------------------
                        Flag_validity_L1B_wvfm_20_ku_name = 'Flag_quality_20_ku';
                        id_aux = netcdf.defVar(ncid,Flag_validity_L1B_wvfm_20_ku_name,int8_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,127);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Flag to asses quality of retrievals at 20 Hz');
                        netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1');
                        netcdf.putAtt(ncid,id_aux,flag_desc_att,'0: bad retrieval; 1: good or valid retrieval');
                        netcdf.putAtt(ncid,id_aux,comment_att,strcat('The criteria to set Flag to 0 is: 1) the input L1B waveform is all set to zeros or NaN, or 2) Misfit (root mean square error) is above ',num2str(cnf_p.quality_flag_misfit_th)));
                        Flag_quality_20_ku=int8(out.RETRACKER.ANALYTICAL(i_index_analytical).Flag_quality);
                        
                        retracked_range_analytical_20_ku_name = 'retracked_range_20_ku';
                        id_aux = netcdf.defVar(ncid,retracked_range_analytical_20_ku_name,double_type, ku_rec_dimension);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range at 20 Hz for Ku band');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Corrected reference range ("range_20_ku") by the retracker offset, the reference range "range_20_ku" includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections). No geo-corrections applied.');
                        retracked_range_analytical_SWH_MSSfixed_20_ku=double((out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range));
                        
                       
                        swh_analytical_20_ku_name = 'SWH_20_ku';
                        id_aux = netcdf.defVar(ncid,swh_analytical_20_ku_name,double_type, ku_rec_dimension);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Significant waveheight at 20 Hz for Ku band');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,char(strcat('Fitted significant wave height.')));
                        swh_analytical_SWH_MSSfixed_20_ku=double(out.RETRACKER.ANALYTICAL(i_index_analytical).Hs);
                                                                      
                        
                        ssh_analytical_20_ku_name = 'SSH_20_ku';
                        id_aux = netcdf.defVar(ncid,ssh_analytical_20_ku_name,double_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Surface height at 20 Hz for Ku band (no geo-corrections applied)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range). No geo-corrections applied.');
                        ssh_analytical_SWH_MSSfixed_20_ku=double(out.RETRACKER.ANALYTICAL(i_index_analytical).SSH);
                        
                        sig0_analytical_20_ku_name = 'sig0_20_ku';
                        id_aux = netcdf.defVar(ncid,sig0_analytical_20_ku_name,double_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,32767);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient at 20 Hz for Ku band');
                        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
                        sig0_analytical_SWH_MSSfixed_20_ku=double(out.RETRACKER.ANALYTICAL(i_index_analytical).sigma0);
                                                
                        
                        Pearson_corr_analytical_20_ku_name = 'Pearson_corr_20_ku';
                        id_aux = netcdf.defVar(ncid,Pearson_corr_analytical_20_ku_name,double_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Pearson correlation coefficient at 20 Hz for Ku band');
                        netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Pearson correlation coefficient as percentage indicating the goodness of fitting between the real waveform and the fitted one.');
                        Pearson_corr_analytical_SWH_MSSfixed_20_ku=double(out.RETRACKER.ANALYTICAL(i_index_analytical).corr_coeff);
                        
                        Misfit_analytical_SWH_MSSfixed_20_ku_name = 'Misfit_20_ku';
                        id_aux = netcdf.defVar(ncid,Misfit_analytical_SWH_MSSfixed_20_ku_name,double_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Misfit for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,number_units);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Misfit between the real waveform and the fitted one.');
                        Misfit_analytical_SWH_MSSfixed_20_ku=double(out.RETRACKER.ANALYTICAL(i_index_analytical).misfit);
                                                                      
                case 'MSS'
                        
                case '2step'
                                      
            end            
                       
            
        case {'THRESHOLD'}
           
            
            
        case {'OCOG'}
           
            
    end
end


%------------------------ Surface Type ------------------------------------
if isfield(out,'surf_type_flag')
    surf_type_20_ku_name = 'surf_type_20_ku';
    id_aux = netcdf.defVar(ncid,surf_type_20_ku_name,int8_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,127);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter surface type');
    netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1,2,3');
    netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean or semi-enclosed_seas, enclosed_seas or lakes, continental_ice, land');
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    surf_type_20_ku =int8(out.surf_type_flag);
end

%------------------------ Sigma0 attenuation correction -------------------
if cnf_p.atm_att_correction_flag
    if isfield(out,'COR_sig0')
        % atmospheric attenuation correction on sigma0
        sig0_atmos_cor_20_ku_name = 'sig0_atmos_cor_20_ku';
        id_aux = netcdf.defVar(ncid,sig0_atmos_cor_20_ku_name,double_type, ku_rec_dimension);
        netcdf.putAtt(ncid,id_aux,long_name_att,'Atmospheric attenuation correction on sigma0');
        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
        netcdf.putAtt(ncid,id_aux,comment_att,'Atmospheric attenuation correction on backscatter coefficient (sigma0) obtained from the NCEP GFS model.');
        sig0_atmos_cor_20_ku=out.COR_sig0.Sig0AtmCorr;
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
    if isfield(out.GLOBAL_ATT,'DATA_FILE_INFO')
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'altimeter_sensor_name')
            netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'gnss_sensor_name')
            netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'doris_sensor_name')
            netcdf.putAtt(ncid,id_aux,'doris_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'acq_station_name')
            netcdf.putAtt(ncid,id_aux,'acq_station_name',out.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name);
        end
        %netcdf.putAtt(ncid,id_aux,'doris_sensor_name',acq_station_name);
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'first_meas_time')
            netcdf.putAtt(ncid,id_aux,'first_meas_time',out.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'last_meas_time')
            netcdf.putAtt(ncid,id_aux,'last_meas_time',out.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_level0')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_level0',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_level1b')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_level1b',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_orbit')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_orbit',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_doris_USO')
            netcdf.putAtt(ncid,id_aux,'xref_doris_USO',out.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_ltm_sar_cal1')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_sar_cal1',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_ltm_ku_cal2')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_ku_cal2',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_ltm_c_cal2')
            netcdf.putAtt(ncid,id_aux,'xref_altimeter_ltm_c_cal2',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'xref_altimeter_characterisation')
        netcdf.putAtt(ncid,id_aux,'xref_altimeter_characterisation',out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'semi_major_ellipsoid_axis')
        netcdf.putAtt(ncid,id_aux,'semi_major_ellipsoid_axis',out.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis);
        end
        if isfield(out.GLOBAL_ATT.DATA_FILE_INFO,'ellipsoid_flattening')
            netcdf.putAtt(ncid,id_aux,'ellipsoid_flattening',out.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening);
        end
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


if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_zero_padding',num2str(cnf_p.ZP));
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_window_azimuth',cnf_p.window_type_a);    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_window_range',cnf_p.window_type_r);
    switch cnf_p.window_type_a
        case {'Adaptive','Adaptive_S3'}
            %disp(strcat('Range PTR approx:',{''},'Adaptive'));
        otherwise
            netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Azimuth_PTR_approx',num2str(sqrt(1.0/(2.0*chd_p.alpha_ga_chd))));
    end
    switch cnf_p.window_type_r
        case {'Adaptive','Adaptive_S3'}
            %disp(strcat('Range PTR approx:',{''},'Adaptive'));
        otherwise
            netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Range_PTR_approx',num2str(sqrt(1.0/(2.0*chd_p.alpha_gr_chd))));
    end    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_use_zeros',num2str(cnf_p.use_zeros_cnf));
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Antenna_compensation_along_track',num2str(cnf_p.antenna_compensation_al));
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_preproc_active',num2str(cnf_p.pre_processing));
    if cnf_p.pre_processing
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_thresh_retracker_percentage',num2str(cnf_p.percent_leading_edge));
    end
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_floor_flag',num2str(cnf_p.Thn_flag));
    if cnf_p.Thn_flag
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_floor_method',num2str(cnf_p.Thn_estimation_method));
        switch lower(cnf_p.Thn_estimation_method)
            case 'external'
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_floor_external',num2str(cnf_p.external_Thn_value));
            case 'fixed_window'
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_1st_sample',num2str(cnf_p.Thn_w_first));
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_win_size',num2str(cnf_p.Thn_w_width));
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_over',(cnf_p.Thn_ML_SL_method));
            case 'adaptive'
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_thresh',num2str(cnf_p.threshold_noise));
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_max_iter',num2str(cnf_p.max_iter_noise));
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_mult_fact_thresh',num2str(cnf_p.factor_increase_noise_iter));
                netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_over',(cnf_p.Thn_ML_SL_method));
        end
    end
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_noise_over',cnf_p.Thn_ML_SL_method);
    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Indexation_method',cnf_p.looks_index_method);
    switch cnf_p.looks_index_method
        case 'Look_angle'
            netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Look_angle_indx_method',cnf_p.look_ang_method);
        case 'Doppler_freq'
            netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Doppler_freq_indx_method',cnf_p.fd_method);
    end
    %netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Roughness_fitting_active',num2str(cnf_p.rou_flag));
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_power_wfm_model',cnf_p.power_wfm_model);
    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Doppler_mask',cnf_p.Doppler_mask_cons_option);
    if strcmpi(cnf_p.Doppler_mask_cons_option,'internal')
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_Doppler_mask_type',cnf_p.Doppler_mask_cons_internal);
    end
    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_LUT_active',num2str(cnf_p.lut_flag));
    if cnf_p.lut_flag
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_LUT_ximin',num2str(cnf_p.LUT_ximin));
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_LUT_ximax',num2str(cnf_p.LUT_ximax));
        netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_LUT_step',num2str(cnf_p.LUT_step));
    end 
%     netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_initial_epoch',num2str(cnf_p.ini_Epoch));
%     netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_initial_SWH',num2str(cnf_p.ini_Hs));
%     netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_initial_Pu',num2str(cnf_p.ini_Pu));
%     netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_initial_MSS',num2str(cnf_p.rou_flag));
    
    netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_fitting_procedure',cnf_p.fitting_fun_type);
    switch cnf_p.fitting_fun_type
        case 'lsq'
            netcdf.putAtt(ncid,id_aux,'ANALYTICAL_retracker_lsq_minimization_algorithm',cnf_p.lsq_algorithm);
        case 'fmin'
    end
        
 
end

netcdf.putAtt(ncid,id_aux,'geo_corr_active',num2str(cnf_p.geo_corr_application_flag));
netcdf.putAtt(ncid,id_aux,'geo_corr_sametype_all_records',num2str(cnf_p.force_geocorr_surf_type));
if cnf_p.force_geocorr_surf_type
    netcdf.putAtt(ncid,id_aux,'geo_corr_type_surf',cnf_p.product_type_surface);
end
%based on the acutal correction we flag the output to account for the fact
%that even when flag is set for those missing maps of correcion of sigma0
%we will not be able to flag it back
if any(out.COR_sig0.Sig0AtmCorr)
    atm_att_correction_flag=1;
else
    atm_att_correction_flag=0;
end
netcdf.putAtt(ncid,id_aux,'atm_att_correction_flag',num2str(atm_att_correction_flag));
netcdf.putAtt(ncid,id_aux,'sigma0_bias',out.bias_sigma0(1));

netcdf.endDef(ncid);
%% --------------------- PAKCING L2 ---------------------------------------
% ----------------------- TIME/POSITION -----------------------------------
% TIME
var_id=netcdf.inqVarID(ncid,'time_20_ku');
netcdf.putVar(ncid,var_id,time_20_ku);


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
% %---------------------------- SCALINGS ------------------------------------
% %--------------------------------------------------------------------------
% var_id=netcdf.inqVarID(ncid,'s0_scale_factor_20_ku');
% netcdf.putVar(ncid,var_id,s0_scale_factor_20_ku);


%--------------------------------------------------------------------------
%----------------------------- FLAGS --------------------------------------
%--------------------------------------------------------------------------
var_id=netcdf.inqVarID(ncid,'Flag_quality_20_ku');
netcdf.putVar(ncid,var_id,Flag_quality_20_ku);

%--------------------------------------------------------------------------
%---------------------- RETRACKERS RESULTS --------------------------------
%--------------------------------------------------------------------------
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            %---------- ANALYTICAL RETRACKER  -----------------------------------------
            switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                case 'SWH'
                    var_id=netcdf.inqVarID(ncid,'retracked_range_20_ku');
                    netcdf.putVar(ncid,var_id,retracked_range_analytical_SWH_MSSfixed_20_ku);
                    
                    
                    var_id=netcdf.inqVarID(ncid,'SWH_20_ku');
                    netcdf.putVar(ncid,var_id,swh_analytical_SWH_MSSfixed_20_ku);
                    
                    
                    var_id=netcdf.inqVarID(ncid,'SSH_20_ku');
                    netcdf.putVar(ncid,var_id,ssh_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'sig0_20_ku');
                    netcdf.putVar(ncid,var_id,sig0_analytical_SWH_MSSfixed_20_ku);
                    
                    
                    var_id=netcdf.inqVarID(ncid,'Pearson_corr_20_ku');
                    netcdf.putVar(ncid,var_id,Pearson_corr_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Misfit_20_ku');
                    netcdf.putVar(ncid,var_id,Misfit_analytical_SWH_MSSfixed_20_ku);
                    
                    
                case 'MSS'
                   
                case '2step'
                    
            end
         
            
        case 'OCOG'
            
            
        case 'THRESHOLD'

            
        otherwise
            error(strcat(char(cnf_p.retracker_name(i_retracker)),' retracker not valid or available'));
    end
end



if isfield(out,'surf_type_flag')
    var_id=netcdf.inqVarID(ncid,'surf_type_20_ku');
    netcdf.putVar(ncid,var_id,surf_type_20_ku);
end

if cnf_p.atm_att_correction_flag
    if isfield(out,'COR_sig0')
        %------------- atmospheric attenuation correction ---------------------
        var_id=netcdf.inqVarID(ncid,'sig0_atmos_cor_20_ku');
        netcdf.putVar(ncid,var_id,sig0_atmos_cor_20_ku);
    end
end

netcdf.close(ncid);



end

