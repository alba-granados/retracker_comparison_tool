function prepare_NetCDF_L2(filename_L2,out,cnf_p,cst_p,chd_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code implements CODING & PACKING algorithm for L2 products into
% netCDF format
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 20/06/2016
% This software is built with internal funding 
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
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: No tracking changes until v2.0
% v2.0: Include explicitly the option to run analytical
% retracker in a 2step using cnf_p.analytical_type_of_fitting; include the
% different options for noise and wvfm selection portion in the global
% attributes


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




%% --------------------- CODING L2 ----------------------------------------
switch cnf_p.mission    
    %some issues with the L1B data set to zeroes lat and long 
    case {'CS2','CR2','S3'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case {'ESA','DeDop'}
                idx_errors=(out.TAI.days==0);        
        end
        mean_offset_height_range = 700000.0; %to save range and altitude 
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
%             netcdf.putAtt(ncid,id_aux,long_name_att,'UTC Seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)');
            netcdf.putAtt(ncid,id_aux,long_name_att,'Time in TAI: it contains the seconds since 2000-01-01 00:00:00.0+00:00 (Ku-band)'); % Alba
            netcdf.putAtt(ncid,id_aux,calendar_name_att,'Gregorian');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
			netcdf.putAtt(ncid,id_aux,comment_att,'time at surface of the SAR measurement(multilooked waveform).');
time_20_ku=double(out.TAI.total);            
switch cnf_p.mission    
    %some issues with the L1B data set to zeroes lat and long 
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                idx_int=find(~idx_errors,2,'last');
                diff_time=time_20_ku(idx_int(2))-time_20_ku(idx_int(1));
                time_20_ku(idx_errors)=time_20_ku(idx_int(2))+(1:length(find(idx_errors))).*diff_time;
        end
end



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
switch cnf_p.mission    
    %some issues with the L1B data set to zeroes lat and long 
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                alt_20_ku(idx_errors)=int32(0);        
        end
end


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
switch cnf_p.mission    
    %some issues with the L1B data set to zeroes lat and long 
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                range_20_ku(idx_errors)=int32(0);        
        end
end


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
switch cnf_p.mission    
    %some issues with the L1B data set to zeroes lat and long 
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                s0_scale_factor_20_ku(idx_errors)=int16(0);        
        end
end


%--------------------------------------------------------------------------
%----------------------------- FLAGS --------------------------------------
%--------------------------------------------------------------------------
Flag_validity_L1B_wvfm_20_ku_name = 'Flag_validity_L1B_wvfm_20_ku';
id_aux = netcdf.defVar(ncid,Flag_validity_L1B_wvfm_20_ku_name,int8_type, ku_rec_dimension);
%             netcdf.defVarFill(ncid,id_aux,false,127);
netcdf.putAtt(ncid,id_aux,long_name_att,'L1B waveform validity flag');
netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1');
netcdf.putAtt(ncid,id_aux,flag_desc_att,'0: L1B waveform not valid; 1: L1B waveform valid');
netcdf.putAtt(ncid,id_aux,comment_att,'Flag indicating whether the L1B waveform is valid or not to be used in the retracking process. One of the reasons of being not valid might be related to the noise filtering (outer looks) in L1B, such that no beam is contributing to the final stack.');
Flag_validity_L1B_wvfm_20_ku=int8(out.flag_L1B);


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
                        retracked_range_analytical_20_ku_name = 'retracked_range_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,retracked_range_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Corrected range by the retracker offset, the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
                        retracked_range_analytical_SWH_MSSfixed_20_ku=int32((out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range-mean_offset_height_range).*1e4);
                        
                        epoch_analytical_20_ku_name = 'epoch_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,epoch_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the analytical retracker (fitting SWH).');
                        epoch_analytical_SWH_MSSfixed_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor.*2/cst_p.c_cst.*1e15);
                        
                        swh_analytical_20_ku_name = 'swh_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,swh_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Significant waveheight for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,char(strcat('Fitted significant waveheight using the analytical retracker (MSS fixed to',{' '},num2str(cnf_p.rou,'%6.2e'),'). This corresponds to 4 times the fitted standard deviation of the surface height.')));
                        swh_analytical_SWH_MSSfixed_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).Hs.*1e3);
                                                                      
                        
                        ssh_analytical_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_analytical_SWH_MSSfixed_20_ku');
                        id_aux = netcdf.defVar(ncid,ssh_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Surface height for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
                        ssh_analytical_SWH_MSSfixed_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).SSH.*1e4);
                        
                        sig0_analytical_20_ku_name = 'sig0_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,sig0_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,32767);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
                        sig0_analytical_SWH_MSSfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).sigma0.*1e2);
                        
                        Pu_analytical_20_ku_name = 'Pu_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Pu_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the analytical retracker.');
                        Pu_analytical_SWH_MSSfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).Pu.*1e2);
                        
                        Pearson_corr_analytical_20_ku_name = 'Pearson_corr_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Pearson_corr_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Pearson coefficient for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Pearson correlation coefficient as percentage indicating the goodness of fitting between the real waveform and the fitted one.');
                        Pearson_corr_analytical_SWH_MSSfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).corr_coeff.*1e2);
                        
                        Misfit_analytical_SWH_MSSfixed_20_ku_name = 'Misfit_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Misfit_analytical_SWH_MSSfixed_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Misfit for Ku band (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Misfit between the real waveform and the fitted one as percentage.');
                        Misfit_analytical_SWH_MSSfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).misfit.*1e2);
                        
                        
                        Flag_fitting_analytical_20_ku_name = 'Flag_fitting_analytical_SWH_MSSfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Flag_fitting_analytical_20_ku_name,int8_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,127);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Exit flag fitting (analytical retracker fitting SWH)');
                        netcdf.putAtt(ncid,id_aux,flag_values_att,'1,2,3,4,0,-1,-2');
                        netcdf.putAtt(ncid,id_aux,flag_desc_att,'1: Function converged to a solution x; 2: Change in x was less than the specified tolerance; 3: Change in the residual was less than the specified tolerance; 4: Magnitude of search direction was smaller than the specified tolerance; 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations; -1: Output function terminated the algorithm; -2: Problem is infeasible: the bounds lb and ub are inconsistent');
                        netcdf.putAtt(ncid,id_aux,comment_att,'Flag on the reason the LSE solver stopped.');
                        Flag_fitting_analytical_SWH_MSSfixed_20_ku=int8(out.RETRACKER.ANALYTICAL(i_index_analytical).flag_fitting);
                    case 'MSS'
                        retracked_range_analytical_20_ku_name = 'retracked_range_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,retracked_range_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Corrected range by the retracker offset, the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
                        retracked_range_analytical_MSS_SWHfixed_20_ku=int32((out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range-mean_offset_height_range).*1e4);
                        
                        epoch_analytical_20_ku_name = 'epoch_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,epoch_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the analytical retracker (fitting SWH). This corresponds to zero-padded sample value.');
                        epoch_analytical_MSS_SWHfixed_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor.*2/cst_p.c_cst.*1e15);
                        
                        
                        MSS_analytical_20_ku_name = 'MSS_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,MSS_analytical_20_ku_name,int64_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Mean-square slope (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,'-');
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,char(strcat('Fitted Mean-Square Slope (SWH fixed to ',{' '},num2str(cnf_p.Hs,'%6.2e'),'): parameter modifying the RCS of surface sigma(theta)=sigma0*exp(-tan^2(theta)/MSS).')));
                        MSS_analytical_MSS_SWHfixed_20_ku=int64(out.RETRACKER.ANALYTICAL(i_index_analytical).rou.*1e15);
                        
                        ssh_analytical_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_analytical_MSS_SWHfixed_20_ku');
                        id_aux = netcdf.defVar(ncid,ssh_analytical_20_ku_name,int32_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Surface height for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-5);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
                        ssh_analytical_MSS_SWHfixed_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).SSH.*1e5);
                        
                        sig0_analytical_20_ku_name = 'sig0_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,sig0_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,32767);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
                        sig0_analytical_MSS_SWHfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).sigma0.*1e2);
                        
                        Pu_analytical_20_ku_name = 'Pu_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Pu_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the analytical retracker.');
                        Pu_analytical_MSS_SWHfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).Pu.*1e2);
                        
                        Pearson_corr_analytical_20_ku_name = 'Pearson_corr_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Pearson_corr_analytical_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Pearson coefficient for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Pearson correlation coefficient as percentage indicating the goodness of fitting between the real waveform and the fitted one.');
                        Pearson_corr_analytical_MSS_SWHfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).corr_coeff.*1e2);
                        
                        Misfit_analytical_MSS_SWHfixed_20_ku_name = 'Misfit_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Misfit_analytical_MSS_SWHfixed_20_ku_name,int16_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Misfit for Ku band (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                        netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                        netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                        netcdf.putAtt(ncid,id_aux,comment_att,'Misfit between the real waveform and the fitted one as percentage.');
                        Misfit_analytical_MSS_SWHfixed_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).misfit.*1e2);
                        
                        Flag_fitting_analytical_20_ku_name = 'Flag_fitting_analytical_MSS_SWHfixed_20_ku';
                        id_aux = netcdf.defVar(ncid,Flag_fitting_analytical_20_ku_name,int8_type, ku_rec_dimension);
                        %             netcdf.defVarFill(ncid,id_aux,false,127);
                        netcdf.putAtt(ncid,id_aux,long_name_att,'Exit flag fitting (analytical retracker fitting MSS)');
                        netcdf.putAtt(ncid,id_aux,flag_values_att,'1,2,3,4,0,-1,-2');
                        netcdf.putAtt(ncid,id_aux,flag_desc_att,'1: Function converged to a solution x; 2: Change in x was less than the specified tolerance; 3: Change in the residual was less than the specified tolerance; 4: Magnitude of search direction was smaller than the specified tolerance; 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations; -1: Output function terminated the algorithm; -2: Problem is infeasible: the bounds lb and ub are inconsistent');
                        netcdf.putAtt(ncid,id_aux,comment_att,'Flag on the reason the LSE solver stopped.');
                        Flag_fitting_analytical_MSS_SWHfixed_20_ku=int8(out.RETRACKER.ANALYTICAL(i_index_analytical).flag_fitting);
                case '2step'
                    retracked_range_analytical_20_ku_name = 'retracked_range_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,retracked_range_analytical_20_ku_name,int32_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Corrected range by the retracker offset, the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
                    retracked_range_analytical_2step_20_ku=int32((out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range-mean_offset_height_range).*1e4);
                    
                    epoch_analytical_20_ku_name = 'epoch_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,epoch_analytical_20_ku_name,int32_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the 2-step analytical retracker. This corresponds to zero-padded sample value.');
                    epoch_analytical_2step_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor.*2/cst_p.c_cst.*1e15);
                    
                    swh_analytical_20_ku_name = 'swh_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,swh_analytical_20_ku_name,int16_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Significant waveheight for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-3);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Fitted significant waveheight. This corresponds to 4 times the fitted standard deviation of the surface height.');
                    swh_analytical_2step_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).Hs.*1e3);
                    
                    MSS_analytical_20_ku_name = 'MSS_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,MSS_analytical_20_ku_name,int64_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Mean-square slope (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,'-');
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Fitted Mean-Square Slope (variance of surface slope): measure of the parameter modifying the RCS of surface sigma(theta)=sigma0*exp(-tan^2(theta)/MSS).');
                    MSS_analytical_2step_20_ku=int64(out.RETRACKER.ANALYTICAL(i_index_analytical).rou.*1e15);
                    
                    
                    ssh_analytical_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_analytical_2step_20_ku');
                    id_aux = netcdf.defVar(ncid,ssh_analytical_20_ku_name,int32_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Surface height for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
                    ssh_analytical_2step_20_ku=int32(out.RETRACKER.ANALYTICAL(i_index_analytical).SSH.*1e4);
                    
                    sig0_analytical_20_ku_name = 'sig0_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,sig0_analytical_20_ku_name,int16_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,32767);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
                    sig0_analytical_2step_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).sigma0.*1e2);
                    
                    Pu_analytical_20_ku_name = 'Pu_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,Pu_analytical_20_ku_name,int16_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform.');
                    Pu_analytical_2step_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).Pu.*1e2);
                    
                    Pearson_corr_analytical_20_ku_name = 'Pearson_corr_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,Pearson_corr_analytical_20_ku_name,int16_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Pearson coefficient for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Pearson correlation coefficient as percentage indicating the goodness of fitting between the real waveform and the fitted one.');
                    Pearson_corr_analytical_2step_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).corr_coeff.*1e2);
                    
                    Misfit_analytical_2step_20_ku_name = 'Misfit_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,Misfit_analytical_2step_20_ku_name,int16_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Misfit for Ku band (2-step analytical retracker)');
                    netcdf.putAtt(ncid,id_aux,units_att,percent_units);
                    netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
                    netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
                    netcdf.putAtt(ncid,id_aux,comment_att,'Misfit between the real waveform and the fitted one as percentage.');
                    Misfit_analytical_2step_20_ku=int16(out.RETRACKER.ANALYTICAL(i_index_analytical).misfit.*1e2);
                    
                    Flag_fitting_analytical_20_ku_name = 'Flag_fitting_analytical_2step_20_ku';
                    id_aux = netcdf.defVar(ncid,Flag_fitting_analytical_20_ku_name,int8_type, ku_rec_dimension);
                    %             netcdf.defVarFill(ncid,id_aux,false,127);
                    netcdf.putAtt(ncid,id_aux,long_name_att,'Exit flag fitting (2-step analytical fitting)');
                    netcdf.putAtt(ncid,id_aux,flag_values_att,'1,2,3,4,0,-1,-2');
                    netcdf.putAtt(ncid,id_aux,flag_desc_att,'1: Function converged to a solution x; 2: Change in x was less than the specified tolerance; 3: Change in the residual was less than the specified tolerance; 4: Magnitude of search direction was smaller than the specified tolerance; 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations; -1: Output function terminated the algorithm; -2: Problem is infeasible: the bounds lb and ub are inconsistent');
                    netcdf.putAtt(ncid,id_aux,comment_att,'Flag on the reason the LSE solver stopped (for the best fitting either SWH or MSS).');
                    Flag_fitting_analytical_2step_20_ku=int8(out.RETRACKER.ANALYTICAL(i_index_analytical).flag_fitting);                    
            end            
                       
            
        case {'THRESHOLD'}
            %---------- ANALYTICAL RETRACKER  -----------------------------------------
            retracked_range_threshold_20_ku_name = 'retracked_range_threshold_20_ku';
            id_aux = netcdf.defVar(ncid,retracked_range_threshold_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range for Ku band (threshold retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
            netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
            netcdf.putAtt(ncid,id_aux,comment_att,'Corrected range by the retracker offset (using analytical retracker), the reference range includes instrumental corrections (already the USO frequency drift and the internal/instrument corrections).');
            retracked_range_threshold_20_ku=int32((out.RETRACKER.THRESHOLD.tracker_range-mean_offset_height_range).*1e4);
            switch cnf_p.mission
                %some issues with the L1B data set to zeroes lat and long
                case {'CS2','CR2'}
                    %% --------------------- CroySAT-2 ------------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            retracked_range_threshold_20_ku(idx_errors)=int32(0);
                    end
            end
            
            
            epoch_threshold_20_ku_name = 'epoch_threshold_20_ku';
            id_aux = netcdf.defVar(ncid,epoch_threshold_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for threshold retracker and Ku band');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the threshold retracker. This corresponds to zero-padded sample value.');
            epoch_threshold_20_ku=int32(out.RETRACKER.THRESHOLD.retracking_cor.*2/cst_p.c_cst.*1e15);
            switch cnf_p.mission
                %some issues with the L1B data set to zeroes lat and long
                case {'CS2','CR2'}
                    %% --------------------- CroySAT-2 ------------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            epoch_threshold_20_ku(idx_errors)=int32(0);
                    end
            end
            
            
            
            ssh_threshold_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_threshold_20_ku');
            id_aux = netcdf.defVar(ncid,ssh_threshold_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Sea surface height for Ku band (threshold retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Sea surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
            ssh_threshold_20_ku=int32(out.RETRACKER.THRESHOLD.SSH.*1e4);
            switch cnf_p.mission
                %some issues with the L1B data set to zeroes lat and long
                case {'CS2','CR2'}
                    %% --------------------- CroySAT-2 ------------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            ssh_threshold_20_ku(idx_errors)=int32(0);
                    end
            end
            
            
            
            sig0_threshold_20_ku_name = 'sig0_threshold_20_ku';
            id_aux = netcdf.defVar(ncid,sig0_threshold_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (threshold retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
            sig0_threshold_20_ku=int16(out.RETRACKER.THRESHOLD.sigma0.*1e2);
            switch cnf_p.mission
                %some issues with the L1B data set to zeroes lat and long
                case {'CS2','CR2'}
                    %% --------------------- CroySAT-2 ------------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            sig0_threshold_20_ku(idx_errors)=int32(0);
                    end
            end
            
            
            
            Pu_threshold_20_ku_name = 'Pu_threshold_20_ku';
            id_aux = netcdf.defVar(ncid,Pu_threshold_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (threshold retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the threshold retracker.');
            Pu_threshold_20_ku=int16(out.RETRACKER.THRESHOLD.Pu.*1e2);  
            switch cnf_p.mission
                %some issues with the L1B data set to zeroes lat and long
                case {'CS2','CR2'}
                    %% --------------------- CroySAT-2 ------------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            Pu_threshold_20_ku(idx_errors)=int32(0);
                    end
            end
            
            
        case {'OCOG'}
            %---------- ANALYTICAL RETRACKER  -----------------------------------------
            retracked_range_OCOG_20_ku_name = 'retracked_range_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,retracked_range_OCOG_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Retracked range for Ku band (OCOG retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
            netcdf.putAtt(ncid,id_aux,add_offset_att,mean_offset_height_range);
            netcdf.putAtt(ncid,id_aux,comment_att,'Corrected range by the retracker offset (using OCOG retracker), the reference range includes already the USO frequency drift and the internal/instrument corrections');
            retracked_range_OCOG_20_ku=int32((out.RETRACKER.OCOG.tracker_range-mean_offset_height_range).*1e4);
            
            epoch_OCOG_20_ku_name = 'epoch_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,epoch_OCOG_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Epoch for Ku band (OCOG retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Estimated epoch in seconds w.r.t center of the window (window delay is given to the center of the window) using the OCOG retracker. This corresponds to zero-padded sample value.');
            epoch_OCOG_20_ku=int32(out.RETRACKER.OCOG.retracking_cor.*2/cst_p.c_cst.*1e15);
            
            
            cog_OCOG_20_ku_name = 'cog_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,cog_OCOG_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Center of gravity for OCOG retracker');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Estimated center of gravity in seconds using the OCOG retracker. This corresponds to zero-padded sample value.');
            cog_OCOG_20_ku=int32(out.RETRACKER.OCOG.COG.*2/cst_p.c_cst.*1e15);
            
            width_OCOG_20_ku_name = 'width_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,width_OCOG_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Window width for OCOG retracker and Ku band');
            netcdf.putAtt(ncid,id_aux,units_att,seconds_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-15);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Estimated window width in seconds using the OCOG retracker. This corresponds to zero-padded sample value.');
            width_OCOG_20_ku=int32(out.RETRACKER.OCOG.COG.*2/cst_p.c_cst.*1e15);
            
            
            ssh_OCOG_20_ku_name = strcat(cnf_p.nc_name_surface_height,'_OCOG_20_ku');
            id_aux = netcdf.defVar(ncid,ssh_OCOG_20_ku_name,int32_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Sea surface height for Ku band (OCOG ice retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,meters_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-4);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Sea surface heigth above the elliposid of reference and extracted using the orbital height and the corrected range (retracked range).');
            ssh_OCOG_20_ku=int32(out.RETRACKER.OCOG.SSH.*1e4);
            
            sig0_OCOG_20_ku_name = 'sig0_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,sig0_OCOG_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,32767);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Backscattering coefficient for Ku band (OCOG retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Backscattering coefficient extracted as the fitted peak power once corrected by the sigma0 scaling factor.');
            sig0_OCOG_20_ku=int16(out.RETRACKER.OCOG.sigma0.*1e2);
            
            Pu_OCOG_20_ku_name = 'Pu_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,Pu_OCOG_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Fitted peak power for Ku band (OCOG retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the OCOG retracker.');
            Pu_OCOG_20_ku=int16(out.RETRACKER.OCOG.Pu.*1e2);
            
            A_OCOG_20_ku_name = 'A_OCOG_20_ku';
            id_aux = netcdf.defVar(ncid,A_OCOG_20_ku_name,int16_type, ku_rec_dimension);
            %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
            netcdf.putAtt(ncid,id_aux,long_name_att,'Amplitude OCOG for Ku band (OCOG retracker)');
            netcdf.putAtt(ncid,id_aux,units_att,dB_units);
            netcdf.putAtt(ncid,id_aux,scale_factor_att,1.e-2);
            netcdf.putAtt(ncid,id_aux,add_offset_att,0.0);
            netcdf.putAtt(ncid,id_aux,comment_att,'Peak power of the fitted waveform using the threshold retracker.');
            A_OCOG_20_ku=int16(out.RETRACKER.OCOG.A.*1e2);
            
    end
end



%--------------------------------------------------------------------------
% -----------------------GEOPHYSICAL CORRECTIONS --------------------------
%--------------------------------------------------------------------------
if isfield(out,'COR')
    dry_tropo_correction_name = 'dry_tropo_correction_20_ku';
    id_aux = netcdf.defVar(ncid,dry_tropo_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Dry Tropospheric Correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    wet_tropo_correction_name = 'wet_tropo_correction_20_ku';
    id_aux = netcdf.defVar(ncid,wet_tropo_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Wet Tropospheric correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    inverse_baro_correction_name = 'inverse_baro_correction_20_ku';
    id_aux = netcdf.defVar(ncid,inverse_baro_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Inverse Barometric Correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    Dynamic_atmospheric_correction_name = 'Dynamic_atmospheric_correction_20_ku';
    id_aux = netcdf.defVar(ncid,Dynamic_atmospheric_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Dynamic Atmospheric Correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    GIM_iono_correction_name = 'GIM_iono_correction_20_ku';
    id_aux = netcdf.defVar(ncid,GIM_iono_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'GIM Ionospheric Correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    model_iono_correction_name = 'model_iono_correction_20_ku';
    id_aux = netcdf.defVar(ncid,model_iono_correction_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Model Ionospheric Correction');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    ocean_equilibrium_tide_name = 'ocean_equilibrium_tide_20_ku';
    id_aux = netcdf.defVar(ncid,ocean_equilibrium_tide_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Equilibrium Tide');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    long_period_tide_height_name = 'long_period_tide_20_ku';
    id_aux = netcdf.defVar(ncid,long_period_tide_height_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Long Period Ocean Tide');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    ocean_loading_tide_name = 'ocean_loading_tide_20_ku';
    id_aux = netcdf.defVar(ncid,ocean_loading_tide_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Ocean Loading Tide');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    solid_earth_tide_name = 'solid_earth_tide_20_ku';
    id_aux = netcdf.defVar(ncid,solid_earth_tide_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Solid Earth Tide');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    geocentric_polar_tide_name = 'geocentric_polar_tide_20_ku';
    id_aux = netcdf.defVar(ncid,geocentric_polar_tide_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Geocentric Polar Tide');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    
    total_geo_corr_name = 'total_geo_corr_20_ku';
    id_aux = netcdf.defVar(ncid,total_geo_corr_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Total geophysical corrections');
    netcdf.putAtt(ncid,id_aux,units_att,meters_units);
    netcdf.putAtt(ncid,id_aux,scale_factor_att,1e-3);
    netcdf.putAtt(ncid,id_aux,comment_att,'It corresponds to the aggregation of the different geophysical corrections contributing, depending on the surface being imaged, either assuming/forcing same surface for all records or using the surface type flag information.');
    
    dry_tropo_correction_20_ku=int32(out.COR.dry_trop.*1e3);
    wet_tropo_correction_20_ku=int32(out.COR.wet_trop.*1e3);
    inverse_baro_correction_20_ku=int32(out.COR.inv_bar.*1e3);
    Dynamic_atmospheric_correction_20_ku=int32(out.COR.dac.*1e3);
    GIM_iono_correction_20_ku=int32(out.COR.gim_ion.*1e3);
    model_iono_correction_20_ku=int32(out.COR.model_ion.*1e3);
    ocean_equilibrium_tide_20_ku=int32(out.COR.ocean_equilibrium_tide.*1e3);
    long_period_tide_20_ku=int32(out.COR.ocean_longperiod_tide.*1e3);
    ocean_loading_tide_20_ku=int32(out.COR.ocean_loading_tide.*1e3);
    solid_earth_tide_20_ku=int32(out.COR.solidearth_tide.*1e3);
    geocentric_polar_tide_20_ku=int32(out.COR.geocentric_polar_tide.*1e3);
    total_geo_corr_20_ku=int32(out.COR.total_geo_corr.*1e3);
    
    %-------------------- FLAGS ON GEOPHYSICAL CORRECTIONS --------------------
    geo_corr_flags=zeros(1,out.N_records); %using a int16
    bit11to15='00000'; %last five bits are reserved
    for i_record=1:out.N_records
        bit0=dec2bin(out.COR.flags.dry_trop(i_record));
        bit1=dec2bin(out.COR.flags.wet_trop(i_record));
        bit2=dec2bin(out.COR.flags.inv_bar(i_record));
        bit3=dec2bin(out.COR.flags.dac(i_record));
        bit4=dec2bin(out.COR.flags.gim_ion(i_record));
        bit5=dec2bin(out.COR.flags.model_ion(i_record));
        bit6=dec2bin(out.COR.flags.ocean_equilibrium_tide(i_record));
        bit7=dec2bin(out.COR.flags.ocean_longperiod_tide(i_record));
        bit8=dec2bin(out.COR.flags.ocean_loading_tide(i_record));
        bit9=dec2bin(out.COR.flags.solidearth_tide(i_record));
        bit10=dec2bin(out.COR.flags.geocentric_polar_tide(i_record));
        geo_corr_flags(i_record)=bin2dec([bit0 bit1 bit2 bit3 bit4 bit5 bit6 bit7 bit8 bit9 bit10 bit11to15]);
    end
    flags_geo_corr_20_ku=int32(geo_corr_flags);
    clear geo_corr_flags;
    
    geo_corr_flags_name = 'flags_geo_corr_20_ku';
    id_aux = netcdf.defVar(ncid,geo_corr_flags_name,int32_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,2147483647);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Geophysical corrections application');
    netcdf.putAtt(ncid,id_aux,flag_values_att,'0 or 1');
    netcdf.putAtt(ncid,id_aux,flag_desc_att,'Geophysical correction applied (1) or not applied (0)');
    netcdf.putAtt(ncid,id_aux,comment_att,'The int32 value needs to be converted to binary and the first 11 bits correspond to the flags of each individual correction: (1) dry tropo, (2) wet tropo, (3) inverse barometric, (4) dynamic atmospheric correction, (5) GIM ionospheric correction, (6) Model ionospheric correction , (7) Ocean tide, (8) Ocean long period tide, (9) Ocean loading tide, (10) Solid earth tide and (11) Geocentric polar tide.');
    
   
end
if isfield(out,'COR_sig0')
    % atmospheric attenuation correction on sigma0
    sig0_atmos_cor_20_ku_name = 'sig0_atmos_cor_20_ku';
    id_aux = netcdf.defVar(ncid,sig0_atmos_cor_20_ku_name,double_type, ku_rec_dimension);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Atmospheric attenuation correction on sigma0');
    netcdf.putAtt(ncid,id_aux,units_att,dB_units);
    netcdf.putAtt(ncid,id_aux,comment_att,'Atmospheric attenuation correction on backscatter coefficient (sigma0) obtained from the NCEP GFS model.');
    sig0_atmos_cor_20_ku=out.COR_sig0.Sig0AtmCorr;
end


%------------------------ Surface Type ------------------------------------
if isfield(out,'surf_type_flag')
    surf_type_20_ku_name = 'surf_type_20_ku';
    id_aux = netcdf.defVar(ncid,surf_type_20_ku_name,int8_type, ku_rec_dimension);
    %             netcdf.defVarFill(ncid,id_aux,false,127);
    netcdf.putAtt(ncid,id_aux,long_name_att,'Altimeter surface type');
    netcdf.putAtt(ncid,id_aux,flag_values_att,'0,1,2,3');
    netcdf.putAtt(ncid,id_aux,flag_desc_att,'open_ocean or semi-enclosed_seas, enclosed_seas or lakes, continental_ice, land,Transponder');
    netcdf.putAtt(ncid,id_aux,comment_att,'Value the closest in time to the reference measurement: for the burst just above the surface');
    surf_type_20_ku =int8(out.surf_type_flag);
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
        case {'CS2','CR2'}
            netcdf.putAtt(ncid,id_aux,'altimeter_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name);            
            netcdf.putAtt(ncid,id_aux,'gnss_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name);
            netcdf.putAtt(ncid,id_aux,'doris_sensor_name',out.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name);
            netcdf.putAtt(ncid,id_aux,'acq_station_name',out.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name);
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
            %--------------- add the attributes related to intermediate product--------
            netcdf.putAtt(ncid,id_aux,'orbit_phase_code',out.GLOBAL_ATT.ORBIT_INFO.orbit_phase_code);
            netcdf.putAtt(ncid,id_aux,'orbit_cycle_num',out.GLOBAL_ATT.ORBIT_INFO.orbit_cycle_num);
            netcdf.putAtt(ncid,id_aux,'orbit_REL_Orbit',out.GLOBAL_ATT.ORBIT_INFO.orbit_REL_Orbit);
            netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Start',out.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Start);
            netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Start',out.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Start);
            netcdf.putAtt(ncid,id_aux,'orbit_ABS_Orbit_Stop',out.GLOBAL_ATT.ORBIT_INFO.orbit_ABS_Orbit_Stop);
            netcdf.putAtt(ncid,id_aux,'orbit_Rel_Time_ASC_Node_Stop',out.GLOBAL_ATT.ORBIT_INFO.orbit_Rel_Time_ASC_Node_Stop);
            netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Time',out.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Time);
            netcdf.putAtt(ncid,id_aux,'orbit_Equator_Cross_Long',out.GLOBAL_ATT.ORBIT_INFO.orbit_Equator_Cross_Long);
            netcdf.putAtt(ncid,id_aux,'orbit_Ascending_Flag',out.GLOBAL_ATT.ORBIT_INFO.orbit_Ascending_Flag);
            netcdf.putAtt(ncid,id_aux,'Start_Lat',out.GLOBAL_ATT.DATA_FILE_INFO.Start_Lat);
            netcdf.putAtt(ncid,id_aux,'Start_Long',out.GLOBAL_ATT.DATA_FILE_INFO.Start_Long);
            netcdf.putAtt(ncid,id_aux,'Stop_Lat',out.GLOBAL_ATT.DATA_FILE_INFO.Stop_Lat);
            netcdf.putAtt(ncid,id_aux,'Stop_Long',out.GLOBAL_ATT.DATA_FILE_INFO.Stop_Long);
            
        case {'S3','S3A','S3B'}
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
            netcdf.putAtt(ncid,id_aux,'cycle_number',out.GLOBAL_ATT.DATA_FILE_INFO.cycle_number);
            netcdf.putAtt(ncid,id_aux,'pass_number',out.GLOBAL_ATT.DATA_FILE_INFO.pass_number);
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
if any(strcmp(cnf_p.retracker_name,'OCOG'))
   netcdf.putAtt(ncid,id_aux,'OCOG_retracker_percentage_power',num2str(cnf_p.OCOG_retracker.percentage_pow_OCOG));
   netcdf.putAtt(ncid,id_aux,'OCOG_retracker_first_zp_sample',num2str(cnf_p.OCOG_retracker.n1));
   netcdf.putAtt(ncid,id_aux,'OCOG_retracker_last_zp_sample',num2str(cnf_p.OCOG_retracker.n2));
   netcdf.putAtt(ncid,id_aux,'OCOG_retracker_offset',num2str(cnf_p.OCOG_retracker.offset));   
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
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            %---------- ANALYTICAL RETRACKER  -----------------------------------------
            switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                case 'SWH'
                    var_id=netcdf.inqVarID(ncid,'retracked_range_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,retracked_range_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'epoch_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,epoch_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'swh_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,swh_analytical_SWH_MSSfixed_20_ku);
                    
                    
                    var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_analytical_SWH_MSSfixed_20_ku'));
                    netcdf.putVar(ncid,var_id,ssh_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'sig0_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,sig0_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pu_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Pu_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pearson_corr_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Pearson_corr_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Misfit_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Misfit_analytical_SWH_MSSfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Flag_fitting_analytical_SWH_MSSfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Flag_fitting_analytical_SWH_MSSfixed_20_ku);
                    
                case 'MSS'
                    var_id=netcdf.inqVarID(ncid,'retracked_range_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,retracked_range_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'epoch_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,epoch_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'MSS_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,MSS_analytical_MSS_SWHfixed_20_ku);
                    
                    
                    var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_analytical_MSS_SWHfixed_20_ku'));
                    netcdf.putVar(ncid,var_id,ssh_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'sig0_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,sig0_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pu_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Pu_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pearson_corr_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Pearson_corr_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Misfit_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Misfit_analytical_MSS_SWHfixed_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Flag_fitting_analytical_MSS_SWHfixed_20_ku');
                    netcdf.putVar(ncid,var_id,Flag_fitting_analytical_MSS_SWHfixed_20_ku);
                case '2step'
                    var_id=netcdf.inqVarID(ncid,'retracked_range_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,retracked_range_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'epoch_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,epoch_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'swh_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,swh_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'MSS_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,MSS_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_analytical_2step_20_ku'));
                    netcdf.putVar(ncid,var_id,ssh_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'sig0_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,sig0_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pu_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,Pu_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Pearson_corr_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,Pearson_corr_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Misfit_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,Misfit_analytical_2step_20_ku);
                    
                    var_id=netcdf.inqVarID(ncid,'Flag_fitting_analytical_2step_20_ku');
                    netcdf.putVar(ncid,var_id,Flag_fitting_analytical_2step_20_ku);
            end
         
            
        case 'OCOG'
            %----------------- OCOG_ICE RETRACKER -------------------------
            var_id=netcdf.inqVarID(ncid,'retracked_range_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,retracked_range_OCOG_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'epoch_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,epoch_OCOG_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'cog_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,cog_OCOG_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'width_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,width_OCOG_20_ku);
                                  
            var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_OCOG_20_ku'));
            netcdf.putVar(ncid,var_id,ssh_OCOG_20_ku);
                        
            var_id=netcdf.inqVarID(ncid,'sig0_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,sig0_OCOG_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'Pu_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,Pu_OCOG_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'A_OCOG_20_ku');
            netcdf.putVar(ncid,var_id,A_OCOG_20_ku);
            
        case 'THRESHOLD'
            %---------- THRESHOLD RETRACKER  -----------------------------------------
            var_id=netcdf.inqVarID(ncid,'retracked_range_threshold_20_ku');
            netcdf.putVar(ncid,var_id,retracked_range_threshold_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'epoch_threshold_20_ku');
            netcdf.putVar(ncid,var_id,epoch_threshold_20_ku);
                                  
            var_id=netcdf.inqVarID(ncid,strcat(cnf_p.nc_name_surface_height,'_threshold_20_ku'));
            netcdf.putVar(ncid,var_id,ssh_threshold_20_ku);
                        
            var_id=netcdf.inqVarID(ncid,'sig0_threshold_20_ku');
            netcdf.putVar(ncid,var_id,sig0_threshold_20_ku);
            
            var_id=netcdf.inqVarID(ncid,'Pu_threshold_20_ku');
            netcdf.putVar(ncid,var_id,Pu_threshold_20_ku);
            
        otherwise
            error(strcat(char(cnf_p.retracker_name(i_retracker)),' retracker not valid or available'));
    end
end

% -----------------------GEOPHYSICAL CORRECTIONS --------------------------
if isfield(out,'COR')
    var_id=netcdf.inqVarID(ncid,'dry_tropo_correction_20_ku');
    netcdf.putVar(ncid,var_id,dry_tropo_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'wet_tropo_correction_20_ku');
    netcdf.putVar(ncid,var_id,wet_tropo_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'inverse_baro_correction_20_ku');
    netcdf.putVar(ncid,var_id,inverse_baro_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'Dynamic_atmospheric_correction_20_ku');
    netcdf.putVar(ncid,var_id,Dynamic_atmospheric_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'GIM_iono_correction_20_ku');
    netcdf.putVar(ncid,var_id,GIM_iono_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'model_iono_correction_20_ku');
    netcdf.putVar(ncid,var_id,model_iono_correction_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'ocean_equilibrium_tide_20_ku');
    netcdf.putVar(ncid,var_id,ocean_equilibrium_tide_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'long_period_tide_20_ku');
    netcdf.putVar(ncid,var_id,long_period_tide_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'ocean_loading_tide_20_ku');
    netcdf.putVar(ncid,var_id,ocean_loading_tide_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'solid_earth_tide_20_ku');
    netcdf.putVar(ncid,var_id,solid_earth_tide_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'geocentric_polar_tide_20_ku');
    netcdf.putVar(ncid,var_id,geocentric_polar_tide_20_ku);
    
    var_id=netcdf.inqVarID(ncid,'total_geo_corr_20_ku');
    netcdf.putVar(ncid,var_id,total_geo_corr_20_ku);
    
    %---------------------- flags of geophysical correction -------------------
    var_id=netcdf.inqVarID(ncid,'flags_geo_corr_20_ku');
    netcdf.putVar(ncid,var_id,flags_geo_corr_20_ku);
       
    
end

if isfield(out,'COR_sig0')
     %------------- atmospheric attenuation correction ---------------------
    var_id=netcdf.inqVarID(ncid,'sig0_atmos_cor_20_ku');
    netcdf.putVar(ncid,var_id,sig0_atmos_cor_20_ku);
end


if isfield(out,'surf_type_flag')
    var_id=netcdf.inqVarID(ncid,'surf_type_20_ku');
    netcdf.putVar(ncid,var_id,surf_type_20_ku);
end

netcdf.close(ncid);



end

