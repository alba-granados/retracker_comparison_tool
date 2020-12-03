function data = readL1B_S6_ISD (filename_L1B,cnf_p,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from L1B Sentinel-6 products
% processed by ISD
%
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 15/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       MANDATORY:
%           -filename_L1B    =   L1B filename with the fullpath information
%           -cnf_p           =   Structure with configuration parameters of the
%                                L2 processor
%       OPTIONAL:
%           -filename_L1BS   =   L1B-S filename with the fullpath name
%       
% 
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Based on read_alt_data.m
% v1.1: Include global attributes sensor name, acq. station, acq date,... for
% the reference L1B product & orbital parameter info in a similar fashion
% as the global attribtues provided for output L1B product according to
% Sentinel-3 format and within the SEOMs projects. It is considered only
% for the final netCDF L1B products as defined in the PSF issue 1.3 
%(not all fields considered in SEOMS Sentinel-3 are available in the Sentinel-6 product format netcdf)
% v2.0: Updated to read new format of netcdf from Sentinel-6 based on
% groups

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(3+1*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
p.parse(varargin{:});
filename_L1BS=char(p.Results.filename_L1BS);
clear p;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
[~,name_file,ext]=fileparts(filename_L1B);
ext=lower(ext);
name='HR';
switch ext
    case {'.nc','.NC'}
        %% --------------------- NetCDF ---------------------------------------
	
        % Based on the L1B product generated for DeDop a la Sentinel-3
        % open the netCDF file
        ncid=netcdf.open(filename_L1B,'NC_NOWRITE');
        
        % Get dimensions
%         dimid=netcdf.inqDimID(ncid,'data_20/ku/time');
%         [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
        num_surfaces = length(ncread(filename_L1B,'data_20/ku/time_tai'));

%         dimid=netcdf.inqDimID(ncid,'ns');
        [~,filename_L1B_nopath,~]=fileparts(filename_L1B);
        name_split = split(filename_L1B_nopath,'_');
        if any(strcmp(name_split, 'range'))
            dimid=netcdf.inqDimID(ncid,'samples_no');
        elseif any(strcmp(name_split, 'HR'))
            dimid=netcdf.inqDimID(ncid,'samples_ov');
        else 
            dimid=netcdf.inqDimID(ncid,'samples'); 
        end
        [~,N_samples]=netcdf.inqDim(ncid,dimid);
        
		switch cnf_p.mode
			case {'RAW','RMC','HR'}
				        dimid=netcdf.inqDimID(ncid,'nl');
			case {'LR-RMC'}
			        dimid=netcdf.inqDimID(ncid,'nb');
		end

        [~,N_looks_conf_stack]=netcdf.inqDim(ncid,dimid); %total number of looks configured for a stack
        netcdf.close(ncid);
        
        switch cnf_p.mode
            case {'RAW', 'SAR'}
                %RMC mode
                data.N_samples=N_samples;
            case {'RMC','LR-RMC'}
                %RMC mode
                %data.N_samples=N_samples/2; %(only half of them selected to fit waveform)
                data.N_samples=N_samples;
        end
        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------

        data.GEO.TAI.total             =   ncread(filename_L1B,'data_20/ku/time_tai').';
        data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
        data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
        data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
        data.GEO.LAT                   =   ncread(filename_L1B,'data_20/ku/latitude').';
        data.GEO.LON                   =   wrapTo180(ncread(filename_L1B,'data_20/ku/longitude').'); %defined between [0,360]        
        data.GEO.H_rate                =   ncread(filename_L1B,'data_20/ku/altitude_rate').'; % m/s
        velocity=double(ncread(filename_L1B,'data_20/ku/velocity_vector')).'; %vector of vx,vx,vz
        data.GEO.V                     =   sqrt(velocity(:,1).^2+velocity(:,2).^2+velocity(:,3).^2).';
        clear velocity;
        data.GEO.H                     =   ncread(filename_L1B,'data_20/ku/altitude').';
        data.GEO.pitch                 =   (ncread(filename_L1B,'data_20/ku/off_nadir_pitch_angle_pf').').*pi/180;
        data.GEO.roll                  =   (ncread(filename_L1B,'data_20/ku/off_nadir_roll_angle_pf').').*pi/180;
        data.GEO.yaw                   =   (ncread(filename_L1B,'data_20/ku/off_nadir_yaw_angle_pf').').*pi/180;
        
        
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = double(ncread(filename_L1B,'data_20/ku/tracker_range_calibrated')).'*2.0/cst_p.c_cst;
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
%         % Currently is not foreseen to be included in the L1B-product of
%         % Sentinel-6
%         data.COR.surf_type_flag=ncread(filename_L1B,'flag_surface_classification').'; %{0: open ocean, 1: land, 2: continental_water, 3: acquatic vegetation, 4: continental_ice_snow, 5: floating_ice, 6: salted_basin}
%         

        %---------------------------------------------------------------
        % SWH/SIGMA0 instrumental correction
        %---------------------------------------------------------------
        data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
        data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
        
        % -----------------------------------------------------------------
        % HRM: High-resolution mode: Waveforms
        % -----------------------------------------------------------------
        % ------------------ Waveforms ------------------------------------
        %if ~isempty(strfind(name,'HR'))            
            switch cnf_p.mode 
                case {'SAR','RAW'}
                    i2q2_meas=double(ncread(filename_L1B,'data_20/ku/power_waveform')).';
                    scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
                    %apply scaling factor to waveforms
                    data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
                    clear i2q2_meas scale_factor;
                case {'RMC'}
                    %consider only the half of the waveform for fitting
                    i2q2_meas=double(ncread(filename_L1B,'data_20/ku/hr_power_waveform')).';
                    scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
                    %apply scaling factor to waveforms
                    data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
                    clear i2q2_meas scale_factor;
                    %data.HRM.power_wav=data.HRM.power_wav(1:data.N_samples,:); 
                case {'LR-RMC'}
                    i2q2_meas=double(ncread(filename_L1B,'data_20/ku/lr_rmc_power_waveform')).';
                    scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
                    %apply scaling factor to waveforms
                    data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
                    clear i2q2_meas scale_factor;					
            end
            
            
            switch cnf_p.mode 
                case {'SAR','RAW','RMC','HR'}
					data.HRM.Neff       =   ncread(filename_L1B,'data_20/ku/num_looks_start_stop').'; %effective number of beams that form the stack including possuible looks that are set entirely to zero
				case {'LR-RMC'}
					data.HRM.Neff       =   64*ones(1,data.N_records); %a single waveform look will be considered
            end
            data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
            data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
            %             data.HRM.ThN        =   zeros(1,num_surfaces);
            %             data.HRM.wfm_count  =   1:1:num_surfaces;
            
            
            % ----sigma0 scaling factor for conversion from Power to sigma0
            %units
            data.HRM.s0_sf=ncread(filename_L1B,'data_20/ku/sig0_scaling_factor').';
            
            
            %--------------------------------------------------------------
            %------------- Stack characterization parameters --------------
            %--------------------------------------------------------------
            %----- Dopppler mask ------------------------------------------
			switch cnf_p.mode 
                case {'SAR','RAW','RMC','HR'}
					data.HRM.Doppler_mask   =   (ncread(filename_L1B,'data_20/ku/stack_mask_start_stop')+1); %due to the way the stack mask vector is saved in the netcdf the mask is not exactly the same as original
				case {'LR-RMC'}
					data.HRM.Doppler_mask   =   (ncread(filename_L1B,'data_20/ku/burst_mask')+1); 
			end
            % internally when using the ncread function the scaling factor
            % of ZP is considered
            data.HRM.pri_surf=ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'; %in seconds
            
            %clock from L1B
            data.HRM.fs_clock_ku_surf=ncread(filename_L1B,'data_20/ku/altimeter_clock').'; %in Hz
            
            %------------- Geometry-related parameters ------------------------
            if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
				switch cnf_p.mode
					case {'HR','SAR','RAW','RMC'}
						switch cnf_p.looks_index_method
							case {'Doppler'}
								%This option Does not provide the correct information but it is kept for comparison purposes
								switch cnf_p.fd_method
									case 'exact'
										% exact method would require precise information for
										% each beam within the stack (not currently developed for L1B products for SEOMs)
										% Dopp_ang= 90-Beam_ang
										data.HRM.beam_ang_stack   =  pi/2-ncread(filename_L1BS,'doppler_angle'); %TBD in radians
										% it is not foreseen to provide the specific
										% PRI/velocity for each look coming from a different
										% burst: use the same pri for that surface
										data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*...
											(ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'); %TBD %eventually for sentinel-6 the PRF could change
										data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V;
									case 'approximate'
										data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'data_20/ku/pointing_angle_start').'); % in radians
										data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'data_20/ku/pointing_angle_stop').'); % in radians
										%corresponds to the complementary of the
										%Beam angle
										data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'data_20/ku/doppler_angle_start').'); % in radians
										data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'data_20/ku/doppler_angle_stop').'); % in radians
								end
							case {'Look_angle'}
								switch cnf_p.look_ang_method
									case 'exact'
										% exact method would require precise information for
										% each beam within the stack (not currently developed for L1B products for SEOMs)
										data.HRM.look_ang_stack   =  ncread(filename_L1B,'data_20/ku/look_angle'); % in radians
										% it is not foreseen to provide the specific
										% PRI/velocity for each look coming from a different
										% burst: use the same pri for that surface
										data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*(ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'); %TBD %eventually for sentinel-6 the PRF could change
										data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V; % TBD norm of velocity for each beam within stack
									case 'approximate'
										data.HRM.look_ang_start_surf=(ncread(filename_L1B,'data_20/ku/look_angle_start').'); % in radians
										data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'data_20/ku/look_angle_stop').'); % in radians
								end
								
						end
				
				end
                
            end
            
            
            %------------------------------------------------------------------
            %-------------------GLOBAL ATTRIBUTES -----------------------------
            %------------------------------------------------------------------            
            data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=ncreadatt(filename_L1B,'/','mission_name');
            data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
            data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
            data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';            
            data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=ncreadatt(filename_L1B,'/','first_measurement_time');
            data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=ncreadatt(filename_L1B,'/','last_measurement_time');
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=ncreadatt(filename_L1B,'/','xref_altimeter_level0');
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=ncreadatt(filename_L1B,'/','xref_orbit');
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=ncreadatt(filename_L1B,'/','xref_doris_uso');
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available';
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available';
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available';
            data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=ncreadatt(filename_L1B,'/','xref_characterization');
%             data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(ncreadatt(filename_L1B,'/','semi_major_ellipsoid_axis'));
            data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(ncreadatt(filename_L1B,'/','ellipsoid_semi_major_axis'));
            data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(ncreadatt(filename_L1B,'/','ellipsoid_flattening'));
            data.GLOBAL_ATT.DATA_FILE_INFO.cycle_number = ncreadatt(filename_L1B,'/','cycle_number');
            data.GLOBAL_ATT.DATA_FILE_INFO.pass_number = ncreadatt(filename_L1B,'/','pass_number');

    case '.mat'
        %% ------------------------ Matlab data -------------------------------
        load(filename_L1B);
        switch cnf_p.mode
            case {'RAW','RMC','HR'}
                
                s=size(wfm_cor_i2q2_sar_ku);
                N_samples=s(2);
                num_surfaces=s(1);
                data.N_samples=N_samples;
                data.N_records=num_surfaces;
                
                % -----------------------------------------------------------------
                % GEO: geographical information
                % -----------------------------------------------------------------
                data.GEO.TAI.total             =   time_surf;
                data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
                data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
                data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
                data.GEO.LAT                   =   lat_surf; %degrees
                data.GEO.LON                   =   wrapTo180(lon_surf); %degrees defined between [-180,180]
                data.GEO.H_rate                =   alt_rate_sat; % m/s
                data.GEO.V                     =   sqrt(x_vel_sat.^2+y_vel_sat.^2+z_vel_sat.^2); %m/s
                data.GEO.H                     =   alt_sat; % m
                data.GEO.pitch                 =   pitch_surf; % rad
                data.GEO.roll                  =   roll_surf; % rad
                data.GEO.yaw                   =   yaw_surf; % rad
                
                % -----------------------------------------------------------------
                % MEA: measurements
                % -----------------------------------------------------------------
                data.MEA.win_delay = win_delay_surf;
                
                % ---------------------------------------------------------------------
                % COR: Geophysical corrections
                % ---------------------------------------------------------------------
                % Currently is not foreseen to be included in the L1B-product of
                % Sentinel-6
                
                
                %---------------------------------------------------------------
                % SWH/SIGMA0 instrumental correction
                %---------------------------------------------------------------
                data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
                data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
                
                
                
                % -----------------------------------------------------------------
                % HRM: High-resolution mode: Waveforms
                % -----------------------------------------------------------------
                % ------------------ Waveforms ------------------------------------
%                 if ~isempty(strfind(name,'HR'))
                    % ------------ SAR mode ---------------------------------------
                    data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
                    
                    
                    data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
                    data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
                    data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
                    %             data.HRM.ThN        =   zeros(1,num_surfaces);
                    %             data.HRM.wfm_count  =   1:1:num_surfaces;
                    
                    
                    % ----sigma0 scaling factor for conversion from Power to sigma0
                    %units
                    data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
                    
                    
                    %--------------------------------------------------------------
                    %------------- Stack characterization parameters --------------
                    %--------------------------------------------------------------
                    %----- Dopppler mask ------------------------------------------
                    data.HRM.Doppler_mask   =   stack_mask_vector_start_stop.';
                    data.HRM.pri_surf=pri_sar_surf; %in seconds
                    data.HRM.fs_clock_ku_surf=1./T0_sar_surf; %in Hz
                    
                    %------------- Geometry-related parameters ------------------------
                    if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
                        switch cnf_p.looks_index_method
                            case {'Doppler'}
                                %This option Does not provide the correct information but it is kept for comparison purposes
                                switch cnf_p.fd_method
                                    case 'exact'
                                        % exact method would require precise information for
                                        % each beam within the stack (not currently developed for L1B products for SEOMs)
                                        data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
                                        data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
                                        data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
                                    case 'approximate'
                                        data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
                                        data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
                                        data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
                                        data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
                                end
                            case {'Look_angle'}
                                switch cnf_p.look_ang_method
                                    case 'exact'
                                        % exact method would require precise information for
                                        % each beam within the stack (not currently developed for L1B products for SEOMs)
                                        data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
                                        data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
                                        data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
                                    case 'approximate'
                                        data.HRM.look_ang_start_surf=start_look_angle; % in radians
                                        data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
                                end
                                
                        end
                    end
                    
                    %------------------------------------------------------------------
                    %-------------------GLOBAL ATTRIBUTES -----------------------------
                    %------------------------------------------------------------------

                    data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
                    data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
                    data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
                    data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
                                        
                    TAI_2012 = (12*365+4+181)*3600*24 + 35;
                    TAI_2015 = (15*365+4+181)*3600*24 + 36;
                    time_UTC_init = double(min(data.GEO.TAI.total));
                    time_UTC_end  = double(max(data.GEO.TAI.total));
                    
                    % add leap seconds to the TAI time. Only valid for
                    if(time_UTC_init < TAI_2012)
                        time_UTC_init = time_UTC_init - 34;                        
                    elseif(time_UTC_init > TAI_2015)
                        time_UTC_init = time_UTC_init - 36;
                    else
                        time_UTC_init = time_UTC_init - 35;
                    end
                    if(time_UTC_end < TAI_2012)
                        time_UTC_end = time_UTC_end - 34;
                    elseif(time_UTC_end > TAI_2015)
                        time_UTC_end = time_UTC_end - 36;
                    else
                        time_UTC_end = time_UTC_end - 35;
                    end
                    
                    data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
                    data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
                    data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
                                        
                    
%                 elseif ~isempty(strfind(name,'LR'))
%                     % ----------- LR: low resolution data extraction --------------
%                 end
                
            case {'LR-RMC'}

                s=size(wfm_i2q2_sar_ku);
                N_samples=s(2);
                num_surfaces=s(1);
                data.N_samples=N_samples;
                data.N_records=num_surfaces;
                
                % -----------------------------------------------------------------
                % GEO: geographical information
                % -----------------------------------------------------------------
                data.GEO.TAI.total             =   time_aux_seconds;
                data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
                data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
                data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
                data.GEO.LAT                   =   latitude_ku; %degrees
                data.GEO.LON                   =   wrapTo180(longitude_ku); %degrees defined between [-180,180]
                data.GEO.H_rate                =   com_altitude_rate_ku; % m/s
                data.GEO.V                     =   sqrt(vel_x.^2+vel_y.^2+vel_z.^2); %m/s
                data.GEO.H                     =   com_altitude_ku; % m
                data.GEO.pitch                 =   off_nadir_pitch_angle_pf_20_ku; % rad
                data.GEO.roll                  =   off_nadir_roll_angle_pf_20_ku; % rad
                data.GEO.yaw                   =   off_nadir_yaw_angle_pf_20_ku; % rad
                
                % -----------------------------------------------------------------
                % MEA: measurements
                % -----------------------------------------------------------------
                data.MEA.win_delay = altimeter_range_calibrated_ku*2/cst_p.c_cst;
                
                % ---------------------------------------------------------------------
                % COR: Geophysical corrections
                % ---------------------------------------------------------------------
                % Currently is not foreseen to be included in the L1B-product of
                % Sentinel-6
                
                
                %---------------------------------------------------------------
                % SWH/SIGMA0 instrumental correction
                %---------------------------------------------------------------
                data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
                data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
               
               
                               
                
                % -----------------------------------------------------------------
                % HRM: High-resolution mode: Waveforms
                % -----------------------------------------------------------------
                % ------------------ Waveforms ------------------------------------
                
                % ------------ SAR mode ---------------------------------------
                data.HRM.power_wav=(double(wfm_i2q2_sar_ku).*double(repmat(wfmScalingFactor.',[1,N_samples]))).';
                %i2q2_meas.*repmat(scale_factor.',N_samples,1);
                
                
                data.HRM.Neff       =   N_beams_start_stop_combined; %effective number of beams that form the stack including possuible looks that are set entirely to zero
                data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
                data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
                %             data.HRM.ThN        =   zeros(1,num_surfaces);
                %             data.HRM.wfm_count  =   1:1:num_surfaces;
                
                
                % ----sigma0 scaling factor for conversion from Power to sigma0
                %units
                data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
                
                
                %--------------------------------------------------------------
                %------------- Stack characterization parameters --------------
                %--------------------------------------------------------------
                %----- Dopppler mask ------------------------------------------
                data.HRM.Doppler_mask   =   stack_mask_vector_start_stop_combined.';
                data.HRM.pri_surf=double(pri_lrm_l1b_ku); %in seconds
                data.HRM.fs_clock_ku_surf=1./double(altimeter_clock_ku); %in Hz
                
                %------------- Geometry-related parameters ------------------------
                if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
                    switch cnf_p.looks_index_method
                        case {'Doppler'}
                            
                        case {'Look_angle'}
                            switch cnf_p.look_ang_method
                                case 'exact'
                                    % exact method would require precise information for
                                    % each beam within the stack (not currently developed for L1B products for SEOMs)
                                    data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
                                    data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
                                    data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
                                case 'approximate'
                                    data.HRM.look_ang_start_surf=start_look_angle_combined; % in radians
                                    data.HRM.look_ang_stop_surf=stop_look_angle_combined; % in radians
                            end
                            
                    end
                end
                
                %------------------------------------------------------------------
                %-------------------GLOBAL ATTRIBUTES -----------------------------
                %------------------------------------------------------------------
                
                data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
                data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
                data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
                data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
                
                TAI_2012 = (12*365+4+181)*3600*24 + 35;
                TAI_2015 = (15*365+4+181)*3600*24 + 36;
                time_UTC_init = double(min(data.GEO.TAI.total));
                time_UTC_end  = double(max(data.GEO.TAI.total));
                
                % add leap seconds to the TAI time. Only valid for
                if(time_UTC_init < TAI_2012)
                    time_UTC_init = time_UTC_init - 34;
                elseif(time_UTC_init > TAI_2015)
                    time_UTC_init = time_UTC_init - 36;
                else
                    time_UTC_init = time_UTC_init - 35;
                end
                if(time_UTC_end < TAI_2012)
                    time_UTC_end = time_UTC_end - 34;
                elseif(time_UTC_end > TAI_2015)
                    time_UTC_end = time_UTC_end - 36;
                else
                    time_UTC_end = time_UTC_end - 35;
                end
                
                data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
                data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
                data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
                data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
                
                
          case {'FF-RAW','FF-RMC'}
                
                s=size(pow_ML_FF);
                N_samples=s(2);
                num_surfaces=s(1);
                data.N_samples=N_samples;
                data.N_records=num_surfaces;
                
                % -----------------------------------------------------------------
                % GEO: geographical information
                % -----------------------------------------------------------------
                data.GEO.TAI.total             =   time_surf_ML;
                data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
                data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
                data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
                data.GEO.LAT                   =   lat_surf_ML; %degrees
                data.GEO.LON                   =   wrapTo180(lon_surf_ML); %degrees defined between [-180,180]
                data.GEO.H_rate                =   alt_rate_sat_surf_ML; % m/s
                data.GEO.V                     =   sqrt(x_vel_sat_surf_ML.^2+y_vel_sat_surf_ML.^2+z_vel_sat_surf_ML.^2); %m/s
                data.GEO.H                     =   alt_sat_surf_ML; % m
                data.GEO.pitch                 =   pitch_sat_surf_ML; % rad
                data.GEO.roll                  =   roll_sat_surf_ML; % rad
                data.GEO.yaw                   =   yaw_sat_surf_ML; % rad
                
                % -----------------------------------------------------------------
                % MEA: measurements
                % -----------------------------------------------------------------
                data.MEA.win_delay = win_delay_surf_ML;
                
                % ---------------------------------------------------------------------
                % COR: Geophysical corrections
                % ---------------------------------------------------------------------
                % Currently is not foreseen to be included in the L1B-product of
                % Sentinel-6
                
                
                %---------------------------------------------------------------
                % SWH/SIGMA0 instrumental correction
                %---------------------------------------------------------------
                data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
                data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
                                
                
                % -----------------------------------------------------------------
                % HRM: High-resolution mode: Waveforms
                % -----------------------------------------------------------------
                % ------------------ Waveforms ------------------------------------
%                 if ~isempty(strfind(name,'HR'))
                    % ------------ SAR mode ---------------------------------------
                    data.HRM.power_wav=pow_ML_FF.';
                    
                    
                    data.HRM.Neff       =   470*ones(1,num_surfaces); %Fix the value to a maximum view on sentinel-6 from simulations
                    data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
                    data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
                    %             data.HRM.ThN        =   zeros(1,num_surfaces);
                    %             data.HRM.wfm_count  =   1:1:num_surfaces;
                    
                    
                    % ----sigma0 scaling factor for conversion from Power to sigma0
                    %units
                    data.HRM.s0_sf=zeros(1,num_surfaces); %force 0 dB uni value until able to provide a
                    
                    
                    %--------------------------------------------------------------
                    %------------- Stack characterization parameters --------------
                    %--------------------------------------------------------------
                    %----- Dopppler mask ------------------------------------------
                    %data.HRM.Doppler_mask   =   data.N_samples.*ones(1,num_surfaces);
                    data.HRM.pri_surf=pri_surf_ML; %in seconds
                    data.HRM.fs_clock_ku_surf=1./T0_surf_ML; %in Hz
                    
                    %------------- Geometry-related parameters ------------------------
                    if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
                        switch cnf_p.looks_index_method
                           
                            case {'Look_angle'}
                                switch cnf_p.look_ang_method
                                    case 'exact'
                                    case 'approximate'
                                        data.HRM.look_ang_start_surf=zeros(1,num_surfaces); % in radians
                                        data.HRM.look_ang_stop_surf=zeros(1,num_surfaces); % in radians
                                end
                                
                        end
                    end
                    
                    % FIx them based on what is observed from the HR-DDP values
                    
                    data.HRM.doppler_angle_start_surf =  0.40*pi/180*ones(1,num_surfaces);
                    data.HRM.doppler_angle_stop_surf = -0.42*pi/180*ones(1,num_surfaces);
                    
                    %------------------------------------------------------------------
                    %-------------------GLOBAL ATTRIBUTES -----------------------------
                    %------------------------------------------------------------------
                    
                    data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
                    data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
                    data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
                    data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
                                        
                    TAI_2012 = (12*365+4+181)*3600*24 + 35;
                    TAI_2015 = (15*365+4+181)*3600*24 + 36;
                    time_UTC_init = double(min(data.GEO.TAI.total));
                    time_UTC_end  = double(max(data.GEO.TAI.total));
                    
                    % add leap seconds to the TAI time. Only valid for
                    if(time_UTC_init < TAI_2012)
                        time_UTC_init = time_UTC_init - 34;                        
                    elseif(time_UTC_init > TAI_2015)
                        time_UTC_init = time_UTC_init - 36;
                    else
                        time_UTC_init = time_UTC_init - 35;
                    end
                    if(time_UTC_end < TAI_2012)
                        time_UTC_end = time_UTC_end - 34;
                    elseif(time_UTC_end > TAI_2015)
                        time_UTC_end = time_UTC_end - 36;
                    else
                        time_UTC_end = time_UTC_end - 35;
                    end
                    
                    data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
                    data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
                    data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
                    data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);      
        end              
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end

% LR-RMC for single waveforms at burst level not averaged by RC
%                 s=size(wfm_cor_i2q2_sar_ku);
%                 N_samples=s(2);
%                 num_surfaces=s(1);
%                 data.N_samples=N_samples;
%                 data.N_records=num_surfaces;
%                 
%                 % -----------------------------------------------------------------
%                 % GEO: geographical information
%                 % -----------------------------------------------------------------
%                 data.GEO.TAI.total             =   time_sar_ku;
%                 data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%                 data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%                 data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%                 data.GEO.LAT                   =   lat_sar_sat.'; %degrees
%                 data.GEO.LON                   =   wrapTo180(lon_sar_sat.'); %degrees defined between [-180,180]
%                 data.GEO.H_rate                =   alt_rate_sar_sat; % m/s
%                 data.GEO.V                     =   sqrt(x_vel_sat_sar.^2+y_vel_sat_sar.^2+z_vel_sat_sar.^2); %m/s
%                 data.GEO.H                     =   alt_sar_sat.'; % m
%                 data.GEO.pitch                 =   pitch_sar; % rad
%                 data.GEO.roll                  =   roll_sar; % rad
%                 data.GEO.yaw                   =   yaw_sar; % rad
%                 
%                 % -----------------------------------------------------------------
%                 % MEA: measurements
%                 % -----------------------------------------------------------------
%                 data.MEA.win_delay = win_delay_sar_ku;
%                 
%                 % ---------------------------------------------------------------------
%                 % COR: Geophysical corrections
%                 % ---------------------------------------------------------------------
%                 % Currently is not foreseen to be included in the L1B-product of
%                 % Sentinel-6
%                 
%                 
%                 %---------------------------------------------------------------
%                 % SWH/SIGMA0 instrumental correction
%                 %---------------------------------------------------------------
%                 data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%                 data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%                 
%                 
%                 % -----------------------------------------------------------------
%                 % HRM: High-resolution mode: Waveforms
%                 % -----------------------------------------------------------------
%                 % ------------------ Waveforms ------------------------------------
%                 
%                 % ------------ SAR mode ---------------------------------------
%                 data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
%                 
%                 
%                 data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
%                 data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%                 data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%                 %             data.HRM.ThN        =   zeros(1,num_surfaces);
%                 %             data.HRM.wfm_count  =   1:1:num_surfaces;
%                 
%                 
%                 % ----sigma0 scaling factor for conversion from Power to sigma0
%                 %units
%                 data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
%                 
%                 
%                 %--------------------------------------------------------------
%                 %------------- Stack characterization parameters --------------
%                 %--------------------------------------------------------------
%                 %----- Dopppler mask ------------------------------------------
%                 data.HRM.Doppler_mask   =   stack_mask_vector_start_stop.';
%                 data.HRM.pri_surf=pri_sar_pre_dat; %in seconds
%                 data.HRM.fs_clock_ku_surf=1./T0_sar_pre_dat; %in Hz
%                 
%                 %------------- Geometry-related parameters ------------------------
%                 if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
%                     switch cnf_p.looks_index_method
%                         case {'Doppler'}
%                             %This option Does not provide the correct information but it is kept for comparison purposes
%                             switch cnf_p.fd_method
%                                 case 'exact'
%                                     % exact method would require precise information for
%                                     % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                     data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
%                                     data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                     data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
%                                 case 'approximate'
%                                     data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
%                                     data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
%                                     data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
%                                     data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
%                             end
%                         case {'Look_angle'}
%                             switch cnf_p.look_ang_method
%                                 case 'exact'
%                                     % exact method would require precise information for
%                                     % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                     data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
%                                     data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                     data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
%                                 case 'approximate'
%                                     data.HRM.look_ang_start_surf=start_look_angle; % in radians
%                                     data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
%                             end
%                             
%                     end
%                 end




% %% ---------------- Handling input variables ------------------------------
% if(nargin<2 || nargin>(3+1*2))
%     error('Wrong number of input parameters');   
% end
% %option to include the L1B_S product to read the exact info from Look to
% %Look within the stack
% p = inputParser;
% p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
% p.parse(varargin{:});
% filename_L1BS=char(p.Results.filename_L1BS);
% clear p;
% 
% %% ------------------------------------------------------------------------- 
% % Loading data L1B
% % ------------------------------------------------------------------------- 
% [~,name_file,ext]=fileparts(filename_L1B);
% ext=lower(ext);
% name='HR';
% switch ext
%     case '.nc'
%         %% --------------------- NetCDF ---------------------------------------
% 	
%         % Based on the L1B product generated for DeDop a la Sentinel-3
%         % open the netCDF file
%         ncid=netcdf.open(filename_L1B,'NC_NOWRITE');
%         
%         % Get dimensions
% %         dimid=netcdf.inqDimID(ncid,'data_20/ku/time');
% %         [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
%         num_surfaces = length(ncread(filename_L1B,'data_20/ku/time_tai'));
% 
%         dimid=netcdf.inqDimID(ncid,'ns');
%         [~,N_samples]=netcdf.inqDim(ncid,dimid);
%         
% 		switch cnf_p.mode
% 			case {'RAW','RMC','HR'}
% 				        dimid=netcdf.inqDimID(ncid,'nl');
% 			case {'LR-RMC'}
% 			        dimid=netcdf.inqDimID(ncid,'nb');
% 		end
% 
%         [~,N_looks_conf_stack]=netcdf.inqDim(ncid,dimid); %total number of looks configured for a stack
%         netcdf.close(ncid);
%         
%         switch cnf_p.mode
%             case {'RAW'}
%                 %RMC mode
%                 data.N_samples=N_samples;
%             case {'RMC','LR-RMC'}
%                 %RMC mode
%                 %data.N_samples=N_samples/2; %(only half of them selected to fit waveform)
%                 data.N_samples=N_samples;
%         end
%         data.N_records=num_surfaces;
%         
%         % -----------------------------------------------------------------
%         % GEO: geographical information
%         % -----------------------------------------------------------------
% 
%         data.GEO.TAI.total             =   ncread(filename_L1B,'data_20/ku/time_tai').';
%         data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%         data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%         data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%         data.GEO.LAT                   =   ncread(filename_L1B,'data_20/ku/latitude').';
%         data.GEO.LON                   =   wrapTo180(ncread(filename_L1B,'data_20/ku/longitude').'); %defined between [0,360]        
%         data.GEO.H_rate                =   ncread(filename_L1B,'data_20/ku/com_altitude_rate').'; % m/s
%         velocity=double(ncread(filename_L1B,'data_20/ku/com_velocity_vector')).'; %vector of vx,vx,vz
%         data.GEO.V                     =   sqrt(velocity(:,1).^2+velocity(:,2).^2+velocity(:,3).^2).';
%         clear velocity;
%         data.GEO.H                     =   ncread(filename_L1B,'data_20/ku/com_altitude').';
%         data.GEO.pitch                 =   (ncread(filename_L1B,'data_20/ku/off_nadir_pitch_angle_pf').').*pi/180;
%         data.GEO.roll                  =   (ncread(filename_L1B,'data_20/ku/off_nadir_roll_angle_pf').').*pi/180;
%         data.GEO.yaw                   =   (ncread(filename_L1B,'data_20/ku/off_nadir_yaw_angle_pf').').*pi/180;
%         
%         
%         % -----------------------------------------------------------------
%         % MEA: measurements
%         % -----------------------------------------------------------------
%         data.MEA.win_delay = double(ncread(filename_L1B,'data_20/ku/tracker_range_calibrated')).'*2.0/cst_p.c_cst;
%         
%         % ---------------------------------------------------------------------
%         % COR: Geophysical corrections
%         % ---------------------------------------------------------------------
% %         % Currently is not foreseen to be included in the L1B-product of
% %         % Sentinel-6
% %         data.COR.surf_type_flag=ncread(filename_L1B,'flag_surface_classification').'; %{0: open ocean, 1: land, 2: continental_water, 3: acquatic vegetation, 4: continental_ice_snow, 5: floating_ice, 6: salted_basin}
% %         
% 
%         %---------------------------------------------------------------
%         % SWH/SIGMA0 instrumental correction
%         %---------------------------------------------------------------
%         data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%         data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%         
%         % -----------------------------------------------------------------
%         % HRM: High-resolution mode: Waveforms
%         % -----------------------------------------------------------------
%         % ------------------ Waveforms ------------------------------------
%         %if ~isempty(strfind(name,'HR'))            
%             switch cnf_p.mode 
%                 case {'SAR','RAW'}
%                     i2q2_meas=double(ncread(filename_L1B,'data_20/ku/hr_power_waveform')).';
%                     scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
%                     %apply scaling factor to waveforms
%                     data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     clear i2q2_meas scale_factor;
%                 case {'RMC'}
%                     %consider only the half of the waveform for fitting
%                     i2q2_meas=double(ncread(filename_L1B,'data_20/ku/hr_power_waveform')).';
%                     scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
%                     %apply scaling factor to waveforms
%                     data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     clear i2q2_meas scale_factor;
%                     %data.HRM.power_wav=data.HRM.power_wav(1:data.N_samples,:); 
%                 case {'LR-RMC'}
%                     i2q2_meas=double(ncread(filename_L1B,'data_20/ku/lr_rmc_power_waveform')).';
%                     scale_factor=ncread(filename_L1B,'data_20/ku/waveform_scale_factor');
%                     %apply scaling factor to waveforms
%                     data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     clear i2q2_meas scale_factor;					
%             end
%             
%             
%             switch cnf_p.mode 
%                 case {'SAR','RAW','RMC','HR'}
% 					data.HRM.Neff       =   ncread(filename_L1B,'data_20/ku/num_looks_start_stop').'; %effective number of beams that form the stack including possuible looks that are set entirely to zero
% 				case {'LR-RMC'}
% 					data.HRM.Neff       =   64*ones(1,data.N_records); %a single waveform look will be considered
%             end
%             data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%             data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%             %             data.HRM.ThN        =   zeros(1,num_surfaces);
%             %             data.HRM.wfm_count  =   1:1:num_surfaces;
%             
%             
%             % ----sigma0 scaling factor for conversion from Power to sigma0
%             %units
%             data.HRM.s0_sf=ncread(filename_L1B,'data_20/ku/sig0_scaling_factor').';
%             
%             
%             %--------------------------------------------------------------
%             %------------- Stack characterization parameters --------------
%             %--------------------------------------------------------------
%             %----- Dopppler mask ------------------------------------------
% 			switch cnf_p.mode 
%                 case {'SAR','RAW','RMC','HR'}
% 					data.HRM.Doppler_mask   =   (ncread(filename_L1B,'data_20/ku/stack_mask_start_stop')+1); %due to the way the stack mask vector is saved in the netcdf the mask is not exactly the same as original
% 				case {'LR-RMC'}
% 					data.HRM.Doppler_mask   =   (ncread(filename_L1B,'data_20/ku/burst_mask')+1); 
% 			end
%             % internally when using the ncread function the scaling factor
%             % of ZP is considered
%             data.HRM.pri_surf=ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'; %in seconds
%             
%             %clock from L1B
%             data.HRM.fs_clock_ku_surf=ncread(filename_L1B,'data_20/ku/altimeter_clock').'; %in Hz
%             
%             %------------- Geometry-related parameters ------------------------
%             if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
% 				switch cnf_p.mode
% 					case {'HR','SAR','RAW','RMC'}
% 						switch cnf_p.looks_index_method
% 							case {'Doppler'}
% 								%This option Does not provide the correct information but it is kept for comparison purposes
% 								switch cnf_p.fd_method
% 									case 'exact'
% 										% exact method would require precise information for
% 										% each beam within the stack (not currently developed for L1B products for SEOMs)
% 										% Dopp_ang= 90-Beam_ang
% 										data.HRM.beam_ang_stack   =  pi/2-ncread(filename_L1BS,'doppler_angle'); %TBD in radians
% 										% it is not foreseen to provide the specific
% 										% PRI/velocity for each look coming from a different
% 										% burst: use the same pri for that surface
% 										data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*...
% 											(ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'); %TBD %eventually for sentinel-6 the PRF could change
% 										data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V;
% 									case 'approximate'
% 										data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'data_20/ku/pointing_angle_start').'); % in radians
% 										data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'data_20/ku/pointing_angle_stop').'); % in radians
% 										%corresponds to the complementary of the
% 										%Beam angle
% 										data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'data_20/ku/doppler_angle_start').'); % in radians
% 										data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'data_20/ku/doppler_angle_stop').'); % in radians
% 								end
% 							case {'Look_angle'}
% 								switch cnf_p.look_ang_method
% 									case 'exact'
% 										% exact method would require precise information for
% 										% each beam within the stack (not currently developed for L1B products for SEaOMs)
% 										data.HRM.look_ang_stack   =  ncread(filename_L1B,'data_20/ku/look_angle'); % in radians
% 										% it is not foreseen to provide the specific
% 										% PRI/velocity for each look coming from a different
% 										% burst: use the same pri for that surface
% 										data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*(ncread(filename_L1B,'data_20/ku/pulse_repetition_interval').'); %TBD %eventually for sentinel-6 the PRF could change
% 										data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V; % TBD norm of velocity for each beam within stack
% 									case 'approximate'
% 										data.HRM.look_ang_start_surf=(ncread(filename_L1B,'data_20/ku/look_angle_start').'); % in radians
% 										data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'data_20/ku/look_angle_stop').'); % in radians
% 								end
% 								
% 						end
% 				
% 				end
%                 
%             end
%             
%             
%             %------------------------------------------------------------------
%             %-------------------GLOBAL ATTRIBUTES -----------------------------
%             %------------------------------------------------------------------            
%             data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=ncreadatt(filename_L1B,'/','mission_name');
%             data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%             data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%             data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';            
%             data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=ncreadatt(filename_L1B,'/','first_measurement_time');
%             data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=ncreadatt(filename_L1B,'/','last_measurement_time');
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=ncreadatt(filename_L1B,'/','xref_altimeter_level0');
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=ncreadatt(filename_L1B,'/','xref_orbit');
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=ncreadatt(filename_L1B,'/','xref_doris_uso');
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available';
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available';
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available';
%             data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=ncreadatt(filename_L1B,'/','xref_characterization');
%             data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(ncreadatt(filename_L1B,'/','semi_major_ellipsoid_axis'));
%             data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(ncreadatt(filename_L1B,'/','ellipsoid_flattening'));
% 
%     case '.mat'
%         %% ------------------------ Matlab data -------------------------------
%         load(filename_L1B);
%         switch cnf_p.mode
%             case {'RAW','RMC','HR'}
%                 
%                 s=size(wfm_cor_i2q2_sar_ku);
%                 N_samples=s(2);
%                 num_surfaces=s(1);
%                 data.N_samples=N_samples;
%                 data.N_records=num_surfaces;
%                 
%                 % -----------------------------------------------------------------
%                 % GEO: geographical information
%                 % -----------------------------------------------------------------
%                 data.GEO.TAI.total             =   time_surf;
%                 data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%                 data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%                 data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%                 data.GEO.LAT                   =   lat_surf; %degrees
%                 data.GEO.LON                   =   wrapTo180(lon_surf); %degrees defined between [-180,180]
%                 data.GEO.H_rate                =   alt_rate_sat; % m/s
%                 data.GEO.V                     =   sqrt(x_vel_sat.^2+y_vel_sat.^2+z_vel_sat.^2); %m/s
%                 data.GEO.H                     =   alt_sat; % m
%                 data.GEO.pitch                 =   pitch_surf; % rad
%                 data.GEO.roll                  =   roll_surf; % rad
%                 data.GEO.yaw                   =   yaw_surf; % rad
%                 
%                 % -----------------------------------------------------------------
%                 % MEA: measurements
%                 % -----------------------------------------------------------------
%                 data.MEA.win_delay = win_delay_surf;
%                 
%                 % ---------------------------------------------------------------------
%                 % COR: Geophysical corrections
%                 % ---------------------------------------------------------------------
%                 % Currently is not foreseen to be included in the L1B-product of
%                 % Sentinel-6
%                 
%                 
%                 %---------------------------------------------------------------
%                 % SWH/SIGMA0 instrumental correction
%                 %---------------------------------------------------------------
%                 data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%                 data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%                 
%                 
%                 % -----------------------------------------------------------------
%                 % HRM: High-resolution mode: Waveforms
%                 % -----------------------------------------------------------------
%                 % ------------------ Waveforms ------------------------------------
%                 if ~isempty(strfind(name,'HR'))
%                     % ------------ SAR mode ---------------------------------------
%                     data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
%                     
%                     
%                     data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
%                     data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%                     data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%                     %             data.HRM.ThN        =   zeros(1,num_surfaces);
%                     %             data.HRM.wfm_count  =   1:1:num_surfaces;
%                     
%                     
%                     % ----sigma0 scaling factor for conversion from Power to sigma0
%                     %units
%                     data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
%                     
%                     
%                     %--------------------------------------------------------------
%                     %------------- Stack characterization parameters --------------
%                     %--------------------------------------------------------------
%                     %----- Dopppler mask ------------------------------------------
%                     data.HRM.Doppler_mask   =   stack_mask_vector_start_stop.';
%                     data.HRM.pri_surf=pri_sar_surf; %in seconds
%                     data.HRM.fs_clock_ku_surf=1./T0_sar_surf; %in Hz
%                     
%                     %------------- Geometry-related parameters ------------------------
%                     if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
%                         switch cnf_p.looks_index_method
%                             case {'Doppler'}
%                                 %This option Does not provide the correct information but it is kept for comparison purposes
%                                 switch cnf_p.fd_method
%                                     case 'exact'
%                                         % exact method would require precise information for
%                                         % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                         data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
%                                         data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                         data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
%                                     case 'approximate'
%                                         data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
%                                         data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
%                                         data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
%                                         data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
%                                 end
%                             case {'Look_angle'}
%                                 switch cnf_p.look_ang_method
%                                     case 'exact'
%                                         % exact method would require precise information for
%                                         % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                         data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
%                                         data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                         data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
%                                     case 'approximate'
%                                         data.HRM.look_ang_start_surf=start_look_angle; % in radians
%                                         data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
%                                 end
%                                 
%                         end
%                     end
%                     
%                     %------------------------------------------------------------------
%                     %-------------------GLOBAL ATTRIBUTES -----------------------------
%                     %------------------------------------------------------------------
% 
%                     data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
%                                         
%                     TAI_2012 = (12*365+4+181)*3600*24 + 35;
%                     TAI_2015 = (15*365+4+181)*3600*24 + 36;
%                     time_UTC_init = double(min(data.GEO.TAI.total));
%                     time_UTC_end  = double(max(data.GEO.TAI.total));
%                     
%                     % add leap seconds to the TAI time. Only valid for
%                     if(time_UTC_init < TAI_2012)
%                         time_UTC_init = time_UTC_init - 34;                        
%                     elseif(time_UTC_init > TAI_2015)
%                         time_UTC_init = time_UTC_init - 36;
%                     else
%                         time_UTC_init = time_UTC_init - 35;
%                     end
%                     if(time_UTC_end < TAI_2012)
%                         time_UTC_end = time_UTC_end - 34;
%                     elseif(time_UTC_end > TAI_2015)
%                         time_UTC_end = time_UTC_end - 36;
%                     else
%                         time_UTC_end = time_UTC_end - 35;
%                     end
%                     
%                     data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
%                     data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
%                     data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
%                     data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
%                                         
%                     
%                 elseif ~isempty(strfind(name,'LR'))
%                     % ----------- LR: low resolution data extraction --------------
%                 end
%                 
%             case {'LR-RMC'}
% 
%                 s=size(wfm_i2q2_sar_ku_combined);
%                 N_samples=s(2);
%                 num_surfaces=s(1);
%                 data.N_samples=N_samples;
%                 data.N_records=num_surfaces;
%                 
%                 % -----------------------------------------------------------------
%                 % GEO: geographical information
%                 % -----------------------------------------------------------------
%                 data.GEO.TAI.total             =   time_aux_seconds;
%                 data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%                 data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%                 data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%                 data.GEO.LAT                   =   latitude_ku; %degrees
%                 data.GEO.LON                   =   wrapTo180(longitude_ku); %degrees defined between [-180,180]
%                 data.GEO.H_rate                =   com_altitude_rate_ku; % m/s
%                 data.GEO.V                     =   sqrt(vel_x.^2+vel_y.^2+vel_z.^2); %m/s
%                 data.GEO.H                     =   com_altitude_ku; % m
%                 data.GEO.pitch                 =   off_nadir_pitch_angle_pf_20_ku; % rad
%                 data.GEO.roll                  =   off_nadir_roll_angle_pf_20_ku; % rad
%                 data.GEO.yaw                   =   off_nadir_yaw_angle_pf_20_ku; % rad
%                 
%                 % -----------------------------------------------------------------
%                 % MEA: measurements
%                 % -----------------------------------------------------------------
%                 data.MEA.win_delay = altimeter_range_calibrated_ku*2/cst_p.c_cst;
%                 
%                 % ---------------------------------------------------------------------
%                 % COR: Geophysical corrections
%                 % ---------------------------------------------------------------------
%                 % Currently is not foreseen to be included in the L1B-product of
%                 % Sentinel-6
%                 
%                 
%                 %---------------------------------------------------------------
%                 % SWH/SIGMA0 instrumental correction
%                 %---------------------------------------------------------------
%                 data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%                 data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%                 
%                 
%                 % -----------------------------------------------------------------
%                 % HRM: High-resolution mode: Waveforms
%                 % -----------------------------------------------------------------
%                 % ------------------ Waveforms ------------------------------------
%                 
%                 % ------------ SAR mode ---------------------------------------
%                 data.HRM.power_wav=wfm_i2q2_sar_ku.';
%                 
%                 
%                 data.HRM.Neff       =   N_beams_start_stop_combined; %effective number of beams that form the stack including possuible looks that are set entirely to zero
%                 data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%                 data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%                 %             data.HRM.ThN        =   zeros(1,num_surfaces);
%                 %             data.HRM.wfm_count  =   1:1:num_surfaces;
%                 
%                 
%                 % ----sigma0 scaling factor for conversion from Power to sigma0
%                 %units
%                 data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
%                 
%                 
%                 %--------------------------------------------------------------
%                 %------------- Stack characterization parameters --------------
%                 %--------------------------------------------------------------
%                 %----- Dopppler mask ------------------------------------------
%                 data.HRM.Doppler_mask   =   stack_mask_vector_start_stop_combined.';
%                 data.HRM.pri_surf=double(pri_lrm_l1b_ku); %in seconds
%                 data.HRM.fs_clock_ku_surf=1./double(data_20/ku/altimeter_clock_ku); %in Hz
%                 
%                 %------------- Geometry-related parameters ------------------------
%                 if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
%                     switch cnf_p.looks_index_method
%                         case {'Doppler'}
%                             
%                         case {'Look_angle'}
%                             switch cnf_p.look_ang_method
%                                 case 'exact'
%                                     % exact method would require precise information for
%                                     % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                     data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
%                                     data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                     data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
%                                 case 'approximate'
%                                     data.HRM.look_ang_start_surf=start_look_angle_combined; % in radians
%                                     data.HRM.look_ang_stop_surf=stop_look_angle_combined; % in radians
%                             end
%                             
%                     end
%                 end
%                 
%                 %------------------------------------------------------------------
%                 %-------------------GLOBAL ATTRIBUTES -----------------------------
%                 %------------------------------------------------------------------
%                 
%                 data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
%                 
%                 TAI_2012 = (12*365+4+181)*3600*24 + 35;
%                 TAI_2015 = (15*365+4+181)*3600*24 + 36;
%                 time_UTC_init = double(min(data.GEO.TAI.total));
%                 time_UTC_end  = double(max(data.GEO.TAI.total));
%                 
%                 % add leap seconds to the TAI time. Only valid for
%                 if(time_UTC_init < TAI_2012)
%                     time_UTC_init = time_UTC_init - 34;
%                 elseif(time_UTC_init > TAI_2015)
%                     time_UTC_init = time_UTC_init - 36;
%                 else
%                     time_UTC_init = time_UTC_init - 35;
%                 end
%                 if(time_UTC_end < TAI_2012)
%                     time_UTC_end = time_UTC_end - 34;
%                 elseif(time_UTC_end > TAI_2015)
%                     time_UTC_end = time_UTC_end - 36;
%                 else
%                     time_UTC_end = time_UTC_end - 35;
%                 end
%                 
%                 data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
%                 data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
%                 data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
%                 data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
%                 
%                 
%                 
%         end              
%         
%     otherwise
%         error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
% end
% 
% end
% 
% 
% % function data = readL1B_S6_ISD (filename_L1B,cnf_p,cst_p,varargin)
% % % -------------------------------------------------------------------------
% % % Created by isardSAT S.L. 
% % % -------------------------------------------------------------------------
% % % This code allows for reading altimetry data from L1B Sentinel-6 products
% % % processed by ISD
% % %
% % % -------------------------------------------------------------------------
% % % 
% % % Author:           Eduard Makhoul / isardSAT
% % %
% % % Reviewer:         M�nica Roca / isardSAT
% % %
% % % Last revision:    Eduard Makhoul / isardSAT V1 15/06/2016
% % % This software is built with internal funding 
% % % -------------------------------------------------------------------------
% % % -------------------------------------------------------------------------
% % % INPUT:
% % %       MANDATORY:
% % %           -filename_L1B    =   L1B filename with the fullpath information
% % %           -cnf_p           =   Structure with configuration parameters of the
% % %                                L2 processor
% % %       OPTIONAL:
% % %           -filename_L1BS   =   L1B-S filename with the fullpath name
% % %       
% % % 
% % % OUTPUT:
% % %       data        =   structure of data as defined by our L2 processor
% % %
% % % -------------------------------------------------------------------------
% % % -------------------------------------------------------------------------
% % % CALLED FUNCTIONS/ROUTINES
% % %
% % %
% % % -------------------------------------------------------------------------
% % % -------------------------------------------------------------------------
% % % COMMENTS/RESTRICTIONS:
% % %
% % % -------------------------------------------------------------------------  
% % % -------------------------------------------------------------------------
% % % Versions control:
% % % v1.0: Based on read_alt_data.m
% % % v1.1: Include global attributes sensor name, acq. station, acq date,... for
% % % the reference L1B product & orbital parameter info in a similar fashion
% % % as the global attribtues provided for output L1B product according to
% % % Sentinel-3 format and within the SEOMs projects. It is considered only
% % % for the final netCDF L1B products as defined in the PSF issue 1.3 
% % %(not all fields considered in SEOMS Sentinel-3 are available in the Sentinel-6 product format netcdf)
% 
% % %% ---------------- Handling input variables ------------------------------
% % if(nargin<2 || nargin>(3+1*2))
%     % error('Wrong number of input parameters');   
% % end
% % %option to include the L1B_S product to read the exact info from Look to
% % %Look within the stack
% % p = inputParser;
% % p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
% % p.parse(varargin{:});
% % filename_L1BS=char(p.Results.filename_L1BS);
% % clear p;
% 
% % %% ------------------------------------------------------------------------- 
% % % Loading data L1B
% % % ------------------------------------------------------------------------- 
% % [~,name_file,ext]=fileparts(filename_L1B);
% % ext=lower(ext);
% % name='HR';
% % switch ext
%     % case '.nc'
%         % %% --------------------- NetCDF ---------------------------------------
% 	
%         % % Based on the L1B product generated for DeDop a la Sentinel-3
%         % % open the netCDF file
%         % ncid=netcdf.open(filename_L1B,'NC_NOWRITE');
%         
%         % % Get dimensions
%         % dimid=netcdf.inqDimID(ncid,'time_20_ku');
%         % [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
%         % dimid=netcdf.inqDimID(ncid,'ns');
%         % [~,N_samples]=netcdf.inqDim(ncid,dimid);
% 		% switch cnf_p.mode
% 			% case {'RAW','RMC','HR'}
% 				        % dimid=netcdf.inqDimID(ncid,'nl');
% 			% case {'LR-RMC'}
% 			        % dimid=netcdf.inqDimID(ncid,'nb');
% 		% end
% 
%         % [~,N_looks_conf_stack]=netcdf.inqDim(ncid,dimid); %total number of looks configured for a stack
%         % netcdf.close(ncid);
%         
%         % switch cnf_p.mode
%             % case {'RAW'}
%                 % %RMC mode
%                 % data.N_samples=N_samples;
%             % case {'RMC','LR-RMC'}
%                 % %RMC mode
%                 % %data.N_samples=N_samples/2; %(only half of them selected to fit waveform)
%                 % data.N_samples=N_samples;
%         % end
%         % data.N_records=num_surfaces;
%         
%         % % -----------------------------------------------------------------
%         % % GEO: geographical information
%         % % -----------------------------------------------------------------
% 
%         % data.GEO.TAI.total             =   ncread(filename_L1B,'time_tai_20_ku').';
%         % data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%         % data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%         % data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%         % data.GEO.LAT                   =   ncread(filename_L1B,'latitude_20_ku').';
%         % data.GEO.LON                   =   wrapTo180(ncread(filename_L1B,'longitude_20_ku').'); %defined between [0,360]        
%         % data.GEO.H_rate                =   ncread(filename_L1B,'com_altitude_rate_20_ku').'; % m/s
%         % velocity=double(ncread(filename_L1B,'com_velocity_vector_20_ku')).'; %vector of vx,vx,vz
%         % data.GEO.V                     =   sqrt(velocity(:,1).^2+velocity(:,2).^2+velocity(:,3).^2).';
%         % clear velocity;
%         % data.GEO.H                     =   ncread(filename_L1B,'com_altitude_20_ku').';
%         % data.GEO.pitch                 =   (ncread(filename_L1B,'off_nadir_pitch_angle_pf_20_ku').').*pi/180;
%         % data.GEO.roll                  =   (ncread(filename_L1B,'off_nadir_roll_angle_pf_20_ku').').*pi/180;
%         % data.GEO.yaw                   =   (ncread(filename_L1B,'off_nadir_yaw_angle_pf_20_ku').').*pi/180;
%         
%         
%         % % -----------------------------------------------------------------
%         % % MEA: measurements
%         % % -----------------------------------------------------------------
%         % data.MEA.win_delay = double(ncread(filename_L1B,'tracker_range_calibrated_20_ku')).'*2.0/cst_p.c_cst;
%         
%         % % ---------------------------------------------------------------------
%         % % COR: Geophysical corrections
%         % % ---------------------------------------------------------------------
% % %         % Currently is not foreseen to be included in the L1B-product of
% % %         % Sentinel-6
% % %         data.COR.surf_type_flag=ncread(filename_L1B,'flag_surface_classification_20_ku ').'; %{0: open ocean, 1: land, 2: continental_water, 3: acquatic vegetation, 4: continental_ice_snow, 5: floating_ice, 6: salted_basin}
% % %         
% 
%         % %---------------------------------------------------------------
%         % % SWH/SIGMA0 instrumental correction
%         % %---------------------------------------------------------------
%         % data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%         % data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%         
%         % % -----------------------------------------------------------------
%         % % HRM: High-resolution mode: Waveforms
%         % % -----------------------------------------------------------------
%         % % ------------------ Waveforms ------------------------------------
%         % %if ~isempty(strfind(name,'HR'))            
%             % switch cnf_p.mode 
%                 % case {'SAR','RAW'}
%                     % i2q2_meas=double(ncread(filename_L1B,'hr_power_waveform_20_ku')).';
%                     % scale_factor=ncread(filename_L1B,'waveform_scale_factor_20_ku');
%                     % %apply scaling factor to waveforms
%                     % data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     % clear i2q2_meas scale_factor;
%                 % case {'RMC'}
%                     % %consider only the half of the waveform for fitting
%                     % i2q2_meas=double(ncread(filename_L1B,'hr_power_waveform_20_ku')).';
%                     % scale_factor=ncread(filename_L1B,'waveform_scale_factor_20_ku');
%                     % %apply scaling factor to waveforms
%                     % data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     % clear i2q2_meas scale_factor;
%                     % %data.HRM.power_wav=data.HRM.power_wav(1:data.N_samples,:); 
%                 % case {'LR-RMC'}
%                     % i2q2_meas=double(ncread(filename_L1B,'lr_rmc_power_waveform_20_ku')).';
%                     % scale_factor=ncread(filename_L1B,'waveform_scale_factor_20_ku');
%                     % %apply scaling factor to waveforms
%                     % data.HRM.power_wav=(i2q2_meas.*repmat(scale_factor,1,N_samples)).';
%                     % clear i2q2_meas scale_factor;					
%             % end
%             
%             
%             % switch cnf_p.mode 
%                 % case {'SAR','RAW','RMC','HR'}
% 					% data.HRM.Neff       =   ncread(filename_L1B,'num_looks_start_stop_20_ku').'; %effective number of beams that form the stack including possuible looks that are set entirely to zero
% 				% case {'LR-RMC'}
% 					% data.HRM.Neff       =   64*ones(1,data.N_records); %a single waveform look will be considered
%             % end
%             % data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%             % data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%             % %             data.HRM.ThN        =   zeros(1,num_surfaces);
%             % %             data.HRM.wfm_count  =   1:1:num_surfaces;
%             
%             
%             % % ----sigma0 scaling factor for conversion from Power to sigma0
%             % %units
%             % data.HRM.s0_sf=ncread(filename_L1B,'sigma0_scaling_factor_20_ku').';
%             
%             
%             % %--------------------------------------------------------------
%             % %------------- Stack characterization parameters --------------
%             % %--------------------------------------------------------------
%             % %----- Dopppler mask ------------------------------------------
% 			% switch cnf_p.mode 
%                 % case {'SAR','RAW','RMC','HR'}
% 					% data.HRM.Doppler_mask   =   (ncread(filename_L1B,'stack_mask_start_stop_20_ku')+1); %due to the way the stack mask vector is saved in the netcdf the mask is not exactly the same as original
% 				% case {'LR-RMC'}
% 					% data.HRM.Doppler_mask   =   (ncread(filename_L1B,'burst_mask_20_ku')+1); 
% 			% end
%             % % internally when using the ncread function the scaling factor
%             % % of ZP is considered
%             % data.HRM.pri_surf=ncread(filename_L1B,'pulse_repetition_interval_20_ku').'; %in seconds
%             
%             % %clock from L1B
%             % data.HRM.fs_clock_ku_surf=ncread(filename_L1B,'altimeter_clock_20_ku').'; %in Hz
%             
%             % %------------- Geometry-related parameters ------------------------
%             % if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
% 				% switch cnf_p.mode
% 					% case {'HR','SAR','RAW','RMC'}
% 						% switch cnf_p.looks_index_method
% 							% case {'Doppler'}
% 								% %This option Does not provide the correct information but it is kept for comparison purposes
% 								% switch cnf_p.fd_method
% 									% case 'exact'
% 										% % exact method would require precise information for
% 										% % each beam within the stack (not currently developed for L1B products for SEOMs)
% 										% % Dopp_ang= 90-Beam_ang
% 										% data.HRM.beam_ang_stack   =  pi/2-ncread(filename_L1BS,'doppler_angle_20_ku'); %TBD in radians
% 										% % it is not foreseen to provide the specific
% 										% % PRI/velocity for each look coming from a different
% 										% % burst: use the same pri for that surface
% 										% data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*...
% 											% (ncread(filename_L1B,'pulse_repetition_interval_20_ku').'); %TBD %eventually for sentinel-6 the PRF could change
% 										% data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V;
% 									% case 'approximate'
% 										% data.HRM.pointing_ang_start_surf=(ncread(filename_L1B,'pointing_angle_start_20_ku').'); % in radians
% 										% data.HRM.pointing_ang_stop_surf=(ncread(filename_L1B,'pointing_angle_stop_20_ku').'); % in radians
% 										% %corresponds to the complementary of the
% 										% %Beam angle
% 										% data.HRM.doppler_angle_start_surf=(ncread(filename_L1B,'doppler_angle_start_20_ku').'); % in radians
% 										% data.HRM.doppler_angle_stop_surf=(ncread(filename_L1B,'doppler_angle_stop_20_ku').'); % in radians
% 								% end
% 							% case {'Look_angle'}
% 								% switch cnf_p.look_ang_method
% 									% case 'exact'
% 										% % exact method would require precise information for
% 										% % each beam within the stack (not currently developed for L1B products for SEOMs)
% 										% data.HRM.look_ang_stack   =  ncread(filename_L1B,'look_angle_20_ku'); % in radians
% 										% % it is not foreseen to provide the specific
% 										% % PRI/velocity for each look coming from a different
% 										% % burst: use the same pri for that surface
% 										% data.HRM.pri_stack        =  ones(N_looks_conf_stack,1)*(ncread(filename_L1B,'pulse_repetition_interval_20_ku').'); %TBD %eventually for sentinel-6 the PRF could change
% 										% data.GEO.V_stack          =  ones(N_looks_conf_stack,1)*data.GEO.V; % TBD norm of velocity for each beam within stack
% 									% case 'approximate'
% 										% data.HRM.look_ang_start_surf=(ncread(filename_L1B,'look_angle_start_20_ku').'); % in radians
% 										% data.HRM.look_ang_stop_surf=(ncread(filename_L1B,'look_angle_stop_20_ku').'); % in radians
% 								% end
% 								
% 						% end
% 				
% 				% end
%                 
%             % end
%             
%             
%             % %------------------------------------------------------------------
%             % %-------------------GLOBAL ATTRIBUTES -----------------------------
%             % %------------------------------------------------------------------            
%             % data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name=ncreadatt(filename_L1B,'/','altimeter_name');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%             % data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%             % data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';            
%             % data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=ncreadatt(filename_L1B,'/','first_measurement_time');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=ncreadatt(filename_L1B,'/','last_measurement_time');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=ncreadatt(filename_L1B,'/','input_file_p4_isp');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit=ncreadatt(filename_L1B,'/','input_file_orbit');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO=ncreadatt(filename_L1B,'/','input_file_uso_drift');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available';
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available';
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available';
%             % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation=ncreadatt(filename_L1B,'/','input_file_characterisation');
%             % data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(ncreadatt(filename_L1B,'/','semi_major_ellipsoid_axis'));
%             % data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(ncreadatt(filename_L1B,'/','ellipsoid_flattening'));
% 
%     % case '.mat'
%         % %% ------------------------ Matlab data -------------------------------
%         % load(filename_L1B);
%         % switch cnf_p.mode
%             % case {'RAW','RMC','HR'}
%                 
%                 % s=size(wfm_cor_i2q2_sar_ku);
%                 % N_samples=s(2);
%                 % num_surfaces=s(1);
%                 % data.N_samples=N_samples;
%                 % data.N_records=num_surfaces;
%                 
%                 % % -----------------------------------------------------------------
%                 % % GEO: geographical information
%                 % % -----------------------------------------------------------------
%                 % data.GEO.TAI.total             =   time_surf;
%                 % data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%                 % data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%                 % data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%                 % data.GEO.LAT                   =   lat_surf; %degrees
%                 % data.GEO.LON                   =   wrapTo180(lon_surf); %degrees defined between [-180,180]
%                 % data.GEO.H_rate                =   alt_rate_sat; % m/s
%                 % data.GEO.V                     =   sqrt(x_vel_sat.^2+y_vel_sat.^2+z_vel_sat.^2); %m/s
%                 % data.GEO.H                     =   alt_sat; % m
%                 % data.GEO.pitch                 =   pitch_surf; % rad
%                 % data.GEO.roll                  =   roll_surf; % rad
%                 % data.GEO.yaw                   =   yaw_surf; % rad
%                 
%                 % % -----------------------------------------------------------------
%                 % % MEA: measurements
%                 % % -----------------------------------------------------------------
%                 % data.MEA.win_delay = win_delay_surf;
%                 
%                 % % ---------------------------------------------------------------------
%                 % % COR: Geophysical corrections
%                 % % ---------------------------------------------------------------------
%                 % % Currently is not foreseen to be included in the L1B-product of
%                 % % Sentinel-6
%                 
%                 
%                 % %---------------------------------------------------------------
%                 % % SWH/SIGMA0 instrumental correction
%                 % %---------------------------------------------------------------
%                 % data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%                 % data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%                 
%                 
%                 % % -----------------------------------------------------------------
%                 % % HRM: High-resolution mode: Waveforms
%                 % % -----------------------------------------------------------------
%                 % % ------------------ Waveforms ------------------------------------
%                 % if ~isempty(strfind(name,'HR'))
%                     % % ------------ SAR mode ---------------------------------------
%                     % data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
%                     
%                     
%                     % data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
%                     % data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%                     % data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%                     % %             data.HRM.ThN        =   zeros(1,num_surfaces);
%                     % %             data.HRM.wfm_count  =   1:1:num_surfaces;
%                     
%                     
%                     % % ----sigma0 scaling factor for conversion from Power to sigma0
%                     % %units
%                     % data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
%                     
%                     
%                     % %--------------------------------------------------------------
%                     % %------------- Stack characterization parameters --------------
%                     % %--------------------------------------------------------------
%                     % %----- Dopppler mask ------------------------------------------
%                     % data.HRM.Doppler_mask   =   stack_mask_vector_start_stop.';
%                     % data.HRM.pri_surf=pri_sar_surf; %in seconds
%                     % data.HRM.fs_clock_ku_surf=1./T0_sar_surf; %in Hz
%                     
%                     % %------------- Geometry-related parameters ------------------------
%                     % if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
%                         % switch cnf_p.looks_index_method
%                             % case {'Doppler'}
%                                 % %This option Does not provide the correct information but it is kept for comparison purposes
%                                 % switch cnf_p.fd_method
%                                     % case 'exact'
%                                         % % exact method would require precise information for
%                                         % % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                         % data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
%                                         % data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                         % data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
%                                     % case 'approximate'
%                                         % data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
%                                         % data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
%                                         % data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
%                                         % data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
%                                 % end
%                             % case {'Look_angle'}
%                                 % switch cnf_p.look_ang_method
%                                     % case 'exact'
%                                         % % exact method would require precise information for
%                                         % % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                         % data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
%                                         % data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                         % data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
%                                     % case 'approximate'
%                                         % data.HRM.look_ang_start_surf=start_look_angle; % in radians
%                                         % data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
%                                 % end
%                                 
%                         % end
%                     % end
%                     
%                     % %------------------------------------------------------------------
%                     % %-------------------GLOBAL ATTRIBUTES -----------------------------
%                     % %------------------------------------------------------------------
% 
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
%                                         
%                     % TAI_2012 = (12*365+4+181)*3600*24 + 35;
%                     % TAI_2015 = (15*365+4+181)*3600*24 + 36;
%                     % time_UTC_init = double(min(data.GEO.TAI.total));
%                     % time_UTC_end  = double(max(data.GEO.TAI.total));
%                     
%                     % % add leap seconds to the TAI time. Only valid for
%                     % if(time_UTC_init < TAI_2012)
%                         % time_UTC_init = time_UTC_init - 34;                        
%                     % elseif(time_UTC_init > TAI_2015)
%                         % time_UTC_init = time_UTC_init - 36;
%                     % else
%                         % time_UTC_init = time_UTC_init - 35;
%                     % end
%                     % if(time_UTC_end < TAI_2012)
%                         % time_UTC_end = time_UTC_end - 34;
%                     % elseif(time_UTC_end > TAI_2015)
%                         % time_UTC_end = time_UTC_end - 36;
%                     % else
%                         % time_UTC_end = time_UTC_end - 35;
%                     % end
%                     
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
%                     % data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
%                                         
%                     
%                 % elseif ~isempty(strfind(name,'LR'))
%                     % % ----------- LR: low resolution data extraction --------------
%                 % end
%                 
%             % case {'LR-RMC'}
% 
%                 % s=size(wfm_i2q2_sar_ku_combined);
%                 % N_samples=s(2);
%                 % num_surfaces=s(1);
%                 % data.N_samples=N_samples;
%                 % data.N_records=num_surfaces;
%                 
%                 % % -----------------------------------------------------------------
%                 % % GEO: geographical information
%                 % % -----------------------------------------------------------------
%                 % data.GEO.TAI.total             =   time_aux_seconds;
%                 % data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
%                 % data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
%                 % data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
%                 % data.GEO.LAT                   =   latitude_ku; %degrees
%                 % data.GEO.LON                   =   wrapTo180(longitude_ku); %degrees defined between [-180,180]
%                 % data.GEO.H_rate                =   com_altitude_rate_ku; % m/s
%                 % data.GEO.V                     =   sqrt(vel_x.^2+vel_y.^2+vel_z.^2); %m/s
%                 % data.GEO.H                     =   com_altitude_ku; % m
%                 % data.GEO.pitch                 =   off_nadir_pitch_angle_pf_20_ku; % rad
%                 % data.GEO.roll                  =   off_nadir_roll_angle_pf_20_ku; % rad
%                 % data.GEO.yaw                   =   off_nadir_yaw_angle_pf_20_ku; % rad
%                 
%                 % % -----------------------------------------------------------------
%                 % % MEA: measurements
%                 % % -----------------------------------------------------------------
%                 % data.MEA.win_delay = altimeter_range_calibrated_ku*2/cst_p.c_cst;
%                 
%                 % % ---------------------------------------------------------------------
%                 % % COR: Geophysical corrections
%                 % % ---------------------------------------------------------------------
%                 % % Currently is not foreseen to be included in the L1B-product of
%                 % % Sentinel-6
%                 
%                 
%                 % %---------------------------------------------------------------
%                 % % SWH/SIGMA0 instrumental correction
%                 % %---------------------------------------------------------------
%                 % data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
%                 % data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
%                 
%                 
%                 % % -----------------------------------------------------------------
%                 % % HRM: High-resolution mode: Waveforms
%                 % % -----------------------------------------------------------------
%                 % % ------------------ Waveforms ------------------------------------
%                 
%                 % % ------------ SAR mode ---------------------------------------
%                 % data.HRM.power_wav=wfm_i2q2_sar_ku.';
%                 
%                 
%                 % data.HRM.Neff       =   N_beams_start_stop_combined; %effective number of beams that form the stack including possuible looks that are set entirely to zero
%                 % data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
%                 % data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
%                 % %             data.HRM.ThN        =   zeros(1,num_surfaces);
%                 % %             data.HRM.wfm_count  =   1:1:num_surfaces;
%                 
%                 
%                 % % ----sigma0 scaling factor for conversion from Power to sigma0
%                 % %units
%                 % data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
%                 
%                 
%                 % %--------------------------------------------------------------
%                 % %------------- Stack characterization parameters --------------
%                 % %--------------------------------------------------------------
%                 % %----- Dopppler mask ------------------------------------------
%                 % data.HRM.Doppler_mask   =   stack_mask_vector_start_stop_combined.';
%                 % data.HRM.pri_surf=double(pri_lrm_l1b_ku); %in seconds
%                 % data.HRM.fs_clock_ku_surf=1./double(altimeter_clock_ku); %in Hz
%                 
%                 % %------------- Geometry-related parameters ------------------------
%                 % if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
%                     % switch cnf_p.looks_index_method
%                         % case {'Doppler'}
%                             
%                         % case {'Look_angle'}
%                             % switch cnf_p.look_ang_method
%                                 % case 'exact'
%                                     % % exact method would require precise information for
%                                     % % each beam within the stack (not currently developed for L1B products for SEOMs)
%                                     % data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
%                                     % data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
%                                     % data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
%                                 % case 'approximate'
%                                     % data.HRM.look_ang_start_surf=start_look_angle_combined; % in radians
%                                     % data.HRM.look_ang_stop_surf=stop_look_angle_combined; % in radians
%                             % end
%                             
%                     % end
%                 % end
%                 
%                 % %------------------------------------------------------------------
%                 % %-------------------GLOBAL ATTRIBUTES -----------------------------
%                 % %------------------------------------------------------------------
%                 
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.altimeter_sensor_name='P4-A';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.gnss_sensor_name='Not available';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.doris_sensor_name='Not available';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.acq_station_name='Not available';
%                 
%                 % TAI_2012 = (12*365+4+181)*3600*24 + 35;
%                 % TAI_2015 = (15*365+4+181)*3600*24 + 36;
%                 % time_UTC_init = double(min(data.GEO.TAI.total));
%                 % time_UTC_end  = double(max(data.GEO.TAI.total));
%                 
%                 % % add leap seconds to the TAI time. Only valid for
%                 % if(time_UTC_init < TAI_2012)
%                     % time_UTC_init = time_UTC_init - 34;
%                 % elseif(time_UTC_init > TAI_2015)
%                     % time_UTC_init = time_UTC_init - 36;
%                 % else
%                     % time_UTC_init = time_UTC_init - 35;
%                 % end
%                 % if(time_UTC_end < TAI_2012)
%                     % time_UTC_end = time_UTC_end - 34;
%                 % elseif(time_UTC_end > TAI_2015)
%                     % time_UTC_end = time_UTC_end - 36;
%                 % else
%                     % time_UTC_end = time_UTC_end - 35;
%                 % end
%                 
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.first_meas_time=datestr(time_UTC_init,'yyyy-mm-dd THH:MM:SS');
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.last_meas_time=datestr(time_UTC_end,'yyyy-mm-dd THH:MM:SS');
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0=strrep(strcat(name_file,ext),'1B','TM');
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=strcat(name_file,ext);
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_orbit='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_doris_USO='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_sar_cal1='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_ku_cal2='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_ltm_c_cal2='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_characterisation='Not available from @Matlab';
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.semi_major_ellipsoid_axis=num2str(cst_p.flat_coeff_cst);
%                 % data.GLOBAL_ATT.DATA_FILE_INFO.ellipsoid_flattening=num2str(cst_p.semi_major_axis_cst);
%                 
%                 
%                 
%         % end              
%         
%     % otherwise
%         % error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
% % end
% 
% % end
% 
% % % LR-RMC for single waveforms at burst level not averaged by RC
% % %                 s=size(wfm_cor_i2q2_sar_ku);
% % %                 N_samples=s(2);
% % %                 num_surfaces=s(1);
% % %                 data.N_samples=N_samples;
% % %                 data.N_records=num_surfaces;
% % %                 
% % %                 % -----------------------------------------------------------------
% % %                 % GEO: geographical information
% % %                 % -----------------------------------------------------------------
% % %                 data.GEO.TAI.total             =   time_sar_ku;
% % %                 data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
% % %                 data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
% % %                 data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
% % %                 data.GEO.LAT                   =   lat_sar_sat.'; %degrees
% % %                 data.GEO.LON                   =   wrapTo180(lon_sar_sat.'); %degrees defined between [-180,180]
% % %                 data.GEO.H_rate                =   alt_rate_sar_sat; % m/s
% % %                 data.GEO.V                     =   sqrt(x_vel_sat_sar.^2+y_vel_sat_sar.^2+z_vel_sat_sar.^2); %m/s
% % %                 data.GEO.H                     =   alt_sar_sat.'; % m
% % %                 data.GEO.pitch                 =   pitch_sar; % rad
% % %                 data.GEO.roll                  =   roll_sar; % rad
% % %                 data.GEO.yaw                   =   yaw_sar; % rad
% % %                 
% % %                 % -----------------------------------------------------------------
% % %                 % MEA: measurements
% % %                 % -----------------------------------------------------------------
% % %                 data.MEA.win_delay = win_delay_sar_ku;
% % %                 
% % %                 % ---------------------------------------------------------------------
% % %                 % COR: Geophysical corrections
% % %                 % ---------------------------------------------------------------------
% % %                 % Currently is not foreseen to be included in the L1B-product of
% % %                 % Sentinel-6
% % %                 
% % %                 
% % %                 %---------------------------------------------------------------
% % %                 % SWH/SIGMA0 instrumental correction
% % %                 %---------------------------------------------------------------
% % %                 data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
% % %                 data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
% % %                 
% % %                 
% % %                 % -----------------------------------------------------------------
% % %                 % HRM: High-resolution mode: Waveforms
% % %                 % -----------------------------------------------------------------
% % %                 % ------------------ Waveforms ------------------------------------
% % %                 
% % %                 % ------------ SAR mode ---------------------------------------
% % %                 data.HRM.power_wav=wfm_cor_i2q2_sar_ku.';
% % %                 
% % %                 
% % %                 data.HRM.Neff       =   N_beams_start_stop; %effective number of beams that form the stack including possuible looks that are set entirely to zero
% % %                 data.HRM.FLAG.mlQ   =   zeros(1,num_surfaces); % if 0 no error if 1 a error ocurred in the stack or multilook
% % %                 data.HRM.FLAG.pQ    =   ones(1,num_surfaces); % if 1 then ok, 0 error in the power
% % %                 %             data.HRM.ThN        =   zeros(1,num_surfaces);
% % %                 %             data.HRM.wfm_count  =   1:1:num_surfaces;
% % %                 
% % %                 
% % %                 % ----sigma0 scaling factor for conversion from Power to sigma0
% % %                 %units
% % %                 data.HRM.s0_sf=wfm_scaling_factor_sar_ku; %dB
% % %                 
% % %                 
% % %                 %--------------------------------------------------------------
% % %                 %------------- Stack characterization parameters --------------
% % %                 %--------------------------------------------------------------
% % %                 %----- Dopppler mask ------------------------------------------
% % %                 data.HRM.Doppler_mask   =   stack_mask_vector_start_stop.';
% % %                 data.HRM.pri_surf=pri_sar_pre_dat; %in seconds
% % %                 data.HRM.fs_clock_ku_surf=1./T0_sar_pre_dat; %in Hz
% % %                 
% % %                 %------------- Geometry-related parameters ------------------------
% % %                 if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
% % %                     switch cnf_p.looks_index_method
% % %                         case {'Doppler'}
% % %                             %This option Does not provide the correct information but it is kept for comparison purposes
% % %                             switch cnf_p.fd_method
% % %                                 case 'exact'
% % %                                     % exact method would require precise information for
% % %                                     % each beam within the stack (not currently developed for L1B products for SEOMs)
% % %                                     data.HRM.beam_ang_stack   =  beam_ang_surf.'; %TBD in radians
% % %                                     data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
% % %                                     data.HRM.V_stack          =  vel_norm_sat_beam_surf.';
% % %                                 case 'approximate'
% % %                                     data.HRM.pointing_ang_start_surf=start_pointing_angle; % in radians
% % %                                     data.HRM.pointing_ang_stop_surf=stop_pointing_angle; % in radians
% % %                                     data.HRM.doppler_angle_start_surf=start_doppler_angle; % in radians
% % %                                     data.HRM.doppler_angle_stop_surf=stop_doppler_angle; % in radians
% % %                             end
% % %                         case {'Look_angle'}
% % %                             switch cnf_p.look_ang_method
% % %                                 case 'exact'
% % %                                     % exact method would require precise information for
% % %                                     % each beam within the stack (not currently developed for L1B products for SEOMs)
% % %                                     data.HRM.look_ang_stack   =  look_ang_surf.'; % in radians
% % %                                     data.HRM.pri_stack        =  pri_sar_surf.'; %TBD %eventually for sentinel-6 the PRF could change
% % %                                     data.HRM.V_stack          =  vel_norm_sat_beam_surf.'; % TBD norm of velocity for each beam within stack
% % %                                 case 'approximate'
% % %                                     data.HRM.look_ang_start_surf=start_look_angle; % in radians
% % %                                     data.HRM.look_ang_stop_surf=stop_look_angle; % in radians
% % %                             end
% % %                             
% % %                     end
% % %                 end
% 
