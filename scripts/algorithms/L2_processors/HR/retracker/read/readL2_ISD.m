function data = readL2_ISD (filename_L2_ISR,retracker_ID,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading L2 data from isardSAT output retracker in
% netcdf or .mat (to be implemneted)
%
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 15/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       MANDATORY:
%           -filename_L2_ISR    =   L1B filename with the fullpath information
%           -cnf_p           =   Structure with configuration parameters of the
%                                L2 processor
%       OPTIONAL:
%           -filename_L2_ISRS   =   L1B-S filename with the fullpath name
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
% v1.0: Only the ntecdf files can be ingested

%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+1*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('sh_name_nc',{'ssh'},@(x)ischar(x));
p.parse(varargin{:});
sh_name_nc=char(p.Results.sh_name_nc);
clear p;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
[~,name,ext]=fileparts(filename_L2_ISR);
ext=lower(ext);
switch ext
    case '.nc'
        %% --------------------- NetCDF ---------------------------------------
	
        % Based on the L1B product generated for DeDop a la Sentinel-3
        % open the netCDF file
        ncid=netcdf.open(filename_L2_ISR,'NC_NOWRITE');
        
        % Get dimensions
        dimid=netcdf.inqDimID(ncid,'time_20_ku');
        [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
        netcdf.close(ncid);
        

        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------
        data.GEO.TAI.total             =   ncread(filename_L2_ISR,'time_20_ku').';
        data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
        data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
        data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
        data.GEO.LAT                   =   ncread(filename_L2_ISR,'lat_20_ku').';
        data.GEO.LON                   =   wrapTo180(ncread(filename_L2_ISR,'lon_20_ku').'); %defined between [0,360]                
        data.GEO.H                     =   ncread(filename_L2_ISR,'alt_20_ku').';
        
        
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = double(ncread(filename_L2_ISR,'range_20_ku')).'*2.0/cst_p.c_cst;
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
        try
            data.COR.dry_trop                   =   double(ncread(filename_L2_ISR,'dry_tropo_correction_20_ku')).';
            data.COR.wet_trop                   =   double(ncread(filename_L2_ISR,'wet_tropo_correction_20_ku')).';
            data.COR.inv_bar                    =   double(ncread(filename_L2_ISR,'inverse_baro_correction_20_ku')).';
            data.COR.dac                        =   double(ncread(filename_L2_ISR,'Dynamic_atmospheric_correction_20_ku')).';
            data.COR.gim_ion                    =   double(ncread(filename_L2_ISR,'GIM_iono_correction_20_ku')).';
            data.COR.model_ion                  =   double(ncread(filename_L2_ISR,'model_iono_correction_20_ku')).';
            data.COR.ocean_equilibrium_tide     =   double(ncread(filename_L2_ISR,'ocean_equilibrium_tide_20_ku')).';
            data.COR.ocean_longperiod_tide      =   double(ncread(filename_L2_ISR,'long_period_tide_20_ku')).';
            data.COR.ocean_loading_tide         =   double(ncread(filename_L2_ISR,'ocean_loading_tide_20_ku')).';
            data.COR.solidearth_tide            =   double(ncread(filename_L2_ISR,'solid_earth_tide_20_ku')).';
            data.COR.geocentric_polar_tide      =   double(ncread(filename_L2_ISR,'geocentric_polar_tide_20_ku')).';
            %---------- Combined corrections --------------------------------------
            data.COR.prop_GIM_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.gim_ion; % not clear??
            data.COR.prop_Mod_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.model_ion; % not clear ??
            data.COR.surf_dac       =   data.COR.dac; % SSB shoudl be added according to Cristina's code??
            data.COR.surf_invb      =   data.COR.inv_bar; % SSB shoudl be added according to Cristina's code??
            data.COR.geop           =   data.COR.ocean_equilibrium_tide + data.COR.ocean_longperiod_tide...
                + data.COR.ocean_loading_tide + data.COR.solidearth_tide + data.COR.geocentric_polar_tide;
        catch 
            disp('No geophysical corrections in the L2')
        end

        %---------------------------------------------------------------
        % SWH/SIGMA0 instrumental correction
        %---------------------------------------------------------------
        data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
        data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
        
        % -----------------------------------------------------------------
        % -------------- Geophysical retrievals ---------------------------
        % -----------------------------------------------------------------
        N_retrackers=length(retracker_ID);
        for i_retracker=1:N_retrackers
            if ~isempty(strfind(lower(char(retracker_ID(i_retracker))),'2-step'))
                data.GR.analytical_2step.SSH = double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_analytical_2step_20_ku'))).';
                data.GR.analytical_2step.epoch = double(ncread(filename_L2_ISR,'epoch_analytical_2step_20_ku')).'; % [seconds]
                data.GR.analytical_2step.retracked_range = double(ncread(filename_L2_ISR,'retracked_range_analytical_2step_20_ku')).'; %[seconds]
                data.GR.analytical_2step.SWH = double(ncread(filename_L2_ISR,'swh_analytical_2step_20_ku')).';
                data.GR.analytical_2step.sig0= double(ncread(filename_L2_ISR,'sig0_analytical_2step_20_ku')).';
                data.GR.analytical_2step.COR = double(ncread(filename_L2_ISR,'Pearson_corr_analytical_2step_20_ku')).';
                try
                    data.GR.analytical_2step.misfit = double(ncread(filename_L2_ISR,'Misfit_analytical_2step_20_ku')).';
                catch
                    disp('No misfit information in the L2');
                end
            elseif ~isempty(strfind(lower(char(retracker_ID(i_retracker))),'mss'))
                data.GR.analytical_MSS.SSH   = double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_analytical_MSS_SWHfixed_20_ku'))).';
                data.GR.analytical_MSS.epoch = double(ncread(filename_L2_ISR,'epoch_analytical_MSS_SWHfixed_20_ku')).'; % [seconds]
                data.GR.analytical_MSS.retracked_range = double(ncread(filename_L2_ISR,'retracked_range_analytical_MSS_SWHfixed_20_ku')).'; %[seconds]
                data.GR.analytical_MSS.sig0  = double(ncread(filename_L2_ISR,'sig0_analytical_MSS_SWHfixed_20_ku')).';
                data.GR.analytical_MSS.COR   = double(ncread(filename_L2_ISR,'Pearson_corr_analytical_MSS_SWHfixed_20_ku')).';
                data.GR.analytical_MSS.SWH   = NaN(1,length(data.GR.analytical_MSS.SSH));
                try
                    data.GR.analytical_MSS.misfit   = double(ncread(filename_L2_ISR,'Misfit_analytical_MSS_SWHfixed_20_ku')).';
                catch
                    disp('No misfit information in the L2');
                end
            elseif ~isempty(strfind(lower(char(retracker_ID(i_retracker))),'swh'))
                data.GR.analytical_SWH.SSH   = double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_analytical_SWH_MSSfixed_20_ku'))).';
                data.GR.analytical_SWH.epoch = double(ncread(filename_L2_ISR,'epoch_analytical_SWH_MSSfixed_20_ku')).'; % [seconds]
                data.GR.analytical_SWH.retracked_range = double(ncread(filename_L2_ISR,'retracked_range_analytical_SWH_MSSfixed_20_ku')).'; %[seconds]
                data.GR.analytical_SWH.SWH   = double(ncread(filename_L2_ISR,'swh_analytical_SWH_MSSfixed_20_ku')).';
                data.GR.analytical_SWH.sig0  = double(ncread(filename_L2_ISR,'sig0_analytical_SWH_MSSfixed_20_ku')).';
                data.GR.analytical_SWH.COR   = double(ncread(filename_L2_ISR,'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).';
                try
                    data.GR.analytical_SWH.misfit   = double(ncread(filename_L2_ISR,'Misfit_analytical_SWH_MSSfixed_20_ku')).';
                catch
                    disp('No misfit information in the L2');
                end
            elseif ~isempty(strfind(lower(char(retracker_ID(i_retracker))),'threshold'))
                data.GR.threshold.SSH   = double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_threshold_20_ku'))).';
                data.GR.threshold.epoch = double(ncread(filename_L2_ISR,'epoch_threshold_20_ku')).'; % [seconds]
                data.GR.threshold.retracked_range = double(ncread(filename_L2_ISR,'retracked_range_threshold_20_ku')).'; %[seconds]
                data.GR.threshold.SWH   = NaN(1,length(data.GR.threshold.SSH));
                data.GR.threshold.sig0  = double(ncread(filename_L2_ISR,'sig0_threshold_20_ku')).';
                data.GR.threshold.COR   = NaN(1,length(data.GR.threshold.SSH));                    
            elseif ~isempty(strfind(lower(char(retracker_ID(i_retracker))),'ocog'))
                data.GR.OCOG.SSH   = double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_OCOG_20_ku'))).';
                data.GR.OCOG.epoch = double(ncread(filename_L2_ISR,'epoch_OCOG_20_ku')).'; % [seconds]
                data.GR.OCOG.retracked_range = double(ncread(filename_L2_ISR,'retracked_range_OCOG_20_ku')).'; %[seconds]
                data.GR.OCOG.SWH   = NaN(1,length(data.GR.OCOG.SSH));
                data.GR.OCOG.sig0  = double(ncread(filename_L2_ISR,'sig0_OCOG_20_ku')).';
                data.GR.OCOG.COR   = NaN(1,length(data.GR.OCOG.SSH));                                    
            else
                disp(strcat('No valid retracker ',char(retracker_ID(i_retracker))));
            end
        end
            
    case '.mat'
        %% ------------------------ Matlab data -------------------------------
        
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end


