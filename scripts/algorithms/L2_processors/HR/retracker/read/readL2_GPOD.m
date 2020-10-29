function data = readL2_GPOD (filename_L2_GPOD,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading L2 data from GPOD
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
%           -filename_L2_ISR    =   L1B filename with the fullpath information
%           -cnf_p           =   Structure with configuration parameters of the
%                                L2 processor
%       OPTIONAL:
%           -
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
if(nargin<2 || nargin>(2+1*2))
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
[~,name,ext]=fileparts(filename_L2_GPOD);
ext=lower(ext);
switch ext
    case '.nc'
        %% --------------------- NetCDF ---------------------------------------
	
        % Based on the L1B product generated for DeDop a la Sentinel-3
        % open the netCDF file
        ncid=netcdf.open(filename_L2_GPOD,'NC_NOWRITE');
        
        % Get dimensions
        dimid=netcdf.inqDimID(ncid,'Meas_Index_20Hz');
        [~,num_surfaces]=netcdf.inqDim(ncid,dimid);
        netcdf.close(ncid);
        

        data.N_records=num_surfaces;
        
        % -----------------------------------------------------------------
        % GEO: geographical information
        % -----------------------------------------------------------------
        data.GEO.TAI.total             =   ncread(filename_L2_GPOD,'TAI_Time_20Hz').';
        time_1Hz                       =   ncread(filename_L2_GPOD,'TAI_Time_1Hz').';
        data.GEO.TAI.days              =   ceil(data.GEO.TAI.total/cst_p.sec_in_day_cst);
        data.GEO.TAI.secs              =   ceil(data.GEO.TAI.total-data.GEO.TAI.days*cst_p.sec_in_day_cst);
        data.GEO.TAI.microsecs         =   ((data.GEO.TAI.total-data.GEO.TAI.secs)*1e6);
        data.GEO.LAT                   =   ncread(filename_L2_GPOD,'latitude_20Hz').';
        data.GEO.LON                   =   wrapTo180(ncread(filename_L2_GPOD,'longitude_20Hz').'); %defined between [0,360]                
        data.GEO.H                     =   ncread(filename_L2_GPOD,'altitude_20Hz').';
        data.GEO.H_rate                =   ncread(filename_L2_GPOD,'altitude_rate_20Hz').'; % m/s
        data.GEO.V                     =   ncread(filename_L2_GPOD,'satellite_velocity_20Hz').'; % m/s
        data.GEO.pitch                 =   ncread(filename_L2_GPOD,'pitch_mispointing_20Hz').'; 
        data.GEO.roll                  =   ncread(filename_L2_GPOD,'roll_mispointing_20Hz').'; 
        data.GEO.yaw                   =   zeros(1,data.N_records); 
        % -----------------------------------------------------------------
        % MEA: measurements
        % -----------------------------------------------------------------
        data.MEA.win_delay = double(ncread(filename_L2_GPOD,'Window_Delay_20Hz'));
        
        % ---------------------------------------------------------------------
        % COR: Geophysical corrections
        % ---------------------------------------------------------------------
        dry_trop                   =   ncread(filename_L2_GPOD,'Dry_Corr_1Hz').';
        wet_trop                   =   ncread(filename_L2_GPOD,'Wet_Corr_1Hz').';
        inv_bar                    =   ncread(filename_L2_GPOD,'IB_Corr_1Hz').';
        dac                        =   ncread(filename_L2_GPOD,'DAC_Corr_1Hz').';
        gim_ion                    =   ncread(filename_L2_GPOD,'GIM_Corr_1Hz').';
        model_ion                  =   zeros(1,length(time_1Hz));
        ocean_equilibrium_tide     =   ncread(filename_L2_GPOD,'OET_Corr_1Hz').';
        ocean_longperiod_tide      =   ncread(filename_L2_GPOD,'OLPT_Corr_1Hz').';
        ocean_loading_tide         =   ncread(filename_L2_GPOD,'OLT_Corr_1Hz').';
        solidearth_tide            =   ncread(filename_L2_GPOD,'SET_Corr_1Hz').';
        geocentric_polar_tide      =   ncread(filename_L2_GPOD,'GPT_Corr_1Hz').';
        data.COR.total_ocean_applied       =   ncread(filename_L2_GPOD,'GEO_Corr_20Hz').';
        data.COR.total_land_applied       =   ncread(filename_L2_GPOD,'GEO_Corr_Land_20Hz').';
        
        %at 20-Hz using closest time
        for i_record=1:data.N_records
            [~,idx_time_closest]=min(abs(data.GEO.TAI.total(i_record)-time_1Hz));
            data.COR.dry_trop(i_record)         =   dry_trop(idx_time_closest(1));
            data.COR.wet_trop(i_record)         =   wet_trop(idx_time_closest(1));
            data.COR.inv_bar(i_record)          =   inv_bar(idx_time_closest(1));
            data.COR.dac(i_record)              =   dac(idx_time_closest(1));
            data.COR.gim_ion(i_record)          =   gim_ion(idx_time_closest(1));
            data.COR.model_ion(i_record)        =   model_ion(idx_time_closest(1));
            data.COR.ocean_equilibrium_tide(i_record)     =   ocean_equilibrium_tide(idx_time_closest(1));
            data.COR.ocean_longperiod_tide(i_record)      =   ocean_longperiod_tide(idx_time_closest(1));
            data.COR.ocean_loading_tide(i_record)         =   ocean_loading_tide(idx_time_closest(1));
            data.COR.solidearth_tide(i_record)            =   solidearth_tide(idx_time_closest(1));
            data.COR.geocentric_polar_tide(i_record)      =   geocentric_polar_tide(idx_time_closest(1));
        end
        %---------- Combined corrections --------------------------------------
        data.COR.prop_GIM_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.gim_ion; % not clear??
        data.COR.prop_Mod_ion   =   data.COR.dry_trop + data.COR.wet_trop + data.COR.model_ion; % not clear ??
        data.COR.surf_dac       =   data.COR.dac; % SSB shoudl be added according to Cristina's code??
        data.COR.surf_invb      =   data.COR.inv_bar; % SSB shoudl be added according to Cristina's code??
        data.COR.geop           =   data.COR.ocean_equilibrium_tide + data.COR.ocean_longperiod_tide...
            + data.COR.ocean_loading_tide + data.COR.solidearth_tide + data.COR.geocentric_polar_tide;

        %---------------------------------------------------------------
        % SWH/SIGMA0 instrumental correction
        %---------------------------------------------------------------
        data.MOD_INSTR_CORR.SWH                   =   zeros(1,data.N_records);
        data.MOD_INSTR_CORR.sig0                  =   zeros(1,data.N_records);
        
        % -----------------------------------------------------------------
        % -------------- Geophysical retrievals ---------------------------
        % -----------------------------------------------------------------
        data.GR.analytical_SWH.SSH   = double(ncread(filename_L2_GPOD,'SSH_20Hz')).';       
        data.GR.analytical_SWH.retracked_range = double(ncread(filename_L2_GPOD,'Range_Unc_20Hz')).'; %[seconds]
        data.GR.analytical_SWH.SWH   = double(ncread(filename_L2_GPOD,'SWH_20Hz')).';
        data.GR.analytical_SWH.sig0  = double(ncread(filename_L2_GPOD,'Sigma0_20Hz')).';
        data.GR.analytical_SWH.misfit   = double(ncread(filename_L2_GPOD,'Misfit_20Hz')).';
        
    otherwise
        error(strcat('File extension ',cnf_p.mission,' is not currently contemplated or not valid'));
end

end


