function output_data_generation_L2_LRM(file,retrackers_results,data,cnf_p,cst_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and generates the output L2 product
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         --- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 07/09/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -file    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -retrackers_results = structure of the fitting procedures for different retrackers and 
%                           different records
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% -prepare_NetCDF_L2_LRM: generate the L2 product with the necessary and
% required output product information on a netcdf file
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
% Code adapted from /retracker/algorithms/LRM/output_data_generation_L2.m
% by Eduard Makhoul
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:


%% -------------- ORGANIZE THE OUTPUT DATA --------------------------------
%--------------------------------------------------------------------------
%strcture to be used for .mat validation
out.N_records   =   data.N_records; 
out.TAI         =   data.GEO.TAI; %seconds since first janaury 2000, days and seconds in a day structure
out.lat         =   data.GEO.LAT;
out.lon         =   data.GEO.LON;
out.H_orb       =   data.GEO.H; %include the height information
out.range       =   data.MEA.win_delay*cst_p.c_cst/2.0;
out.s0_sf       =   data.HRM.s0_sf; %sigma scale factor copied from L1B in dBs

%--------------------------------------------------------------------------
%------------------------- RETRACKERS OUTPUTS -----------------------------
%--------------------------------------------------------------------------
for i_retracker=1: length(cnf_p.retracker_name)    
    switch char(cnf_p.retracker_name(i_retracker))
        case {'LR'}
            out.RETRACKER.LR.Epoch           =   retrackers_results.LR.Epoch;
            out.RETRACKER.LR.Pu              =   retrackers_results.LR.Pu;   
            out.RETRACKER.LR.SSH             =   retrackers_results.LR.SSH;     
            out.RETRACKER.LR.sigma0          =   retrackers_results.LR.sigma0;
            out.RETRACKER.LR.COR          =   retrackers_results.LR.COR.*1e2;
            out.RETRACKER.LR.Hs          =   retrackers_results.LR.Hs;
    end
    
end

% -------------------------------------------------------------------------
% ------------------------ PROCESSING OPTIONS -----------------------------
% -------------------------------------------------------------------------
out.PROC_CNF=cnf_p;

% -------------------------------------------------------------------------
% ------------------------ GLOBAL ATTRIBUTES ------------------------------
% -------------------------------------------------------------------------
% Imported directly from L1B: Including both data file related attributes &
% orbital info
if isfield(data,'GLOBAL_ATT')
    out.GLOBAL_ATT=data.GLOBAL_ATT;
    %include the original L1B file used
    out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=[file.filename_L1B_nopath file.fileext_L1B];    
end

%% --------------- WRITE OUTPUT PRODUCT -----------------------------------
%--------------------------------------------------------------------------
if cnf_p.write_output
    
    %-------------- Define the output name for the output file ----------------
    date_creation = datestr(now, '_yyyymmddTHHMMSS_');
    aux=strsplit(date_creation,'_');
    out.date_creation=aux(2);
    clear aux;

    name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','L2'),'_isd');

    if cnf_p.optional_ext_file_flag
        %name_L2_product=strrep(name_L2_product,'isd',strcat(cnf_p.file_ext_string,'_isd'));
        name_L2_product=strcat(name_L2_product,'_',cnf_p.file_ext_string);
    end

    % ---------------------- GENERATION OF OUTPUT PRODUCT --------------------
    switch cnf_p.output_product_format
        case 'mat' % Matlab output
            % .mat file
            save([file.resultPath 'data' filesep name_L2_product '.mat'],'out');        
        case 'nc'  % NetCDF output product
            prepare_NetCDF_L2_LRM([file.resultPath 'data' filesep name_L2_product '.nc'],out,cnf_p); %can have a switch to define the variables if mission is S3 or S6
    end

end


end

