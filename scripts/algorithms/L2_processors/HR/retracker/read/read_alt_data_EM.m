    function [data,exit_flag] = read_alt_data_EM (filename_L1B, cnf_p,...
                                                            cst_p,chd_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from any L1 processor we work
% with and store it in a common structure for the L2 processor
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 13/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       filename_L1B    =   L1B filename with the fullpath information
%       cnf_p       =   configuration parameters
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%       exit_flag   =   flag indicating whether processing succesful 1 or
%       there is an error -1
%  
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - readL1B_CS2_ESA: read Cryosat-2 L1B products produced by ESA 
%                    (currently only developed for .DBL, takes into account netCDF .nc in future)
% - readL1B_CS2_ISD: read Cryosat-2 L1B products produced by ISD 
%                    (currently developed for netCDF .nc files using a la Sentinel-3 format)
% - readL1B_S3_ISD: read Sentinel-3 L1B products produced by ISD 
%                    (exactly the same as the one for CS-2 a la Sentinel-3 format)
% - readL1B_S6_ISD: read Sentinel-6 L1B products produced by ISD 
%                    (reading either the .mat files or the .nc with the final format for Sentinel-6 using PSD issue 1.3)
% - filter_L1B_data: filter the data by regions of interest and_/or number
%                    of looks within the stack 
% - extract_ref_height: extract the reference height (window delay) from a DEM for the surfaces beneath satellite 
%                       (used for locating the part of the waveform coming from nadir/water body) 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Based on read_alt_data.m defined different reading functions for
% the different missions and different potential processors 
% (based on the formating of the output products, .nc, DBL. and .mat)

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(4+11*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
p.addParamValue('filename_L2',{''},@(x)ischar(x));
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('DEM_ref','SRTM',@(x)ischar(x));
p.addParamValue('dir_DEM','.',@(x)ischar(x));
p.addParamValue('path_Results',{''},@(x)ischar(x));
p.addParamValue('input_path_L2_GPOD',{''}); 
p.addParamValue('attitude_extraction_GPOD',0);
p.addParamValue('file_unique_id',{''});
p.addParamValue('Sig0AtmCorr_path',{''});
p.addParamValue('filename_dist_to_coast',{''});

p.parse(varargin{:});
filename_L1BS=char(p.Results.filename_L1BS);
filename_L2=char(p.Results.filename_L2);
filename_mask_KML=char(p.Results.filename_mask_KML);
DEM_ref=p.Results.DEM_ref;
dir_DEM=p.Results.dir_DEM;
path_Results=char(p.Results.path_Results);
input_path_L2_GPOD  = p.Results.input_path_L2_GPOD;
attitude_extraction_GPOD = p.Results.attitude_extraction_GPOD;
file_unique_id = p.Results.file_unique_id;
Sig0AtmCorr_path = p.Results.Sig0AtmCorr_path;
filename_dist_to_coast = char(p.Results.filename_dist_to_coast);
clear p;
exit_flag=1;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
switch cnf_p.mission
    
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                [data]=readL1B_CS2_ESA(filename_L1B,cnf_p,cst_p,chd_p);
            case 'ISD'
                [data]=readL1B_CS2_ISD(filename_L1B,cnf_p,cst_p,chd_p,'filename_L1BS',filename_L1BS);
            case 'GPOD'
                [data]=readL1B_GPOD_ESA(filename_L1B,cnf_p,cst_p,chd_p);
        end
    case 'S3'
    %% --------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case {'ESA'}
                % TBD
                [data]=readL1B_S3_ESA(filename_L1B,cnf_p,cst_p,chd_p,'filename_L1BS',filename_L1BS,'filename_L2',filename_L2);
            case {'DeDop'}
                [data]=readL1B_S3_DeDop(filename_L1B,cnf_p,cst_p,chd_p,'filename_L1BS',filename_L1BS,'filename_L2',filename_L2);
            case 'ISD'
                [data]=readL1B_S3_ISD(filename_L1B,cnf_p,cst_p,chd_p,'filename_L1BS',filename_L1BS);
        end
    case {'S6','JCS'}
    %% --------------------- Sentinel-6 -----------------------------------        
        switch cnf_p.L1proc
            case 'ISD'
                [data]=readL1B_S6_ISD(filename_L1B,cnf_p,cst_p,'filename_L1BS',filename_L1BS);
        end 
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end


%% ------------ Reading Sigma-0 atmospheric attenuation correction---------
if cnf_p.atm_att_correction_flag
    try 
        [data] = readSig0AtmCorr_GPOD(Sig0AtmCorr_path,filename_L1B,data);
    catch
        %if there geophsyical correction map missing set to zero correction and
        %flag
        data.COR_sig0.Sig0AtmCorr = zeros(1,data.N_records);        
    end
else
    data.COR_sig0.Sig0AtmCorr = zeros(1,data.N_records);
end

% inputFile = dir(strcat(char(input_path_L2_GPOD),'*',file_unique_id,'*'));
% if ~isempty(inputFile)
%     filename_L2_GPOD = strcat(char(input_path_L2_GPOD),char(inputFile.name));
%     GPOD_lat_surf = ncread(filename_L2_GPOD,'latitude_20Hz').';  
%     GPOD_lon_surf = ncread(filename_L2_GPOD,'longitude_20Hz').';  
%     Delta_Sigma0_20Hz_GPOD = ncread(filename_L2_GPOD,'Delta_Sigma0_20Hz');
% else
%     disp(strcat('No input path for file with ID ',file_unique_id));
%     exit_flag = -1;
%     return;
% end
% 
% figure; plot(data.GEO.LAT,data.COR.Sig0AtmCorr);
% hold on; plot(GPOD_lat_surf,Delta_Sigma0_20Hz_GPOD);
% xlabel('Latitude [deg]'); ylabel('\Delta\sigma^0 [dB]'); title(strcat('Atmospheric attenuation correction:',file_unique_id),'Interpreter','none');
% legend('ISR','GPOD');
% print('-dpng ',[path_Results,'plots',filesep,file_unique_id,'_att_corr','.png']);
% 
% return;

%% -------------------- Filter data ---------------------------------------
%Bye geographic location using an external .kml or/and depending on the
%size of the associated stack or depending on whether land shall be
%filtered out or not (using land sea mask options)
[data,flag]=filter_L1B_data (data,cnf_p,'filename_mask_KML',filename_mask_KML,...
                            'attitude_extraction_GPOD',attitude_extraction_GPOD,...
                            'input_path_L2_GPOD',input_path_L2_GPOD,...
                            'file_unique_id',file_unique_id,...
                            'filename_dist_to_coast',filename_dist_to_coast);
if flag==-1
    %track is not within the limits of the 
    exit_flag=flag;
    return
end

%% ---------------- Seeding or portion selection --------------------------
if cnf_p.wvfm_portion_selec 
    switch cnf_p.wvfm_portion_selec_type
        case 'ref_height'
            %% ------------------- Load DEM -------------------------------------------
            %interpolate the reference DEM for the location of interest
            %read the DEM
            switch cnf_p.mission
                case {'S3','S3A','S3B'}
                    %% --------------------- Sentinel-3 -----------------------------------
                    switch cnf_p.L1proc
                        case 'ESA'
                            % TBD
                            dumm=strsplit(filename_L1B,{'\','/'});
                            L1B_filename=char(dumm(end-1));
                            clear dumm;
                    end
                otherwise
                    [~,L1B_filename,~]=fileparts(filename_L1B);
            end
            [data]=extract_ref_height(data,cnf_p,cst_p,'DEM_ref',DEM_ref,'dir_DEM',dir_DEM,'L1B_filename',L1B_filename,'path_Results',path_Results);
        case 'CP4O'
            [data]=seed_cp4o(data,cnf_p,chd_p,cst_p);
    end       
end
end

