function [SSH, SWH, sigma0, epoch, Pu, COR, filename_L2, filename_L1, index_inside_mask, flag]=read_L2(filesBulk, i_fileL2_input,  name_bs, cnf_tool, varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code reads geophysical parameters stored in L2 file and equivalent L1B file
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
%       -filesBulk{}    =   cell array of length num. of baselines of structures of input files within the folder where
%       to process the data (including the L1B as well as configuration/characterization files):
%        filesBulk{}.inputPath        --> full path to the input L1B folder
%        filesBulk{}.resultPath       --> full path to results folder save L2 prod
%        filesBulk{}.inputPath_L1BS   --> full path to input L1BS folder
%        filesBulk{}.cnf_chd_cst_path --> full path to cnf/chd/cst/LUTs folder
%        filesBulk{}.CNF_file         --> full filename cnf config file
%        filesBulk{}.CHD_file         --> full filename chd charac file
%        filesBulk{}.CST_file         --> full filename cst const file
%        filesBulk{}.LUT_f0_file      --> full filename LUT for f0 function
%        filesBulk{}.LUT_f1_file      --> full filename LUT for f1 function
%        filesBulk{}.nFilesL1B        --> total number of L1B to be processed
%        filesBulk{}.L1BFiles         --> information of the L1B files to be
%                                       processed: fields--> {name,date,bytes,isdir,datenum}
%       - i_fileL2_input = L2 file index to be processed
%       - name_bs = cell array containing baselines names

%      OPTIONAL
%       
% OUTPUT:
%       -
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% Geo_mask is not ready to use. It needs to be reviewed. Same for index_inside_mask
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:

%% --------------- HANDLING OPTIONAL INPUTS ---------------------------
if(nargin<3)
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('filename_mask_KML','',@(x)ischar(x));

p.parse(varargin{:});
filename_mask_KML=p.Results.filename_mask_KML;
clear p;

%% ------------ CHECKING ISR L2 AVAILABLE PRODUCT for other baselines------------   

% number of baselines to be compared
N_baselines = length(filesBulk);

empty_flag=ones(1,length(filesBulk)); % skip file - exit function - if file does not exist in any baseline (1)

i_baseline = 1;

filename_L2{i_baseline}=char(filesBulk(i_baseline).L2Files(i_fileL2_input).name); % no path
[~,aux,fileext_L2]=fileparts(filename_L2{i_baseline});
data_string=aux(21:21+30); % alba: take sensing time.
filename_L2{i_baseline}=strcat(char(filesBulk(i_baseline).input_path_L2_ISR_bs),filename_L2{i_baseline}); % with path

input_L1_ISR_Files=dir(fullfile(char(filesBulk(i_baseline).input_path_L1_ISR_bs),['*' data_string '*'])); 
if ~isempty(input_L1_ISR_Files)
    empty_flag(i_baseline)=0; % L1 file exists in baseline 1, 
else
    fprintf(char(strcat({'\n L2 file not available in L1 ISR data set '},char(name_bs(i_baseline)),': ',data_string,'\n')));
end

i_fileL1_input=find(~cellfun(@isempty,strfind({filesBulk(i_baseline).L1BFiles(:).name},input_L1_ISR_Files.name)), 1);
filename_L1{i_baseline}=char(filesBulk(i_baseline).L1BFiles(i_fileL1_input).name);
[~,aux,fileext_L1]=fileparts(filename_L1{i_baseline});
filename_L1{i_baseline}=strcat(char(filesBulk(i_baseline).input_path_L1_ISR_bs),filename_L1{i_baseline});

for i_baseline=2:N_baselines % find equivalent product for all baseline based on sensing time string

    input_L2_ISR_Files   = dir(fullfile(char(filesBulk(i_baseline).input_path_L2_ISR_bs),['*' data_string  strcat('*', fileext_L2)]));
    if ~isempty(input_L2_ISR_Files)
        empty_flag(i_baseline)=0;
    else
        empty_flag(i_baseline)=1;
        fprintf(char(strcat({'\n File not available in L2 ISR data set '},char(name_bs(i_baseline)),': ',data_string,'\n')));
    end
    i_fileL2_input=find(~cellfun(@isempty,strfind({filesBulk(i_baseline).L2Files(:).name},input_L2_ISR_Files.name)), 1);
    filename_L2{i_baseline}=char(filesBulk(i_baseline).L2Files(i_fileL2_input).name);
    filename_L2{i_baseline}=strcat(char(filesBulk(i_baseline).input_path_L2_ISR_bs),filename_L2{i_baseline});

    % alba: same for L1B file
    input_L1_ISR_Files = dir(fullfile(char(filesBulk(i_baseline).input_path_L1_ISR_bs),['*' data_string strcat('*', fileext_L1)])); % '*.NC'])); % alba: .nc?
    if ~isempty(input_L1_ISR_Files)
        empty_flag(i_baseline)=0;
    else
        empty_flag(i_baseline)=1;
        fprintf(char(strcat({'\n File not available in L1 ISR data set '},char(name_bs(i_baseline)),': ',data_string,'\n')));
    end       
    i_fileL1_input=find(~cellfun(@isempty,strfind({filesBulk(i_baseline).L1BFiles(:).name},input_L1_ISR_Files.name)), 1);
    filename_L1{i_baseline}=char(filesBulk(i_baseline).L1BFiles(i_fileL1_input).name);
    filename_L1{i_baseline}=strcat(char(filesBulk(i_baseline).input_path_L1_ISR_bs),filename_L1{i_baseline});

end

% fprintf('Total number of L2 files (%s) to be processed: %.0f\n', char(name_bs(i_baseline)), filesBulk(i_baseline).nFilesL2);


flag = 1;
% skip function if file does not exist in any baseline
if any(empty_flag)
    fprintf(char(strcat({'\n Missing files in baseline '},char(name_bs(find(empty_flag,1))),': ',data_string,'\n')));
    flag = 0;
    return;
end


%% ------------ LOAD THE GEOGRAPHICAL MASK --------------------------------
if isempty(filename_mask_KML)
    geo_mask  = [];
else
    geo_mask  = kml2lla(filename_mask_KML);
end


for i_baseline=1:N_baselines
    
    mode = 'HR';
    if any(contains(filename_L2{i_baseline},'LROS'))
        mode = 'LR';
    elseif any(contains(filename_L2{i_baseline},'LR'))
        mode = 'LR';
    end
    
    switch cnf_tool.L2proc{i_baseline}
        case {'GPP'}
            ncid = netcdf.open(filename_L2{i_baseline},'NOWRITE');            
            L2_num_surfaces = length(ncread(filename_L2{i_baseline},'data_20/ku/time'));
            lat_surf_L2=double(ncread(filename_L2{i_baseline},'data_20/ku/latitude')).';
            lon_surf_L2=double(ncread(filename_L2{i_baseline},'data_20/ku/longitude')).';
            if ~isempty(geo_mask)
                idx_inside_mask=inpolygon(lon_surf_L2,lat_surf_L2,geo_mask.coord(:,1),geo_mask.coord(:,2));
            else
                idx_inside_mask=logical(ones(1,L2_num_surfaces));
            end
            
            lat_surf_L2 = lat_surf_L2(idx_inside_mask);
            index_inside_mask{i_baseline} = idx_inside_mask;

%             SSH_dumm=double(ncread(filename_L2{i_baseline},'data_20/ku/ssha')).'; % NaNs for HR EUMETSAT!! 'Sea surface height anomaly = Altitude of satellite (altitude) - Ku band corrected ocean altimeter range (range_ocean) - geo corrections (ocean_geo_corrections) - mean sea surface (mean_sea_surface_sol1)'
            range_ocean=double(ncread(filename_L2{i_baseline},'data_20/ku/range_ocean')).'; % '20Hz_Ku: Range computed from MLE4 retracking (LR, LR-OS) or SAMOSA retracking (HR,LR-RMC). It includes all instrumental corrections from net_instr_cor_range_ocean'
            epoch_ocean=double(ncread(filename_L2{i_baseline},'data_20/ku/epoch_ocean')).'; % '2-way epoch derived from the ocean retracker'
            altitude=double(ncread(filename_L2{i_baseline},'data_20/ku/altitude')).';
            ocean_geo_corrections=double(ncread(filename_L2{i_baseline},'data_20/ku/ocean_geo_corrections')).'; % NaNs!!
            mean_sea_surface_sol1=double(ncread(filename_L2{i_baseline},'data_20/ku/mean_sea_surface_sol1')).'; % 'mean sea surface height solution 1 (CNES-CLS15) above WGS84 ellipsoid' 
            SSH_dumm=altitude - range_ocean - mean_sea_surface_sol1; % skip geo_corr -> NaNs
            SWH_dumm=double(ncread(filename_L2{i_baseline},'data_20/ku/swh_ocean')).';
            sigma0_dumm=double(ncread(filename_L2{i_baseline},'data_20/ku/sig0_ocean')).';  
            epoch_dumm=double(ncread(filename_L2{i_baseline},'data_20/ku/range_ocean')).'; 
%             COR_dumm=double(ncread(filename_L2{i_baseline},'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).';
%             Pu_dumm=10.^(double(ncread(filename_L2{i_baseline},'Pu_analytical_SWH_MSSfixed_20_ku')).'./10);

            netcdf.close(ncid);
            clear ncid;
            
            SSH{i_baseline}=SSH_dumm(idx_inside_mask);
            SWH{i_baseline}=SWH_dumm(idx_inside_mask);
            sigma0{i_baseline}=sigma0_dumm(idx_inside_mask);
            epoch{i_baseline}=epoch_dumm(idx_inside_mask);
            COR{i_baseline}=[];
%         %     Misfit{i_baseline}=Misfit_dumm(idx_inside_mask);
            Pu{i_baseline}=[]; 

            clear SSH_dumm SWH_dumm sigma0_dumm COR_dumm nb_dumm % Misfit_dumm;
            
            fprintf('Num surfaces L2 (%s): %d\n', char(name_bs(i_baseline)), length(lat_surf_L2));
            
        case {'ISD'}
               %---------------- read retrieved parameters: L2 file ------------------
            ncid = netcdf.open(filename_L2{i_baseline},'NOWRITE');
            dimid = netcdf.inqDimID(ncid,'time_20_ku');
            [~, L2_num_surfaces] = netcdf.inqDim(ncid,dimid);              
            lat_surf_L2=double(ncread(filename_L2{i_baseline},'lat_20_ku')).';
            lon_surf_L2=double(ncread(filename_L2{i_baseline},'lon_20_ku')).';
            if ~isempty(geo_mask)
                idx_inside_mask=inpolygon(lon_surf_L2,lat_surf_L2,geo_mask.coord(:,1),geo_mask.coord(:,2));
            else
                idx_inside_mask=logical(ones(1,L2_num_surfaces));
            end

            %Assuming they should have the same surfaces b1 and b2
            lat_surf_L2 = lat_surf_L2(idx_inside_mask);

            index_inside_mask{i_baseline} = idx_inside_mask;

            if strcmp(mode, 'LR')
                SSH_dumm=double(ncread(filename_L2{i_baseline},'ssh_LR_20_ku')).';
                SWH_dumm=double(ncread(filename_L2{i_baseline},'swh_LR_20_ku')).';
                sigma0_dumm=double(ncread(filename_L2{i_baseline},'sig0_LR_20_ku')).';  
                epoch_dumm=double(ncread(filename_L2{i_baseline},'epoch_LR_20_ku')).'; 
                COR_dumm=double(ncread(filename_L2{i_baseline},'Pearson_corr_Pu_LR_20_ku')).';
                Pu_dumm=10.^(double(ncread(filename_L2{i_baseline},'Pu_LR_20_ku')).'./10);
            else
                SSH_dumm=double(ncread(filename_L2{i_baseline},'ssh_analytical_SWH_MSSfixed_20_ku')).';
                SWH_dumm=double(ncread(filename_L2{i_baseline},'swh_analytical_SWH_MSSfixed_20_ku')).';
                sigma0_dumm=double(ncread(filename_L2{i_baseline},'sig0_analytical_SWH_MSSfixed_20_ku')).';  
                epoch_dumm=double(ncread(filename_L2{i_baseline},'epoch_analytical_SWH_MSSfixed_20_ku')).'; 
                COR_dumm=double(ncread(filename_L2{i_baseline},'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).';
                Pu_dumm=10.^(double(ncread(filename_L2{i_baseline},'Pu_analytical_SWH_MSSfixed_20_ku')).'./10);
            end

            netcdf.close(ncid);
            clear ncid;

            SSH{i_baseline}=SSH_dumm(idx_inside_mask);
            SWH{i_baseline}=SWH_dumm(idx_inside_mask);
            sigma0{i_baseline}=sigma0_dumm(idx_inside_mask);
            epoch{i_baseline}=epoch_dumm(idx_inside_mask);
            COR{i_baseline}=COR_dumm(idx_inside_mask);
        %     Misfit{i_baseline}=Misfit_dumm(idx_inside_mask);
            Pu{i_baseline}=Pu_dumm(idx_inside_mask); 

            clear SSH_dumm SWH_dumm sigma0_dumm COR_dumm nb_dumm % Misfit_dumm;

            %------------ read waveforms information: L1B file ----------------------
            L1B_num_surfaces = length(ncread(filename_L1{i_baseline},'data_20/ku/time_tai'));
            lat_surf=ncread(filename_L1{i_baseline},'data_20/ku/latitude').';
            lon_surf=wrapTo180(ncread(filename_L1{i_baseline},'data_20/ku/longitude').');
            if ~isempty(geo_mask)
                idx_inside_mask_L1B=inpolygon(lon_surf,lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
            else
                idx_inside_mask_L1B=logical(ones(1,L1B_num_surfaces));
            end
            first_indx_L1B = find(idx_inside_mask_L1B==1,1,'first');                            
            [~,dumm_indx] = min(abs(lat_surf(idx_inside_mask_L1B)-lat_surf_L2(1)));
            idx_inside_mask_L1B(1:(first_indx_L1B+dumm_indx-2))= 0; %from beginning up to the surface where coincides
            [~,dumm_indx] = min(abs(lat_surf-lat_surf_L2(end)));
            idx_inside_mask_L1B(dumm_indx+1:end)= 0; %from beginning up to the surface where coincides

            fprintf('Num surfaces L2 (%s): %d\n', char(name_bs(i_baseline)), length(lat_surf_L2));
            fprintf('Num surfaces L1 (%s): %d\n', char(name_bs(i_baseline)), length(find(idx_inside_mask_L1B)));
            
    end

end