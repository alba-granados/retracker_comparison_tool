function filesBulk = split_bsl_tracks_into_pass(filesBulk)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code splits all tracks contained in a baseline - filesBulk - into
% passes in order to avoid duplicated latitudes to perform retrieval
% statistics
%
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 13/12/2020
% This software is built within the Sentinel-6 P4 L1 GPP project 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%         - filesBulk{}    =   cell array of length num. of baselines of structures of input files within the folder where
%       to process the data (including the L1B as well as configuration/characterization files):
%         - filesBulk(i_baseline).input_path_L1_ISR_bs
%         - filesBulk(i_baseline).nFilesL1B
%         - filesBulk(i_baseline).L1BFiles
%         - filesBulk(i_baseline).input_path_L2_ISR_bs
%         - filesBulk(i_baseline).nFilesL2
%         - filesBulk(i_baseline).L2Files
% 
% OUTPUT:
%       - filesBulk{} same structure with passes
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% netcdfSplit.m by Ferran Gibert on 20/2/2020
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:

indx_delete_files_duplicated_lats = [];

for ii=1:filesBulk.nFilesL1B
    lat = ncread([filesBulk.input_path_L1_ISR_bs, filesBulk.L1BFiles(ii).name], 'data_20/ku/latitude').';
    [~, indx_min_lat] = min(lat); 
    [~, indx_max_lat] = max(lat);
    
    if  ~isempty(find(max(indx_min_lat, indx_max_lat) > 1 && max(indx_min_lat, indx_max_lat) < length(lat), 1)) % check whether max/min latitude is not at beginning/end of sensing = track contains more than one pass
        % tracks has duplicated latitudes (contains > 1 pass)
        
        indx_delete_files_duplicated_lats = [indx_delete_files_duplicated_lats ii];
        
        [~,filenopath_L1B,file_ext]=fileparts(filesBulk.L1BFiles(ii).name);
        original_sensing_time=filenopath_L1B(21:21+30); 

        % find equivalent L2 file
        input_L2_ISR_Files   = dir(fullfile(char(filesBulk.input_path_L2_ISR_bs),['*' original_sensing_time '*']));
        if isempty(input_L2_ISR_Files)
            fprintf('\n File not available in L2 ISR data set %s\n',original_sensing_time);
        end
        i_fileL2_input=find(~cellfun(@isempty,strfind({filesBulk.L2Files(:).name},input_L2_ISR_Files.name)), 1);   
        [~,filenopath_L2,file_ext_L2]=fileparts(filesBulk.L2Files(i_fileL2_input).name);

        max_min_lat_indx = max(indx_min_lat, indx_max_lat);
      
        % create first part of L1B and equivalent L2 netcdf (first pass) with new sensing time
        interval_to_select = [1 max_min_lat_indx-1];

        [file_L1B_new]=netcdfSplit_linux(filesBulk.input_path_L1_ISR_bs, filenopath_L1B, interval_to_select, 'data_20/ku/time', file_ext);      
        [file_L2_new]=netcdfSplit_linux(filesBulk.input_path_L2_ISR_bs, filenopath_L2, interval_to_select, 'time_20_ku', file_ext_L2);    
        
        % update L1B/L2 list of products to be compared 
        filesBulk.L1BFiles(end+1) = dir(file_L1B_new);
        filesBulk.L2Files(end+1) = dir(file_L2_new);
        
        % create second part of L1B and equivalent L2 netcdf (first pass) with new sensing time
        interval_to_select = [max_min_lat_indx length(lat)-1];

        [file_L1B_new]=netcdfSplit_linux(filesBulk.input_path_L1_ISR_bs, filenopath_L1B, interval_to_select, 'data_20/ku/time', file_ext);      
        [file_L2_new]=netcdfSplit_linux(filesBulk.input_path_L2_ISR_bs, filenopath_L2, interval_to_select, 'time_20_ku', file_ext_L2);  
        
        % update L1B/L2 list of products to be compared 
        filesBulk.L1BFiles(end+1) = dir(file_L1B_new);
        filesBulk.L2Files(end+1) = dir(file_L2_new);
        
        filesBulk.nFilesL1B = filesBulk.nFilesL1B + 2;
        filesBulk.nFilesL2 = filesBulk.nFilesL2 + 2;
        
    end
    
end

for ii=1:numel(indx_delete_files_duplicated_lats)
    filesBulk.L1BFiles(indx_delete_files_duplicated_lats(ii)) = [];
    filesBulk.L2Files(indx_delete_files_duplicated_lats(ii)) = [];
    filesBulk.nFilesL1B=filesBulk.nFilesL1B-1;
    filesBulk.nFilesL2=filesBulk.nFilesL2-1;
end
    
end