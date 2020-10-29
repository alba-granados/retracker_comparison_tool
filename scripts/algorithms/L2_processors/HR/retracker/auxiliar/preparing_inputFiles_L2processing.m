function  [filename_L1B,filename_L1B_nopath,fileext_L1B,filename_L1BS,filename_L2,cnf_p,file_unique_id]=preparing_inputFiles_L2processing(filesBulk, i_fileL1B_input, cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code sets and prepares the input files definition and untar extraction for L1B
% and L1Bs file used in the L2 processing (used to have a more clear L2_processing.m function)
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 21/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filesBulk    =   structure of input files within the folder where
%       to process the data (including the L1B as well as configuration/characterization files)
%       -i_fileL1B_input = index of the L1B
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       -filename_L1B        =   full filename of the input L1B including
%                                the path
%       -filename_L1B_nopath =   filename of the input L1B with neither
%                                full path nor file extension
%       -fileext_L1B         =   file extension of the input L1B file
%       -filename_L1BS       =   full filename of the input L1BS including
%                                the path 
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
%  It is assumed that within each input folder either for the L1B or for
%  the LBS there is only a single file with the same unique_identifier 
% (inital-final measurement/validty time), i.e., it must be avoided to
% have same file with .nc or .DBL and the same compressed file .TGZ
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 
%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+3*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.parse(varargin{:});
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
clear p;

%% --------------------- Prepare the input files definition ---------------
%--------------------------------------------------------------------------
% ------------- L1B file --------------------------------------------------
% -------------------------------------------------------------------------
switch cnf_p.mission
    case {'S3A','S3B','S3'}
        switch cnf_p.L1proc
            case 'ESA'
                %a measurement netcdf file is included within the 
                %reading directly from L2 enhanced 
                if strfind(filesBulk.L1BFiles(i_fileL1B_input).name,'_SR_1')
                    filename_L1B=char([filesBulk.inputPath filesBulk.L1BFiles(i_fileL1B_input).name filesep 'measurement.nc']);
                    cnf_p.input_type =1;
                elseif strfind(filesBulk.L1BFiles(i_fileL1B_input).name,'_SR_2')
                    filename_L1B=char([filesBulk.inputPath filesBulk.L1BFiles(i_fileL1B_input).name filesep 'enhanced_measurement.nc']);
                    cnf_p.input_type =2;
                end
            case 'DeDop'
                filename_L1B=char([filesBulk.inputPath filesBulk.L1BFiles(i_fileL1B_input).name]);
                cnf_p.input_type =1;
        end
    otherwise
        filename_L1B=char([filesBulk.inputPath filesBulk.L1BFiles(i_fileL1B_input).name]);
        cnf_p.input_type =1;
end

if targz_option_active_L1B
    %untar the file    
    extracted_files=untar(filename_L1B,filesBulk.inputPath);
    if any(~cellfun(@isempty,strfind(extracted_files,'.DBL')))
        %DBL like files
        filename_L1B=strrep(filename_L1B,'TGZ','DBL');
    elseif any(~cellfun(@isempty,strfind(extracted_files,'.nc')))
        filename_L1B=strrep(filename_L1B,'TGZ','.nc');
    elseif any(~cellfun(@isempty,strfind(extracted_files,'.mat')))
        filename_L1B=strrep(filename_L1B,'TGZ','.mat');
    else
        error('Extracted untar files do not contain any valid L1B input file (.DBL or .nc)');
    end    
end

disp(strcat('Modification Date',filesBulk.L1BFiles(i_fileL1B_input).date));

switch cnf_p.mission
    case {'S3A','S3B','S3'}
        switch cnf_p.L1proc
            case 'ESA'
                %a measurement netcdf file is included within the 
                filename_L1B_nopath=char(filesBulk.L1BFiles(i_fileL1B_input).name);
                fileext_L1B='.nc';
            case 'DeDop'
                [~,filename_L1B_nopath,fileext_L1B]=fileparts(filename_L1B);
        end
    otherwise
        [~,filename_L1B_nopath,fileext_L1B]=fileparts(filename_L1B);
end
%-------------------- File unqiue identifier extraction -------------------
%identifier of a given file is the initial and final acq times
switch cnf_p.mission
    case 'CS2'
        %--------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                switch fileext_L1B
                    case '.DBL'
                        file_unique_id=filename_L1B_nopath(20:50);
                    case '.nc'
                        %TBD
                end
            case 'ISD'
                %following a la Sentinel-3 format as per SEOMs projects
                file_unique_id=filename_L1B_nopath(17:47);
            case 'GPOD'
                file_unique_id=filename_L1B_nopath(17+7:47+7);
        end
    case 'S3'
        %--------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                 file_unique_id=filename_L1B_nopath(17:47);
            case 'DeDop'
                file_unique_id=filename_L1B_nopath;
        end

    case 'S6'
        %--------------------- Sentinel-6 -----------------------------------
        switch cnf_p.L1proc
            case 'ISD'
                switch fileext_L1B
                    case {'.nc','.NC'}
                        % according to name convention in IODD "JC-ID-ESA-GP-0175" issue 1.3
                        file_unique_id=filename_L1B_nopath;%filename_L1B_nopath(17:17+35);
                        %for the .mat file we assume the L1BS info is already
                        %saved in the .mat workspace
                    case '.mat'
                        file_unique_id = ''; %a way to overpass the extraction of 
                                             %the L1BS files as for
                                             %Sentinel-6 .mat we wil have
                                             %stack in the .mat workspace
                end
        end
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end

%--------------------------------------------------------------------------
%------------------------- L1BS file --------------------------------------
%--------------------------------------------------------------------------
if isfield(filesBulk,'L1BSFiles') && (~isempty(file_unique_id))
    %if it is not empty means we have provided some specific folder
    % we should look for the specific file linked to the L1B to be
    % processed
    %------------ Check whether file exist ----------------------
    %assume there is a single file definition:    
    idx_int=find(~cellfun(@isempty,strfind({filesBulk.L1BSFiles(:).name},file_unique_id)), 1);
    if ~isempty(idx_int)
        switch cnf_p.mission
            case {'S3A','S3B','S3'}
                switch cnf_p.L1proc
                    case 'ESA'
                        %a measurement netcdf file is included within the
                        filename_L1BS=char([filesBulk.inputPath_L1BS filesBulk.L1BSFiles(idx_int).name filesep 'measurement_l1bs.nc']);                        
                end
            otherwise
                filename_L1BS=char([filesBulk.inputPath filesBulk.L1BSFiles(idx_int).name]);
        end
        
        if targz_option_active_L1BS
            %untar the file
            extracted_files=untar(filename_L1BS,filesBulk.inputPath);
            if any(~cellfun(@isempty,strfind(extracted_files,'.DBL')))
                %DBL like files
                filename_L1BS=strrep(filename_L1BS,'TGZ','DBL');
            elseif any(~cellfun(@isempty,strfind(extracted_files,'.nc')))
                filename_L1BS=strrep(filename_L1BS,'TGZ','.nc');
            elseif any(~cellfun(@isempty,strfind(extracted_files,'.mat')))
                filename_L1BS=strrep(filename_L1BS,'TGZ','.mat');
            else
                error('Extracted untar files do not contain any valid L1B input file (.DBL or .nc)');
            end
        end
        
    else
        filename_L1BS=''; %define the file as empty since there is no such L1BS related to the L1B
    end
else
    filename_L1BS='';
end

%--------------------------------------------------------------------------
%------------------------- L2 file --------------------------------------
%--------------------------------------------------------------------------
if isfield(filesBulk,'L2Files') && (~isempty(file_unique_id))
    %if it is not empty means we have provided some specific folder
    % we should look for the specific file linked to the L1B to be
    % processed
    %------------ Check whether file exist ----------------------
    %assume there is a single file definition:
    
    %need to check the associated reference file of L1B
    L1b_date=datenum([file_unique_id(1:8) file_unique_id(10:15)],'yyyymmddHHMMSS');
    switch cnf_p.mission
        case {'S3A','S3B','S3'}
            switch cnf_p.L1proc
                case {'ESA','DeDop'}
                    %a measurement netcdf file is included within the
                    for i_file_L2=1:filesBulk.nFilesL2
                        filename_L2=char([filesBulk.inputPath_L2 filesBulk.L2Files(i_file_L2).name filesep 'enhanced_measurement.nc']);
                        filename_L1B_ref_L2_all=ncreadatt(filename_L2,'/','xref_altimeter_level1');
                        filename_L1B_ref_L2_cell=strsplit(filename_L1B_ref_L2_all,',');
                        for i_fcell=1:length(filename_L1B_ref_L2_cell)
                            L2_date= datenum([filename_L1B_ref_L2_cell{i_fcell}(17:17+7) filename_L1B_ref_L2_cell{i_fcell}(17+9:17+14)],'yyyymmddHHMMSS');
                            L1b_date=datenum([file_unique_id(1:8) file_unique_id(10:15)],'yyyymmddHHMMSS');
                            if abs(L2_date - L1b_date) < 60*10/86400 % less than 10 minutes between starting dates L1b & L2
                                disp('L2 found');
                                break;
                            else
                                if i_fcell==length(filename_L1B_ref_L2_cell)
                                    filename_L2='';
                                end
                            end
                        end
                        if ~isempty(filename_L2)
                            break;
                        end
                    end
            end
        otherwise
            idx_int=find(~cellfun(@isempty,strfind({filesBulk.L2Files(:).name},file_unique_id(1:15))), 1); %L2 products not have the same naming of acq.
            filename_L2=char([filesBulk.inputPath_L2 filesBulk.L2Files(idx_int).name]);
    end    
    
%     %old version
%     idx_int=find(~cellfun(@isempty,strfind({filesBulk.L2Files(:).name},file_unique_id(1:15))), 1); %L2 products not have the same naming of acq.
%     if ~isempty(idx_int)
%         switch cnf_p.mission
%             case {'S3A','S3B','S3'}
%                 switch cnf_p.L1proc
%                     case 'ESA'
%                         %a measurement netcdf file is included within the
%                         filename_L2=char([filesBulk.inputPath_L2 filesBulk.L2Files(idx_int).name filesep 'enhanced_measurement.nc']);                        
%                 end
%             otherwise
%                 filename_L2=char([filesBulk.inputPath_L2 filesBulk.L2Files(idx_int).name]);
%         end                       
%     else
%         filename_L2=''; %define the file as empty since there is no such L1BS related to the L1B
%     end
else
    filename_L2='';
end

end

