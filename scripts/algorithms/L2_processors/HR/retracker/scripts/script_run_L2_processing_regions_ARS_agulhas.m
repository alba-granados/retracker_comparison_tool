%profile on;
%% -------------------- Number of Pools -----------------------------------
num_pools=1;
%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_GPP/Matlab_Code/';

cd(code_folder_full_path);
FolderInfo=dir(code_folder_full_path);
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat(code_folder_full_path,folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])&~strcmp(folders,['old'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end



%% Agulhas
proc_bsln_id='S3.json';
%NEED TO CHANGE TO YOUR LOCAL PATHS
input_path_L1B='C:\Users\eduard.makhoul\Desktop\temporal\data\agulhas\L1B\s3\';
output_path_L2='C:\Users\eduard.makhoul\Desktop\temporal\data\agulhas\L2\s3\';
cnf_chd_cst_path='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\configuration_L2\';
filename_mask_KML='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\configuration_L2\Agulhas_mask_close_coast.kml';
attitude_extraction_GPOD = 1;
input_path_L2_GPOD = 'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\L2_GPOD_SARVATORE\Agulhas\2013\subset\';
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            'input_path_L2_GPOD',input_path_L2_GPOD)

