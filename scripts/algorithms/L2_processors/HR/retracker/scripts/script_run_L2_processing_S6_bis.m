%profile on;
%% -------------------- Number of Pools -----------------------------------
num_pools=1;
num_threads=4;
%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_GPP/Matlab_Code/';
%code_folder_full_path='/home/emak/Sentinel-6/code/L2_GPP/Matlab_code/';

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



%% Run different parallel processings
common_path='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\FFS\Euribia_v2dot1\S6A_OS10_0001_1400RC_2017\results\';
%common_path='/home/dl380/test_data/Sentinel-6/data/AR/results/';


                                 
input_path_L1B = strcat(common_path,{'L1B\FF\RAW\','L1B\FF\RMC\','L1B\HR_DDP\RAW\','L1B\HR_DDP\RMC\'});                                 
%input_path_L1B = strcat(common_path,{'L1B\HR_DDP\RAW\'});                                 

%output_path_L2={'C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/QR2/results/L2/MAT/INPUT_MAT/02_S6A_OS2_2mSWH/S6A_OS21_0003_lat-40/LR_RMC_from_RAW_adpt_code_AR/'}
output_path_L2=strcat(strrep(input_path_L1B,'L1B','L2'));
%
                                                          
                               
cnf_chd_cst_path=cell(1,length(input_path_L1B));
cnf_chd_cst_path(:)={'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\FFS\cnf_chd_cst_LUTs_L2\'};
%cnf_chd_cst_path(:)={'/home/emak/Sentinel-6/code/L2_GPP/L2_cnf_cst_chd_LUTs/'};
filter_input_file = cell(1,length(input_path_L1B));
filter_input_file(:) = {'.mat'};
filename_mask_KML = cell(1,length(input_path_L1B));
filename_mask_KML(:) = {''}; 
num_pools = num_pools*ones(1,length(input_path_L1B));

proc_bsln_id = cell(1,length(input_path_L1B));
proc_bsln_id(:)={'FF_S6_RAW.json','FF_S6_RMC.json','HRDDP_S6_RAW.json','HRDDP_S6_RMC.json'};
%proc_bsln_id(:)={'HRDDP_S6_RAW.json'};

PTR_along_external = [];
PTR_across_external = [];

version_matlab=version;

if num_threads > 1
	%open pools
	if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_threads);
    else
        matlabpool('open',num_threads);
    end    
	
    parfor i_proc=1:length(input_path_L1B)
        call_L2_processing_parallel(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            proc_bsln_id,filter_input_file,num_pools,...
            filename_mask_KML,...
			PTR_along_external,PTR_across_external,i_proc);
    end
    
	%close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end
else
    for i_proc=1:length(input_path_L1B)
        call_L2_processing_parallel(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            proc_bsln_id,filter_input_file,num_pools,...
            filename_mask_KML,...
			PTR_along_external,PTR_across_external,i_proc);
    end
end
