%profile on;

%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
    %NEED TO CHANGE IT TO YOUR LOCAL PATH
    code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_GPP/Matlab_Code/';
    %code_folder_full_path='/home/emak/Sentinel-6/code/L2_GPP/Matlab_Code/';
    
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
    

%% Input paths
common_path='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/';
%common_path='/home/dl380/test_data/Sentinel-6/data/AR/results/';
    
% input data
input_path_L1B_total = strcat(common_path,{'L1B_GPP/hr_raw/'});
%                                     
% {'L1B/S6A_OS10_0001_1400RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS16_0001_1500RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS20_0001_1400RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS21_0001_2500RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS22_0001_1400RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS23_0001_1400RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS24_0001_1400RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS30_0001_1500RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS40_0001_1500RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS41_0001_1500RC_2017/LR_RMC/',...
%                                         'L1B/S6A_OS10_0001_1400RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS16_0001_1500RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS20_0001_1400RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS21_0001_2500RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS22_0001_1400RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS23_0001_1400RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS24_0001_1400RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS30_0001_1500RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS40_0001_1500RC_2017/LR_RMC_FROM_RAW/',...
%                                         'L1B/S6A_OS41_0001_1500RC_2017/LR_RMC_FROM_RAW/'}                                    
                                    
                                    

for i_loop=1:length(input_path_L1B_total)
    
    %% -------------------- Number of Pools -----------------------------------
    num_pools=1;
    num_threads=1;
    
    
    %% External PTRs
    [PTR_along_external,PTR_across_external]=meshgrid([0.5,1.0],[0.5,1.0]);%meshgrid(0.5:0.125/4:1.0,0.5:0.125/4:1.0);
    PTR_along_external = reshape(PTR_along_external,[1 length(PTR_along_external(:))]);
    PTR_across_external = reshape(PTR_across_external,[1 length(PTR_across_external(:))]);

    
    %% Input paths    
    input_path_L1B=cell(1,length(PTR_along_external));
    input_path_L1B(:)=input_path_L1B_total(i_loop);
    
    % output path
    output_path_L2=strcat(strrep(input_path_L1B,'L1B_GPP','L2_from_L1B_GPP'));
        
    cnf_chd_cst_path=cell(1,length(input_path_L1B));
    cnf_chd_cst_path(:)={'C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/configuration_L2/'};
    %cnf_chd_cst_path(:)={'/home/emak/Sentinel-6/code/L2_GPP/L2_cnf_cst_chd_LUTs/'};
    filter_input_file = cell(1,length(input_path_L1B));
    filter_input_file(:) = {'.nc'};
    filename_mask_KML = cell(1,length(input_path_L1B));
    filename_mask_KML(:) = {''};
    num_pools = num_pools*ones(1,length(input_path_L1B));
    
    proc_bsln_id = cell(1,length(input_path_L1B));
    proc_bsln_id(:)={'S6_zp2.json'};
    
    version_matlab=version;
    
    %% Run different parallel processings
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
end
%exit