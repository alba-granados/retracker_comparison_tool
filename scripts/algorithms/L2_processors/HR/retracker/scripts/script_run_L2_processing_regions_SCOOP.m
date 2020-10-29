ttotal=tic;
%profile on;
%% -------------------- Number of Pools -----------------------------------
num_pools=1;

%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
%code_folder_full_path='/home/emak/SCOOP/L2_GPP/code_delivery/test_24_01_2019/';
code_folder_full_path='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\L2_GPP\Matlab_Code\';


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

%% --------------------------------------- Agulhas ----------------------------------------------------
%NEED TO CHANGE TO YOUR LOCAL PATHS
proc_bsln_id = 'L1B_ISR_CS2.json'
input_path_L1B='C:\Users\eduard.makhoul\Desktop\test_L2\input_L1B\';
output_path_L2='C:\Users\eduard.makhoul\Desktop\test_L2\L2\resampling_index\';
cnf_chd_cst_path='C:\Users\eduard.makhoul\Desktop\test_L2\cnf_cst_chd_LUTs\';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
filename_dist_to_coast='C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\Dist2Coast\dist2coast_ferret_listing.nc'
filename_mask_KML='';
attitude_extraction_GPOD = 0;
input_path_L2_GPOD = '';
Sig0AtmCorr_path = 'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\Sig0AtmCorr\';
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,...
            'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            'input_path_L2_GPOD',input_path_L2_GPOD,...
            'Sig0AtmCorr_path',Sig0AtmCorr_path,...
            'proc_bsln_id',proc_bsln_id,...
            'filename_dist_to_coast',filename_dist_to_coast)




% %% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
% %NEED TO CHANGE IT TO YOUR LOCAL PATH
% %code_folder_full_path='/home/emak/SCOOP/L2_GPP/Matlab_Code/';
% code_folder_full_path='/home/emak/SCOOP/prova/Matlab_Code/';
% 
% 
% cd(code_folder_full_path);
% FolderInfo=dir(code_folder_full_path);
% FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
% aux=struct2cell(FolderInfo); %into a cell array where first row is 
% folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
% clear aux;
% folders=strcat(code_folder_full_path,folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])&~strcmp(folders,['old'])));
% for i_folder=1:length(folders)
%     addpath(genpath(char(folders(i_folder))));
% end
% 
% 
% 
% 
% 
% %% --------------------------------------- prova ----------------------------------------------------
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final_coast/L1B/agulhas/2012/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final_coast/L2/agulhas/2012/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='/home/emak/SCOOP/prova/aux_files/agulhas.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
%             'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
%             'attitude_extraction_GPOD',attitude_extraction_GPOD,...
%             'input_path_L2_GPOD',input_path_L2_GPOD,...
%             'Sig0AtmCorr_path',Sig0AtmCorr_path)
% 			
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final_coast/L1B/agulhas/2013/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final_coast/L2/agulhas/2013/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='/home/emak/SCOOP/prova/aux_files/agulhas.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
%             'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
%             'attitude_extraction_GPOD',attitude_extraction_GPOD,...
%             'input_path_L2_GPOD',input_path_L2_GPOD,...
%             'Sig0AtmCorr_path',Sig0AtmCorr_path)		



			
			
time=toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time (all regions): ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
%exit

% %% --------------------------------------- east_pacific ----------------------------------------------------
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final/L1B/east_pacific/2012/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final/L2/east_pacific/2012/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';%'/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            % 'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            % 'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            % 'input_path_L2_GPOD',input_path_L2_GPOD,...
            % 'Sig0AtmCorr_path',Sig0AtmCorr_path)

% %% --------------------------------------- west_pacific ----------------------------------------------------
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final/L1B/west_pacific/2013/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final/L2/west_pacific/2013/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';%'/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            % 'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            % 'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            % 'input_path_L2_GPOD',input_path_L2_GPOD,...
            % 'Sig0AtmCorr_path',Sig0AtmCorr_path)
			
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final/L1B/west_pacific/2012/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final/L2/west_pacific/2012/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';%'/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            % 'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            % 'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            % 'input_path_L2_GPOD',input_path_L2_GPOD,...
            % 'Sig0AtmCorr_path',Sig0AtmCorr_path)
			
% %% --------------------------------------- central_pacific ----------------------------------------------------
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final/L1B/central_pacific/2013/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final/L2/central_pacific/2013/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';%'/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            % 'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            % 'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            % 'input_path_L2_GPOD',input_path_L2_GPOD,...
            % 'Sig0AtmCorr_path',Sig0AtmCorr_path)
			
% proc_bsln_id='S3_zp2_Hamm.json';
% %NEED TO CHANGE TO YOUR LOCAL PATHS
% input_path_L1B='/data/emak/SCOOP/results/Phase_2_final/L1B/central_pacific/2012/data/';
% output_path_L2='/data/emak/SCOOP/results/Phase_2_final/L2/central_pacific/2012/';
% cnf_chd_cst_path='/home/emak/SCOOP/prova/aux_files/';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/';
% filename_mask_KML='';%'/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';
% attitude_extraction_GPOD = 0;
% input_path_L2_GPOD = '';%'/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';
% Sig0AtmCorr_path = '';
% L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            % 'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id,...
            % 'attitude_extraction_GPOD',attitude_extraction_GPOD,...
            % 'input_path_L2_GPOD',input_path_L2_GPOD,...
            % 'Sig0AtmCorr_path',Sig0AtmCorr_path)	