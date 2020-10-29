ttotal=tic;

%% -------------------- Number of Pools -----------------------------------
num_pools=1;
%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
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

common_path = '\\N8800PRO\share\ISARDS\emak\L2_processor\test_data\';%'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\L2_GPP\test_data\';

%% ----------------------- L1B S-3 from ESA -------------------------------
proc_bsln_id='L1B_ESA_S3.json';
input_path_L1B=strcat(common_path,'inputs\L1B_S3_ESA\');
input_path_L2 =strcat(common_path,'inputs\L2_S3_ESA\');
output_path_L2=strcat(common_path,'results\L1B_S3_ESA\');
cnf_chd_cst_path=strcat(common_path,'inputs\cnf_cst_chd_LUTs\');
filename_mask_KML=strcat(common_path,'\inputs\masks\mask_track_029_256.kml');
filename_dist_to_coast=strcat(common_path,'dist_to_coast\dist2coast_ferret_listing.nc');
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'input_path_L2',input_path_L2,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,...           
            'proc_bsln_id',proc_bsln_id,...
            'filename_dist_to_coast',filename_dist_to_coast);        
                    %'proc_bsln_id',proc_bsln_id,...


%% ----------------------- L1B CS-2 from isardSAT (SCOOP)------------------
proc_bsln_id='L1B_ISR_CS2.json';
input_path_L1B=strcat(common_path,'inputs\L1B_CS2_ISR\');
input_path_L2 ='';
output_path_L2=strcat(common_path,'results\L1B_CS2_ISR\');
cnf_chd_cst_path=strcat(common_path,'inputs\cnf_cst_chd_LUTs\');
filename_mask_KML='';
filename_dist_to_coast='';
Sig0AtmCorr_path = strcat(common_path,'inputs\Sig0AtmCorr\');
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'input_path_L2',input_path_L2,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,...           
            'proc_bsln_id',proc_bsln_id,...
            'filename_dist_to_coast',filename_dist_to_coast,...
            'Sig0AtmCorr_path',Sig0AtmCorr_path); 
        

%% ----------------------- L1B CS-2 from isardSAT (SHAPE)------------------
% to run an example of code to be used for generating .exe and delivery it
% to ESA (so a reduced set of configuration options is input to the processor compared 
% to the case of using the generic cnf of the processor: so ESA can not use other options to process)
% TO do so the cnf file with reduced options shall be input to the processor and the user shall set up 
% reduced_set_cnf.SCOOP=0; reduced_set_cnf.SHAPE=1; in the
% L2_bulk_processing_paralelization.m function

proc_bsln_id='L1B_ISR_CS2_SHAPE.json';
input_path_L1B=strcat(common_path,'inputs\L1B_CS2_ISR_SHAPE\');
input_path_L2 ='';
output_path_L2=strcat(common_path,'results\L1B_CS2_ISR_SHAPE\');
cnf_chd_cst_path=strcat(common_path,'inputs\cnf_cst_chd_LUTs\');
filename_mask_KML=strcat(common_path,'inputs\masks\Region_polygon_Amazon_river.kml');
filename_dist_to_coast='';
Sig0AtmCorr_path = strcat(common_path,'inputs\Sig0AtmCorr\');
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'input_path_L2',input_path_L2,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,...           
            'proc_bsln_id',proc_bsln_id,...
            'filename_dist_to_coast',filename_dist_to_coast,...
            'Sig0AtmCorr_path',Sig0AtmCorr_path); 

        
%% ----------------------- L1B S-6 from isardSAT GPP in C++ ---------------
proc_bsln_id='L1B_ISR_S6.json';
input_path_L1B=strcat(common_path,'inputs\L1B_S6_CPP\');
input_path_L2 ='';
output_path_L2=strcat(common_path,'results\L1B_S6_CPP\');
cnf_chd_cst_path=strcat(common_path,'inputs\cnf_cst_chd_LUTs\');
filter_input_file = '.nc'; %filter
filename_mask_KML='';
filename_dist_to_coast='';
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,...
            'input_path_L2',input_path_L2,...
            'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,...           
            'proc_bsln_id',proc_bsln_id,...
            'filename_dist_to_coast',filename_dist_to_coast,...
            'filter_input_file',filter_input_file);         
        
                    
time=toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time (all regions): ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
%exit
