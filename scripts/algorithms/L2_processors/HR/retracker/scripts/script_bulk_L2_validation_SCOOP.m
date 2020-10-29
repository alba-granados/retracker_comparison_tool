%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
%code_folder_full_path='/home/emak/SCOOP/L2_GPP/Matlab_code/';%
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

% ---------------- options  -----------------------------------------------
num_pools=1;
set_default_plot;
figure_format='png';
res_fig='-r150';
win_size_detrending=20;
step_SWH=0.10;
sh_name_nc='ssh';
define_min_max_SWH=1;
min_SWH=0.5;
max_SWH=10;

flag_outliers_removal=0;
type_outliers_removal='hampel';
smooth_param=1;
plot_downsampling=1;
generate_plots=1;
generate_plot_SSH=1;
generate_plot_SWH=1;
generate_plot_sigma0=1;
generate_plot_nb=0;
generate_plot_COR=1;
generate_plot_Misfit=1;
generate_plot_misspointing=0;
generate_kml=0;
filter_land_surf_type=0;
compute_comparison_performance_tracks=0;
filter_ISR_baselines_mask=0;
remove_GPOD_misfit_thresh = 0;
GPOD_misfit_thresh = 3.5;
remove_ISR_misfit_thresh = 0;
ISR_misfit_thresh = 3.5;


%% ------------------------------ Agulhas ------------------------------------------------------------------------------
%int_path ='/home/dl380/test_data/SCOOP/';
int_path ='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/Phase_2/';

%---------------------- ESA reference path for L2 data --------------------
input_path_L2_ESA='';%
input_path_L1_ESA='';%strcat(int_path,'input/ESA_IPF_CS2/L1B/agulhas/2013/');%

%---------------- Starlab -------------------------------------------------
input_path_L2_STL='';%strcat(int_path,'/data/DATA/SCOOP/STL/L2_S3/agulhas/2013_reprocess_phase1_CP_new/');

%------------------ GPOD --------------------------------------------------
input_path_L2_GPOD='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/inputs/L2_GPOD_SARVATORE/Agulhas/2013/data/';
%input_path_L2_GPOD='/home/dl380/test_data/SCOOP/input/GPOD/agulhas/2013_LUT_PTR/';%strcat(int_path,'input/GPOD/agulhas/2013/');

%--------------- ISR baselines L2 -----------------------------------------
%common_path=strcat(int_path,'results/Phase_2/agulhas/');
common_path=strcat(int_path,'Phase_2_final_nocut/');


input_path_L2_ISR_bs=strcat(common_path,...
        {'L2/agulhas/2013'},'/data/subset/');

name_bs={'ISR'};


%----------------- ISRs L1Bs ----------------------------------------------
input_path_L1_ISR_bs=strcat(common_path,'L1B/agulhas/2013',...
                            {'/data/'});

L1B_processor = {'ISR'};

%------------------ KML mask ----------------------------------------------
%filename_mask_KML='/home/emak/SCOOP/L2_GPP/conf_cst_chd/Agulhas_mask_close_coast.kml';%strcat(int_path,'inputs/configuration_L2_GPOD/Agulhas_mask_close_coast.kml');
filename_mask_KML='';%'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\configuration_L2\Agulhas_mask_close_coast.kml';

%------------------ Path to save results ----------------------------------
path_comparison_results=strcat(common_path,'comparison_performance_comp_load/');
%path_comparison_results=strcat(common_path,'comparison_performance/');



performance_global_ISR_baselines_paralelization_new_GPOD(input_path_L2_ISR_bs,name_bs,path_comparison_results,...
                                'figure_format',figure_format,...
                                'res_fig',res_fig,...
                                'win_size_detrending',win_size_detrending,...
                                'step_SWH',step_SWH,...
                                'sh_name_nc',sh_name_nc,...
                                'flag_outliers_removal',flag_outliers_removal,...
                                'type_outliers_removal',type_outliers_removal,...
                                'smooth_param',smooth_param,...
                                'input_path_L2_ESA',input_path_L2_ESA,...
                                'input_path_L1_ISR_bs',input_path_L1_ISR_bs,...
                                'input_path_L1_ESA',input_path_L1_ESA,...
                                'input_path_L2_STL',input_path_L2_STL,...
                                'input_path_L2_GPOD',input_path_L2_GPOD,...
                                'filename_mask_KML',filename_mask_KML,...
								'num_pools',num_pools,...
                                'generate_plots',generate_plots,...
                                'plot_downsampling',plot_downsampling,...
                                'define_min_max_SWH',define_min_max_SWH,...
                                'min_SWH',min_SWH,...
                                'max_SWH',max_SWH,...
                                'generate_kml',generate_kml,...
                                'filter_land_surf_type',filter_land_surf_type,...
                                'compute_comparison_performance_tracks',compute_comparison_performance_tracks,...
                                'generate_plot_SSH',generate_plot_SSH,'generate_plot_SWH',generate_plot_SWH,...
                                'generate_plot_sigma0',generate_plot_sigma0,'generate_plot_nb',generate_plot_nb,...
                                'generate_plot_COR',generate_plot_COR,...
                                'generate_plot_Misfit',generate_plot_Misfit,...
                                'generate_plot_misspointing',generate_plot_misspointing,...
                                'filter_ISR_baselines_mask',filter_ISR_baselines_mask,...
                                'L1B_processor',L1B_processor,...
                                'remove_GPOD_misfit_thresh',remove_GPOD_misfit_thresh,...
                                'GPOD_misfit_thresh',GPOD_misfit_thresh,...
                                'remove_ISR_misfit_thresh',remove_ISR_misfit_thresh,...
                                'ISR_misfit_thresh',ISR_misfit_thresh);

% 

