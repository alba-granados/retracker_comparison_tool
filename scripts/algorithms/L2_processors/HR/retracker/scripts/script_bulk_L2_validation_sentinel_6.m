num_pools=1;

set_default_plot_ASR;
figure_format='jpg';
res_fig='-r600';

win_size_detrending=20;
step_SWH=0.10;
sh_name_nc='ssh';
define_min_max_SWH=1;
min_SWH=0.5;
max_SWH=10;

flag_outliers_removal=0;
type_outliers_removal='hampel';
smooth_param=1;
plot_downsampling=25;
generate_plots=1;
generate_plot_SSH=1;
generate_plot_SWH=1;
generate_plot_sigma0=1;
generate_plot_nb=1;
generate_plot_COR=1;
generate_plot_misspointing=0;
generate_kml=0;
filter_land_surf_type=0;
compute_comparison_performance_tracks=0;
filter_ISR_baselines_mask=0;

% global filter_ISR_baselines_mask
% filter_ISR_baselines_mask=1;

%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
%NEED TO CHANGE IT TO YOUR LOCAL PATH
code_folder_full_path='C:/Users/eduard.makhoul/isardSAT/projects/SCOOP/processing/L2_GPP/Matlab_Code/';%'/home/emak/feina/projects/SCOOP/processing/L2_GPP/';%

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



%% ------------------------------ Simulation QR2 -----------------
mission='S6'

%---------------------- ESA reference path for L2 data --------------------
input_path_L2_ESA='' 
input_path_L1_ESA=''

%---------------- Starlab -------------------------------------------------
input_path_L2_STL=''

%------------------ GPOD --------------------------------------------------
input_path_L2_GPOD=''


%--------------- ISR baselines L2 -----------------------------------------
common_path='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/QR2/results/L2/'

name_b1='L1B GPP - GR MAT'
input_path_L2_ISR_b1='MAT/INPUT_CPP/02_S6A_OS2_2mSWH/S6A_OS21_0003_lat-40/No_windowing/RAW/data/'

name_b2='L1B GPP - GR Cpp'
input_path_L2_ISR_b2='CPP/02_S6A_OS2_2mSWH/S6A_OS21_0003_lat-40/No_windowing/RAW/data/'

processor_ID = {'MAT','CPP'};

% name_b3='RAW-threshold'
% input_path_L2_ISR_b3='RAW/Z_ML/remove_ambiguities/peak_retracker_wd_ref_256/data/'
% 
% name_b4='RMC-threshold'
% input_path_L2_ISR_b4='RMC/Z_ML/remove_ambiguities/peak_retracker_wd_ref_256/data/'


% input_path_L2_ISR_bs=strcat(common_path,{input_path_L2_ISR_b1,input_path_L2_ISR_b2,...
%                             input_path_L2_ISR_b3,input_path_L2_ISR_b4});
% name_bs={name_b1,name_b2,name_b3,name_b4};

input_path_L2_ISR_bs=strcat(common_path,{input_path_L2_ISR_b1,input_path_L2_ISR_b2});
name_bs={name_b1,name_b2};


%----------------- ISRs L1Bs ----------------------------------------------
input_path_L1_ISR_bs={''};
%input_path_L1_ISR_bs=strcat(common_path,{'L1B/S3/data/subset/'});


%------------------ Path to save results -------------------------------
path_comparison_results=strcat(common_path,'comparison_performance/L1B_GPP_L2_MAT_GPP/02_S6A_OS2_2mSWH/S6A_OS21_0003_lat-40/');

filename_mask_KML='';

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
                                'generate_plot_misspointing',generate_plot_misspointing,...
                                'filter_ISR_baselines_mask',filter_ISR_baselines_mask,...
                                'mission',mission,...
                                'processor_ID',processor_ID);                            

% %% ------------------------------ Simulation lat40 PREQR2 -----------------
% mission='S6'
% 
% %---------------------- ESA reference path for L2 data --------------------
% input_path_L2_ESA='' 
% input_path_L1_ESA=''
% 
% %---------------- Starlab -------------------------------------------------
% input_path_L2_STL=''
% 
% %------------------ GPOD --------------------------------------------------
% input_path_L2_GPOD=''
% 
% 
% %--------------- ISR baselines L2 -----------------------------------------
% common_path='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/PreQR2/results/L2/lat40/'
% %------------------ Path to save results -------------------------------
% path_comparison_results=strcat(common_path,'/comparison_performance_lat40_RAW_RMC/');
% % 
% % name_b1='MAT PTR boxcar'
% % input_path_L2_ISR_b1='MAT/RMC/Z_ML/PTR_boxcar/data/'
% % 
% % name_b2='C++ PTR boxcar'
% % input_path_L2_ISR_b2='CPP/RMC/Z_ML/PTR_boxcar/data/'
% % 
% % name_b3='MAT PTR Hamming'
% % input_path_L2_ISR_b3='MAT/RMC/Z_ML/PTR_0dot54/data/'
% % 
% % name_b4='C++ PTR Hamming'
% % input_path_L2_ISR_b4='CPP/RMC/Z_ML/PTR_0dot54/data/'
% 
% % name_b1='MAT zeros ML'
% % input_path_L2_ISR_b1='MAT/RAW/Z_ML/data/'
% % 
% % name_b2='C++ zeros ML'
% % input_path_L2_ISR_b2='CPP/RAW/Z_ML/data/'
% % 
% % name_b3='MAT no zeros ML'
% % input_path_L2_ISR_b3='MAT/RAW/NZ_ML/PTR_boxcar/data/'
% % 
% % name_b4='C++ no zeros ML'
% % input_path_L2_ISR_b4='CPP/RAW/NZ_ML/PTR_boxcar/data/'
% 
% 
% name_b1='MAT RAW'
% input_path_L2_ISR_b1='MAT/RAW/Z_ML/data/'
% 
% name_b2='C++ RAW'
% input_path_L2_ISR_b2='CPP/RAW/Z_ML/data/'
% 
% name_b3='MAT RMC'
% input_path_L2_ISR_b3='MAT/RMC/Z_ML/PTR_boxcar/data/'
% 
% name_b4='C++ RMC'
% input_path_L2_ISR_b4='CPP/RMC/Z_ML/PTR_boxcar/data/'
% 
% 
% input_path_L2_ISR_bs=strcat(common_path,{input_path_L2_ISR_b1,input_path_L2_ISR_b2,...
%                         input_path_L2_ISR_b3,input_path_L2_ISR_b4});
% name_bs={name_b1,name_b2,name_b3,name_b4};
% 
% 
% %----------------- ISRs L1Bs ----------------------------------------------
% input_path_L1_ISR_bs={''};
% %input_path_L1_ISR_bs=strcat(common_path,{'L1B/S3/data/subset/'});
% 
% 
% filename_mask_KML='';
% 
% performance_global_ISR_baselines_paralelization_new_GPOD(input_path_L2_ISR_bs,name_bs,path_comparison_results,...
%                                 'figure_format',figure_format,...
%                                 'res_fig',res_fig,...
%                                 'win_size_detrending',win_size_detrending,...
%                                 'step_SWH',step_SWH,...
%                                 'sh_name_nc',sh_name_nc,...
%                                 'flag_outliers_removal',flag_outliers_removal,...
%                                 'type_outliers_removal',type_outliers_removal,...
%                                 'smooth_param',smooth_param,...
%                                 'input_path_L2_ESA',input_path_L2_ESA,...
%                                 'input_path_L1_ISR_bs',input_path_L1_ISR_bs,...
%                                 'input_path_L1_ESA',input_path_L1_ESA,...
%                                 'input_path_L2_STL',input_path_L2_STL,...
%                                 'input_path_L2_GPOD',input_path_L2_GPOD,...
%                                 'filename_mask_KML',filename_mask_KML,...
% 								'num_pools',num_pools,...
%                                 'generate_plots',generate_plots,...
%                                 'plot_downsampling',plot_downsampling,...
%                                 'define_min_max_SWH',define_min_max_SWH,...
%                                 'min_SWH',min_SWH,...
%                                 'max_SWH',max_SWH,...
%                                 'generate_kml',generate_kml,...
%                                 'filter_land_surf_type',filter_land_surf_type,...
%                                 'compute_comparison_performance_tracks',compute_comparison_performance_tracks,...
%                                 'generate_plot_SSH',generate_plot_SSH,'generate_plot_SWH',generate_plot_SWH,...
%                                 'generate_plot_sigma0',generate_plot_sigma0,'generate_plot_nb',generate_plot_nb,...
%                                 'generate_plot_COR',generate_plot_COR,...
%                                 'generate_plot_misspointing',generate_plot_misspointing,...
%                                 'filter_ISR_baselines_mask',filter_ISR_baselines_mask,...
%                                 'mission',mission);       

