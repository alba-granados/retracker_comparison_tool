% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This is the main file to run the model fit tool WP1700 within the S6 GPP
% P4 L1 project.
% It runs over the different input L1B and L2 files defined in a given folder, using a given
% set of configuration and characterization files provided in a given
% folder, storing the results in a specified folder. 
% L2 product can be generated in this tool or they can be called in a 
% previously generated folder. Fitted curves and performance graphics are plotted. 
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ----- / isardSAT
% 
% Last revision:    Alba Granados / isardSAT V1 07/09/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%  - (if cnf_tool.run_L2=1) L2_bulk_processing_paralelization: isardSAT retracker
%  - read_L2: reads geophysical parameters stored in L2 file and equivalent L1B file
%  - model_fit_to_power_waveform:  plots the fitted curves from retrieved parameters in L2 files
%  - performance_baselines_S6: comparison of the geophsycial retrievals of different input
%               baselines of isardSAT
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% RESTRICTIONS:
% - It is assumed all the L1B input files correspond
% to the same mission and have been processed with the same configuration
% parameters
% - LR fitted wavefors cannot be plotted from L2 LR netcdf file in this
% version.
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% V1


%% ----------- path to model fit tool WP1700 source code  ---------------

path_to_waveformfittingtool = strcat(pwd, filesep, 'scripts');
fprintf('\nPath to waveform fitting tool (current_directory/scripts/):\n%s\n\n', path_to_waveformfittingtool);

dir_bar = filesep;

FolderInfo=dir(path_to_waveformfittingtool);
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); 
folders=(aux(1,[FolderInfo.isdir]))'; 
clear aux;
folders=strcat(path_to_waveformfittingtool,filesep, folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end

%% ----------- RUN TOOL OPTIONS ----------------------------------

% % wait for user input
% cnf_waveformfittingtool_path = input('Enter path to waveform fitting tool configuration file (leave empty and press enter if testing):\n', 's');
% if isempty(cnf_waveformfittingtool_path)
%     cnf_waveformfittingtool_path  = strcat(pwd, filesep);    
%     cnf_waveformfittingtool_filename  = strcat(cnf_waveformfitting_path, 'cnf_waveformfittingtool.json');    
%     fprintf('Default: %s\n', cnf_waveformfittingtool_path);
% end

% Default path to configuration file
cnf_waveformfittingtool_path  = strcat(pwd, filesep);    
tool_bsln_id = 'S6_LX2_CL_HR_dl20.json';
inputFiles      =   dir(cnf_waveformfittingtool_path);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
if ~isempty(tool_bsln_id)
    cnf_waveformfittingtool_filename=[cnf_waveformfittingtool_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_waveformfittingtool','_',tool_bsln_id])))).name];
else
    cnf_waveformfittingtool_filename=[cnf_waveformfittingtool_path inputFiles(~cellfun(@isempty,strfind(aux,char('cnf_waveformfittingtool')))).name];
end
fprintf('Path to waveform fitting tool configuration file:\n%s\n', cnf_waveformfittingtool_path);
% cnf_waveformfittingtool_filename  = strcat(cnf_waveformfittingtool_path, 'cnf_waveformfittingtool.json');    

[cnf_tool]=read_CNF_tooloptions(cnf_waveformfittingtool_filename);

mission=cnf_tool.mission; 
name_bs=cnf_tool.name_bs; %{'RAW','RMC','LRM'};
fprintf('\nBaselines to be processed:\n');
disp(name_bs);

num_baselines=length(name_bs);

fprintf('Run L2 processor to generate the L2 products for baseline? [1/0] \n');
for b=1:length(cnf_tool.run_L2)
    fprintf('\tBaseline %s: %d\n', char(name_bs(b)), cnf_tool.run_L2(b));
end


% ----------- OUTPUT PLOTS OPTIONS ----------------------------------
% fitted waveforms
figure_format=cnf_tool.figure_format;
res_fig=cnf_tool.res_fig;
color_bs=cnf_tool.color_bs;
plot_downsampling=cnf_tool.plot_downsampling;
text_interpreter =cnf_tool.text_interpreter;
textbox_fontsize =cnf_tool.textbox_fontsize;
legend_fontsize =cnf_tool.legend_fontsize;
default_fontsize =cnf_tool.default_fontsize;
overlay_baselines =cnf_tool.overlay_baselines; 
num_pools = cnf_tool.num_pools; 

% retrieved parameters performance analysis
LineStyle=cnf_tool.LineStyle;
marker_bs=cnf_tool.marker_bs;
win_size_detrending=cnf_tool.win_size_detrending;
flag_outliers_removal=cnf_tool.flag_outliers_removal;
smooth_param=cnf_tool.smooth_param;
generate_plots=cnf_tool.generate_plots;
generate_plot_SSH=cnf_tool.generate_plot_SSH;
generate_plot_SWH=cnf_tool.generate_plot_SWH;
generate_plot_sigma0=cnf_tool.generate_plot_sigma0;
generate_plot_COR=cnf_tool.generate_plot_COR;
filter_ISR_baselines_mask=cnf_tool.filter_ISR_baselines_mask; % not in use
type_outliers_removal=cnf_tool.type_outliers_removal; 


%  Set the default configuration for display
mida = get(0,'ScreenSize');
mida(3:4)=cnf_tool.default_figuresize;

set(0,'defaultFigurePaperUnits','points'); set(0,'defaultFigurePosition',mida);
set(groot,'defaultAxesFontSize', cnf_tool.default_fontsize); set(groot,'defaultTextFontSize', cnf_tool.default_fontsize);
set(groot,'defaultLegendFontSize', cnf_tool.legend_fontsize);
set(groot, 'defaultAxesTickLabelInterpreter',cnf_tool.text_interpreter); set(groot, 'defaultLegendInterpreter',cnf_tool.text_interpreter);
set(0,'defaultFigurePaperPosition', mida);
set(0,'defaultLineLineWidth',cnf_tool.default_linewidth);   % set the default line width
set(0,'DefaultFigureVisible',cnf_tool.visible_figures);
 
% for L2 HR processor, please check /retracker/plotting/set_default_plot.m
% ( or comment call in /retracker/algorithms/L2_processing.m )


%% ------------- DEFINITION OF THE INPUTS/OUTPUT PATHS --------------------
% % --------------- main DATA path S6 GPP model fit tool ---------------------------
output_data_path  = cnf_tool.output_data_path;    
fprintf('\nOutput data path:\n%s\n\n', output_data_path);

% %-------- path to inputs ISRs L1B containing the L1B file/s to be processed --------
fprintf('Path to L1B input data:\n');
input_path_L1_ISR_bs = cnf_tool.input_path_L1_ISR_bs;
for i_baseline=1:num_baselines
    fprintf('%s\n', char(input_path_L1_ISR_bs(i_baseline)));
end

% % ------- path to inputs/outputs ISRs L2 containing the L2 file/s to be processed ------
input_path_L2_ISR_bs = cnf_tool.input_path_L2_ISR_bs;
output_path_L2_ISR_bs_common = strcat(output_data_path, 'L2_processor', filesep);
fprintf('\nPath to L2 input data:\n');
output_path_L2_ISR_bs = cell(1,num_baselines); % output path of L2 product if cnf_tool.run_L2
for i_baseline=1:num_baselines
    output_path_L2_ISR_bs(i_baseline) = cellstr(strcat(output_path_L2_ISR_bs_common, name_bs{i_baseline}, filesep));                 
    if cnf_tool.run_L2(i_baseline) % L2 processsor has run for this baseline
        if ~exist(char(output_path_L2_ISR_bs(i_baseline)), 'dir')
            mkdir(char(output_path_L2_ISR_bs(i_baseline)))
        end
        input_path_L2_ISR_bs{i_baseline} = cellstr(strcat(output_path_L2_ISR_bs(i_baseline), 'data', filesep));   
    else % L2 processor not run
        output_path_L2_ISR_bs{i_baseline} =[];
    end
    fprintf('%s\n', char(input_path_L2_ISR_bs{i_baseline}));
end

fprintf('\nPath to L2 output data:\n');
for i_baseline=1:num_baselines
    if ~isempty(output_path_L2_ISR_bs{i_baseline})
        fprintf('(%s) %s\n', name_bs{i_baseline}, char(output_path_L2_ISR_bs{i_baseline}));
    end
end

%------------------ path to output plots of tool  -------------------------------
% model fit tool / fitted waveforms
path_fitted_waveforms =...
    strcat(output_data_path, 'model_fit_to_power_waveform', filesep, 'fitted_waveforms', filesep);
if ~exist(char(path_fitted_waveforms), 'dir')
    mkdir(char(path_fitted_waveforms));
end
% model fit tool / performace comparison
path_comparison_results=strcat(output_data_path, 'model_fit_to_power_waveform', filesep, 'performance', filesep);
mkdir(path_comparison_results);

fprintf('\nPath to output tool plots:\n%s\n%s\n', path_fitted_waveforms, path_comparison_results);

% % ----------- path to L2 HR processor source code ---------------
path_to_L2_processor  = cnf_tool.path_to_L2_processor; strcat(path_to_waveformfittingtool, filesep, 'L2_processors', filesep, 'HR', filesep, 'retracker', filesep);   
fprintf('\nPath to L2 HR processor:\n%s\n\n', path_to_L2_processor);

FolderInfo=dir(path_to_L2_processor);
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); 
folders=(aux(1,[FolderInfo.isdir]))'; 
clear aux;
folders=strcat(path_to_L2_processor,folders(~strcmp(folders,['.git'])&~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])&~strcmp(folders,['old'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end
proc_bsln_id='L1B_ISR_S6.json';

% % ------- path to config files for L2 processing --------------------------------
% % path containing the different configurationcnf_file.m 
% % (configuration processing parameters for L2), characterization chd_file.m
% % (different missions parameters) and the constant definition cst_file.m; 
% % the Look Up Tables (LUTs) for the f0 and f1 functions to be used by 
% % the single look waveform model are also included in this folder
% % cnf_chd_cst_path = [common_data_path, 'inputs_cnf_chd_cst_LUTs_Euribia', filesep];
% cnf_chd_cst_path = [input_data_path, 'inputs_cnf_chd_cst_LUTs_Euribia', filesep];


cnf_chd_cst_path  = cnf_tool.cnf_chd_cst_path;
fprintf('Path to cnf, chd, cst, LUTs files for L2 HR processor:\n%s\n\n', cnf_chd_cst_path);


% ----------------- KML files ---------------------------------------------
filename_mask_KML='';


%% run L2 processors and define L2 product input path for model fit tool

for i_baseline=1:num_baselines 
    if cnf_tool.run_L2(i_baseline)
        switch char(name_bs(i_baseline))
            case {'LRM', 'LROS-RAW', 'LROS-RMC'} 
                fprintf('\nRunning L2 LR processor...\n');

                if ~exist(char(output_path_L2_ISR_bs(i_baseline)), 'dir')
                    mkdir(char(output_path_L2_ISR_bs(i_baseline)));
                end

                L2_LRM_S6_modelfittool(char(input_path_L1_ISR_bs(i_baseline)), char(output_path_L2_ISR_bs(i_baseline)), cnf_chd_cst_path,  cnf_tool,...
                    'proc_bsln_id',proc_bsln_id, 'MODE', char(name_bs(i_baseline)));
                % redefine L2 product input path for model fit tool is output/data/ of L2 processors
                input_path_L2_ISR_bs{i_baseline} = cellstr(strcat(output_path_L2_ISR_bs(i_baseline), 'data', filesep));   

            otherwise
                fprintf('Running L2 HR processor...\n\n');

                if ~exist(char(output_path_L2_ISR_bs(i_baseline)), 'dir')
                    mkdir(char(output_path_L2_ISR_bs(i_baseline)));
                end

                L2_bulk_processing_paralelization(char(input_path_L1_ISR_bs(i_baseline)),char(output_path_L2_ISR_bs(i_baseline)),cnf_chd_cst_path, ...
                    'proc_bsln_id',proc_bsln_id, 'num_pools',num_pools);
                % redefine L2 product input path for model fit tool is output/data/ of L2 processors
                input_path_L2_ISR_bs(i_baseline) = cellstr(strcat(output_path_L2_ISR_bs(i_baseline), 'data', filesep));   
        end 
    end
end


%% ----------- FIND INPUT L2/L1B and config files -------------------------

fprintf('\nRunning waveform fitting tool....\n');

for i_baseline=1:num_baselines
    
    % --------- Define Paths -------------------------------------------------   
    filesBulk(i_baseline).input_path_L1_ISR_bs       =   char(input_path_L1_ISR_bs(i_baseline));
    filesBulk(i_baseline).resultPath      =   char(path_fitted_waveforms);
%     mkdir(filesBulk(i_baseline).resultPath);        
    
    % --------------- Configuration/characterization/LUTS --------------------
    filesBulk(i_baseline).cnf_chd_cst_path        =   cnf_chd_cst_path;
    inputFiles      =   dir(filesBulk(i_baseline).cnf_chd_cst_path);
    aux=struct2cell(inputFiles); aux=aux(1,:);
    if ~isempty(proc_bsln_id)
        filesBulk(i_baseline).CNF_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file','_',proc_bsln_id])))).name];
    else
        filesBulk(i_baseline).CNF_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file'])))).name];
    end
    filesBulk(i_baseline).CST_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cst_file','.json'])))).name];

    reduced_set_cnf.SCOOP=0;
    reduced_set_cnf.SHAPE=0;
    cnf_p=read_CNF_json(filesBulk(i_baseline).CNF_file,reduced_set_cnf);
    cnf_p.SCOOP_flag =0;
    cnf_p.SHAPE_flag =0;

    filesBulk(i_baseline).CHD_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['chd_file','_',cnf_p.mission,'.json'])))).name];
    cst_p=read_CST_json(filesBulk(i_baseline).CST_file); %CST --> provide the global constant variables
    chd_p=read_CHD_json(filesBulk(i_baseline).CHD_file,cnf_p,cst_p);

    idx_int=~cellfun(@isempty,strfind(aux,'LUT_f0'));
    if any(idx_int)
        filesBulk(i_baseline).LUT_f0_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(idx_int).name];
    end
    idx_int=~cellfun(@isempty,strfind(aux,'LUT_f1'));
    if any(idx_int)
        filesBulk(i_baseline).LUT_f1_file=[filesBulk(i_baseline).cnf_chd_cst_path inputFiles(idx_int).name];
    end
    
    % ------------------ Input L1B files filtering ---------------------------
    inputFiles      =   dir(filesBulk(i_baseline).input_path_L1_ISR_bs);
    aux=struct2cell(inputFiles); aux=aux(1,:); 
    aux_filter=aux(~strcmp(aux,['.'])&~strcmp(aux,['..']));
    accepted_filters = {'.NC', '.nc'};
    for ii=1:numel(aux_filter)
        [~,~,filter]=fileparts(aux_filter{ii});
        if any(find(ismember(accepted_filters,filter)))
            break;
        end
    end

%     filter='.NC';  
    filterDATAFILES=(~cellfun(@isempty,strfind(aux,filter)));
%     filterDATAFILES=(~cellfun(@isempty,cat(2, strfind(aux, '.nc'), strfind(aux, '.NC'))));
    indexFilesL1=find(filterDATAFILES);

    filesBulk(i_baseline).nFilesL1B=length(indexFilesL1);
    filesBulk(i_baseline).L1BFiles=inputFiles(indexFilesL1);
    
    % ------------------ Input L2 files filtering ---------------------------
    filesBulk(i_baseline).input_path_L2_ISR_bs        =   char(input_path_L2_ISR_bs{i_baseline});
    inputFiles      =   dir(filesBulk(i_baseline).input_path_L2_ISR_bs);
    aux=struct2cell(inputFiles); aux=aux(1,:);
    aux_filter=aux(~strcmp(aux,['.'])&~strcmp(aux,['..']));
    accepted_filters = {'.NC', '.nc'};
    for ii=1:numel(aux_filter)
        [~,~,filter]=fileparts(aux_filter{ii});
        if any(find(ismember(accepted_filters,filter)))
            break;
        end
    end
    
%     filter2DATAFILES=(~cellfun(@isempty,strfind(lower(aux),filter)));
    filter2DATAFILES=(~cellfun(@isempty,strfind(aux,filter)));
    indexFilesL2=find(filter2DATAFILES);        
    
    filesBulk(i_baseline).nFilesL2=length(indexFilesL2);
    filesBulk(i_baseline).L2Files=inputFiles(indexFilesL2);
    
    clear inputFiles indexFilesL2 filter2DATAFILES;    
    
    fprintf('Total number of L2 files (%s) to be processed: %.0f\n', char(name_bs(i_baseline)), filesBulk(i_baseline).nFilesL2);
    
end


%% ---------- PLOT FITTED CURVES and PERFORMANCE GRAPHICS --------------------

ref_SSH = 12;
ref_sigma0 = 11;
for i_fileL2_input=1:filesBulk(1).nFilesL2
    
    fprintf('\nReading L2 data of file no. %d...\n\n', i_fileL2_input);

    [SSH, SWH, sigma0, epoch, Pu, COR, filename_L2, filename_L1, index_inside_mask, flag]=read_L2(filesBulk, i_fileL2_input,  ...
        name_bs, cnf_tool, 'filename_mask_KML', filename_mask_KML);

    [~,filenopath,fileext_L1]=fileparts(filename_L1{1});

%       % find scenario from file name, e.g., OS20    
%     aux=strsplit(filenopath, '_');
%     SCENARIO = str2num(strcat(aux{2}(3), '0')); % assume OS10, OS20, ...
%     
%     switch SCENARIO
%         case {10,16}
%             ref_SWH = 1;
%         case {20,21,22,23,24}
%             ref_SWH = 2;
%         case 30
%             ref_SWH = 5;
%         case {40,41}
%             ref_SWH = 8;
%     end
    
    if flag==1        
%         if ~cnf_tool.run_L2 
        fprintf('\nPlotting fitted power waveforms...\n');

        model_fit_to_power_waveform(filesBulk, name_bs, cnf_p, cst_p, chd_p, filename_L1, SWH, sigma0, epoch, Pu, COR, cnf_tool);
%         model_fit_to_power_waveform_dynamic(filesBulk, name_bs, cnf_p, cst_p, chd_p, filename_L1, SWH, SSH, sigma0, epoch, Pu, COR, cnf_tool);

%         end

        % ---------------- RUN RETRIEVED PARAMETERS COMPARISON ------------------
        
        fprintf('\nPlotting comparison between retrieved parameters file %s...\n\n', filenopath);
        performance_baselines_S6(SSH, SWH, sigma0, COR, filename_L2, index_inside_mask, name_bs, path_comparison_results, cnf_tool,...
                        'win_size_detrending',win_size_detrending,...
                        'flag_outliers_removal',flag_outliers_removal,...
                        'type_outliers_removal',type_outliers_removal,...
                        'smooth_param',smooth_param,...
                        'generate_plot_SSH',generate_plot_SSH,'generate_plot_SWH',generate_plot_SWH,...
                        'generate_plot_sigma0',generate_plot_sigma0,...
                        'generate_plot_COR',generate_plot_COR)

    else
      fprintf('Missing files in some baseline \n');
    end
                    
end

fprintf('\nEnd of processing.\n')
                           
