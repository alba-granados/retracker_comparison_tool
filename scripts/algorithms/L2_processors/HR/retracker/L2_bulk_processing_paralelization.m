function L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs a bulk L2 processing either in parallel or sequentially
% over the different input L1B files defined in a given folder, using a given
% set of configuration and characterization files provided in a given
% folder, storing the results in a specified folder
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 20/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -input_path_L1B    =   path containing the L1B file/s to be
%       processed either in parallel or sequentially
%       -output_path_L2 = output folder to save the different L2 products
%       (will be written whenever the corresponding flag cnf_p.write_output (in the cnf_file.m) is activated )
%       -cnf_chd_cst_path: path containing the different configuration
%       cnf_file.m
%       (configuration processing parameters for L2), characterization
%       chd_file.m
%       (different missions parameters) and the constant definition
%       cst_file.m; the Look Up Tables (LUTs) for the f0 and f1 functions
%       to be used by the single look waveform model are also included in this folder
%     OPTIONAL
%       -input_path_L1BS: path to the related L1BS products associated to
%       L1B products in input_path_L1B
%       -filename_mask_KML: fullfilename of a KML file containing the ROI
%       in order to filter out the records of interest: it is used only if
%       the related flag is activated cnf_p.mask_ROI_flag (in cnf_file.m)
%       -num_pools: indicates the number of pools or threads to be run for
%       parallel processing of different tracks: by default is set to 1 if
%       no paralelization is desired (sequential processing)
%       -targz_option_active_L1B: flag to indicate whether the available
%       L1B products within the folder to be processed are compressed and need
%       to be untar
%      -targz_option_active_L1BS: flag to indicate whether the available
%      L1BS products within the folder to be processed are compressed and need
%       to be untar
%       
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%  - L2_processing: runs the L2 processing; geophysical retrievals
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% RESTRICTIONS:
% - It is assumed all the L1B input files in the inputh_path_L1B correspond
% to the same mission and have been processed with the same configuration
% parameters
% - LUTs files are saved as .mat files, need to be updated for a more
% generic case of binary or netcdf files
% - The optional L1B-S files included in a different folder 
% - Further updates to have the possibility of having compressed and
% uncompressed files in the same folder
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 

clear global
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');

%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+12*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('input_path_L1BS',{''},@(x)ischar(x));
p.addParamValue('input_path_L2',{''},@(x)ischar(x));
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('num_pools',1);
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.addParamValue('proc_bsln_id',{''},@(x)ischar(x));
p.addParamValue('log_flag',1);
p.addParamValue('filter_input_file','.nc');
p.addParamValue('PTR_along_external',[]); %reference SSH to compute error for S6 sim
p.addParamValue('PTR_across_external',[]); %reference SWH to compute error for S6 sim
p.addParamValue('input_path_L2_GPOD',{''}); 
p.addParamValue('attitude_extraction_GPOD',0);
p.addParamValue('Sig0AtmCorr_path',{''});
p.addParamValue('filename_dist_to_coast',{''});
p.parse(varargin{:});
input_path_L1BS=char(p.Results.input_path_L1BS);
input_path_L2=char(p.Results.input_path_L2);
filename_mask_KML=char(p.Results.filename_mask_KML);
num_pools=p.Results.num_pools;
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
proc_bsln_id=char(p.Results.proc_bsln_id);
log_flag=p.Results.log_flag;
filter_input_file=p.Results.filter_input_file;
PTR_along_external  = p.Results.PTR_along_external;
PTR_across_external = p.Results.PTR_across_external;
input_path_L2_GPOD  = p.Results.input_path_L2_GPOD;
attitude_extraction_GPOD = p.Results.attitude_extraction_GPOD;
Sig0AtmCorr_path = p.Results.Sig0AtmCorr_path;
filename_dist_to_coast = char(p.Results.filename_dist_to_coast);
clear p;

%define the type of configurations available for L2 processor (set to 1
%only one of the 2 at time): to be used for generating the executable file
%from matlab (delivered to ESA) so we can use a reduced version of cnf file with reduced
%number of options (otherwise forced): if to use SCOOP: set
%reduced_set_cnf.SCOOP =1 and  reduced_set_cnf.SHAPE=0; and viceversa for
%SHAPE; normal operation would be to set reduced_set_cnf.SCOOP=0 and
%reduced_set_cnf.SHAPE=0; ONLY WHEN GENERATING THE .EXE MATLAB TO BE
%DELIVERED TO ESA
reduced_set_cnf.SCOOP=0;
reduced_set_cnf.SHAPE=0

version_matlab=version;
ttotal=tic;

    if(log_flag)
        str_err={''};
    end

%% ----- Include the different subfolders in the search path --------------
FolderInfo=dir;
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat(['.' filesep] ,folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end


%% --------- Define Paths -------------------------------------------------
filesBulk.inputPath       =   input_path_L1B;
filesBulk.resultPath      =   output_path_L2;
mkdir(filesBulk.resultPath);
mkdir([filesBulk.resultPath 'data' filesep]);
mkdir([filesBulk.resultPath 'plots' filesep]);
mkdir([filesBulk.resultPath 'plots' filesep 'fitted_waveforms' filesep]);




%% --------------- Configuration/characterization/LUTS --------------------
%assume all the files in the folder correspond to the same mission and
%processed with the same L1B processor: use the same
%configuration/characterization for Level-2 processing
filesBulk.cnf_chd_cst_path        =   cnf_chd_cst_path;
inputFiles      =   dir(filesBulk.cnf_chd_cst_path);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
if ~isempty(proc_bsln_id)
    filesBulk.CNF_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file','_',proc_bsln_id])))).name];
else
    filesBulk.CNF_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file'])))).name];
end
filesBulk.CST_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cst_file','.json'])))).name];
%filesBulk.CHD_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,'chd_file'))).name];
% run(filesBulk.CNF_file); %CNF --> generate a cnf_p structure
% run(filesBulk.CHD_file); %CHD --> provide the global characterization variables

%new read configuration files based on json files
cnf_p=read_CNF_json(filesBulk.CNF_file,reduced_set_cnf);
% Force specific parameters of the configuration as not all them available
%at the input configuration file
if reduced_set_cnf.SCOOP
    %&&&&&&&&&&&&&&&&&&&& Force some cnf parameters &&&&&&&&&&&&&&&&&&&&&&&&&&&
    call_force_cnf_parameters_L2_GPP_SCOOP; %(not to output all variables)
    %flag in cnf structure
    cnf_p.SCOOP_flag =1;
    cnf_p.SHAPE_flag =0;
elseif  reduced_set_cnf.SHAPE
    call_force_cnf_parameters_L2_GPP_SHAPE; %(not to output all variables)
    %flag in cnf structure
    cnf_p.SCOOP_flag =0;
    cnf_p.SHAPE_flag =1;
else
    cnf_p.SCOOP_flag =0;
    cnf_p.SHAPE_flag =0;
end
%--------------- constants file -------------------------------------------
cst_p=read_CST_json(filesBulk.CST_file); %CST --> provide the global constant variables
% ------------- characterization file -------------------------------------
if  ~strcmp(cnf_p.mode,'SARin')
    filesBulk.CHD_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['chd_file','_',cnf_p.mission,'.json'])))).name];
else
    filesBulk.CHD_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['chd_file','_',cnf_p.mission,'_SIN','.json'])))).name];
end
chd_p=read_CHD_json(filesBulk.CHD_file,cnf_p,cst_p,'PTR_along_external',PTR_along_external,'PTR_across_external',PTR_across_external);





%change name extension to define PTR externals
if ~isempty(PTR_along_external)
    cnf_p.file_ext_string='';
    cnf_p.optional_ext_file_flag=1;
    cnf_p.file_ext_string = strcat(cnf_p.file_ext_string,'_PTR_AL_',num2str(PTR_along_external,'%08f'));
end    

if ~isempty(PTR_across_external)
    if isempty(PTR_along_external)        
        cnf_p.file_ext_string='';
    end
    cnf_p.optional_ext_file_flag=1;
    cnf_p.file_ext_string = strcat(cnf_p.file_ext_string,'_PTR_AC_',num2str(PTR_across_external,'%08f'));
end


%---------------------- LUTs files-----------------------------------------
%assuming .mat files
idx_int=~cellfun(@isempty,strfind(aux,'LUT_f0'));
if any(idx_int)
    filesBulk.LUT_f0_file=[filesBulk.cnf_chd_cst_path inputFiles(idx_int).name];
end
idx_int=~cellfun(@isempty,strfind(aux,'LUT_f1'));
if any(idx_int)
    filesBulk.LUT_f1_file=[filesBulk.cnf_chd_cst_path inputFiles(idx_int).name];
end


%% ------------------ Input L1B files filtering ---------------------------
inputFiles      =   dir(filesBulk.inputPath);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the

if targz_option_active_L1B
    filterDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
    indexFilesL1B=find(filterDATAFILES);
else
    switch cnf_p.mission
        case {'CS2','CR2'}
            % --------------------- CroySAT-2 ---------------------------------
            switch cnf_p.L1proc
                case {'ESA'}
                    filter='.dbl';
                case {'ISD','GPOD'}
                    filter='.nc';
            end
        case {'S6','JCS'}
            % --------------------- Sentinel-6 -----------------------------            
            filter=filter_input_file;%'.nc';
        case {'S3A','S3B','S3'}
            switch cnf_p.L1proc
                case 'ESA'
                    filter='.sen3'; %each netcdf is within a folder ending name as .SEN3
                case {'ISD','DeDop'}
                    filter='.nc';
            end
        otherwise
            error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
    end
    filterDATAFILES=(~cellfun(@isempty,strfind(lower(aux),filter)));
    indexFilesL1B=find(filterDATAFILES);
    
end

filesBulk.nFilesL1B=length(indexFilesL1B);
filesBulk.L1BFiles=inputFiles(indexFilesL1B);

disp('Total number of L1B files to be processed');
disp(num2str(filesBulk.nFilesL1B))


%% --------------- Input L1B-S files --------------------------------------
if ~isempty(input_path_L1BS)
    filesBulk.inputPath_L1BS        =   input_path_L1BS;
    inputFiles      =   dir(filesBulk.inputPath_L1BS);
    aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
    
    if targz_option_active_L1BS
        filterL1BSDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
        indexFilesL1BS=find(filterL1BSDATAFILES);
    else
        filterL1BSDATAFILES=(~cellfun(@isempty,strfind(aux,filter)));
        indexFilesL1BS=find(filterL1BSDATAFILES);        
    end
    
    filesBulk.nFilesL1BS=length(indexFilesL1BS);
    filesBulk.L1BSFiles=inputFiles(indexFilesL1BS);
    clear inputFiles indexFilesL1BS filterL1BSDATAFILES;    
end

%% --------------- Input L2 files --------------------------------------
if ~isempty(input_path_L2)
    filesBulk.inputPath_L2        =   input_path_L2;
    inputFiles      =   dir(filesBulk.inputPath_L2);
    aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
    
    filter2DATAFILES=(~cellfun(@isempty,strfind(lower(aux),filter)));
    indexFilesL2=find(filter2DATAFILES);        

    
    filesBulk.nFilesL2=length(indexFilesL2);
    filesBulk.L2Files=inputFiles(indexFilesL2);
    clear inputFiles indexFilesL2 filter2DATAFILES;    
end


%% --------------- Displaying options of retrackers -----------------------
disp('------------------------------------------------------------------------------------')
disp('--------------- General Parameters & Options retrackers ----------------------------')
disp('------------------------------------------------------------------------------------')
disp(strcat('Filtering Geographical mask',{': '},num2str(cnf_p.mask_ROI_flag)))
disp(strcat('Filtering land surfaces',{': '},num2str(cnf_p.filter_land)));
if cnf_p.filter_land
    disp(strcat('Filtering land surfaces type',{': '},(cnf_p.filter_land_type)));
end
disp(strcat('Filtering Number looks stack',{': '},num2str(cnf_p.mask_looks_flag)))
if cnf_p.mask_looks_flag
    disp(strcat('Minimum Number looks stack (filtering)',{': '},num2str(cnf_p.Neff_thres)))
end
disp(strcat('IFmask',{': '},num2str(cnf_p.IFmask_N)))

disp(strcat('Discard waveform samples',{': '},num2str(cnf_p.wvfm_discard_samples)))
disp(strcat('Number of samples discarded from begining',{': '},num2str(cnf_p.wvfm_discard_samples_begin)))
disp(strcat('Number of samples discarded from end',{': '},num2str(cnf_p.wvfm_discard_samples_end)))

disp(strcat('Waveform portion selection',{': '},num2str(cnf_p.wvfm_portion_selec)))
if cnf_p.wvfm_portion_selec
    disp(strcat('Portion selection type',{': '},(cnf_p.wvfm_portion_selec_type)))
    switch lower(cnf_p.wvfm_portion_selec_type)
        case 'ref_height'
            disp(strcat('Reference DEM',{': '},cnf_p.wvfm_portion_selec_DEM_ref))
            disp(strcat('Path Reference DEM',{': '},cnf_p.wvfm_portion_selec_DEM_path))
            disp(strcat('Number of left samples from ref. position (height)',{': '},num2str(cnf_p.wvfm_portion_selec_l_samples)))
            disp(strcat('Number of right samples from ref. position (height)',{': '},num2str(cnf_p.wvfm_portion_selec_r_samples)))
        case 'peak_win'            
            disp(strcat('Number of left samples from ref. position (peak)',{': '},num2str(cnf_p.wvfm_portion_selec_l_samples)))
            disp(strcat('Number of right samples from ref. position (peak)',{': '},num2str(cnf_p.wvfm_portion_selec_r_samples)))
        case 'peak_thresh'            
            disp(strcat('Threshold value for left samples (w.r.t peak)',{': '},num2str(cnf_p.wvfm_portion_selec_l_thres)))
            disp(strcat('Threshold value for right samples (w.r.t peak)',{': '},num2str(cnf_p.wvfm_portion_selec_r_thres)))
    end
end
disp(strcat('Reference sample of on-board tracker',{': '},num2str(cnf_p.ref_sample_wd)))
disp(strcat('Geophysical corrections applied',{': '},num2str(cnf_p.geo_corr_application_flag)))
disp(strcat('force_geocorr_surf_type',{': '},num2str(cnf_p.force_geocorr_surf_type)))
disp(strcat('atm_att_correction_flag',{': '},num2str(cnf_p.atm_att_correction_flag)))
if cnf_p.force_geocorr_surf_type
    disp(strcat('product_type_surface',{': '},cnf_p.product_type_surface))
end
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of Analytical retracker -----------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Type of analytical retracker',{': '},char(cnf_p.analytical_type_of_fitting(i_index_analytical))))
            disp(strcat('Use zeros in multilooking',{': '},num2str(cnf_p.use_zeros_cnf)))
            disp(strcat('Zero Padding in range',{': '},num2str(cnf_p.ZP)))
            disp(strcat('Antenna compensation along-track',{': '},num2str(cnf_p.antenna_compensation_al)))
            disp(strcat('Doppler mask construction',{': '},cnf_p.Doppler_mask_cons_option));
            if strcmpi(cnf_p.Doppler_mask_cons_option,'internal')
                disp(strcat('Doppler mask construction internal type',{': '},cnf_p.Doppler_mask_cons_internal));
            end
            disp(strcat('window_type_a',{': '},cnf_p.window_type_a));
            disp(strcat('window_type_r',{': '},cnf_p.window_type_r));
            switch cnf_p.window_type_r
                case {'Adaptive','Adaptive_S3'}
                    disp(strcat('Range PTR approx:',{''},'Adaptive'));
                otherwise
                    disp(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*chd_p.alpha_gr_chd)))))
            end
            switch cnf_p.window_type_a
                case {'Adaptive','Adaptive_S3'}
                    disp(strcat('Azimuth PTR approx:',{''},'Adaptive'));
                otherwise
                    disp(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*chd_p.alpha_ga_chd)))))        
            end
            
            disp(strcat('Noise estimation method:',{''},cnf_p.Thn_estimation_method))
            switch lower(cnf_p.Thn_estimation_method)
                case 'external'
                    disp(strcat('External mean noise power:',{''},num2str(cnf_p.external_Thn_value)))
                case 'fixed_window'
                    disp(strcat('First noise sample',{': '},num2str(cnf_p.Thn_w_first)))
                    disp(strcat('Width noise window',{': '},num2str(cnf_p.Thn_w_width)))
                    disp(strcat('Noise estimation over',{': '},num2str(cnf_p.Thn_ML_SL_method)))
                case 'adaptive'
                    disp(strcat('Threshold noise',{': '},num2str(cnf_p.threshold_noise)))
                    disp(strcat('Max numer of iterations',{': '},num2str(cnf_p.max_iter_noise)));
                    disp(strcat('Mult. factor noise thresh.',{': '},num2str(cnf_p.factor_increase_noise_iter)));                    
                    disp(strcat('Noise estimation over',{': '},num2str(cnf_p.Thn_ML_SL_method)))
            end
            
%             disp(strcat('Sign pitch',{': '},num2str(cnf_p.sign_pitch)))
%             disp(strcat('Sign roll',{': '},num2str(cnf_p.sign_roll)))
            
            disp(strcat('Value of fixed roughness for SWH fitting',{': '},num2str(cnf_p.rou)))
            disp(strcat('Value of fixed SWH for MSS fitting',{': '},num2str(cnf_p.Hs)))

            disp(strcat('Indexation method',{': '},(cnf_p.looks_index_method)))
            switch cnf_p.looks_index_method
                case 'Look_angle'
                    disp(strcat('Look angle method',{': '},cnf_p.look_ang_method))
                case 'Doppler_freq'
                    disp(strcat('Doppler freq. method',{': '},cnf_p.fd_method))
            end
            disp(strcat('Power wvfm model',{': '},cnf_p.power_wfm_model))
            disp(strcat('LUT flag',{': '},num2str(cnf_p.lut_flag)))
            if cnf_p.pre_processing
                disp(strcat('Threshold (leading edge pre-processing)',{': '},num2str(cnf_p.percent_leading_edge)))
            end
            switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                case '2step'
                    disp(strcat('Two step fitting COR threshold',{': '},num2str(cnf_p.two_step_fitting_COR_threshold_rou)))
            end
            disp(strcat('initial_param_fit_feedback_flag',{': '},num2str(cnf_p.initial_param_fit_feedback_flag)))
            disp(strcat('ini_Hs_rou_sliding_win_opt',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt)))
            disp(strcat('ini_Hs_rou_sliding_win_opt_discard_std_threshold',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold)))
            disp(strcat('ini_Hs_rou_sliding_win_size',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_size)))
            disp(strcat('Range indexation method',{': '},cnf_p.range_index_method))
            disp(strcat('Fitting method',{': '},cnf_p.fitting_fun_type))
            cnf_p.fitting_options
            
        case 'THRESHOLD'
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of Threshold retracker ------------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Threshold value',{': '},num2str(cnf_p.th_retracker.percentage_peak)))

        case 'OCOG'
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of OCOG retracker -----------------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Percentage OCOG Amplitude',{': '},num2str(cnf_p.OCOG_retracker.percentage_pow_OCOG)))
            disp(strcat('First ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n1)))
            disp(strcat('Last ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n2)))
            disp(strcat('Offset',{': '},num2str(cnf_p.OCOG_retracker.offset)))            
            disp(strcat('OCOG method',{': '},num2str(cnf_p.OCOG_retracker.implementation_method)))
            disp(strcat('OCOG param computation method',{': '},num2str(cnf_p.OCOG_retracker.param_comp_method)))
    end
end

%% --------------- Run parallel processing --------------------------------
if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);
    else
        matlabpool('open',num_pools);
    end    
    
    %% ------------- Loop per file to be processed ------------------------
    parfor i_fileL1B_input=1:filesBulk.nFilesL1B
        try
            exit_flag=L2_processing (filesBulk, i_fileL1B_input,...
                cnf_p,cst_p,chd_p,...
                'filename_mask_KML',filename_mask_KML,...
                'targz_option_active_L1B',targz_option_active_L1B,...
                'targz_option_active_L1BS',targz_option_active_L1BS,...
                'attitude_extraction_GPOD',attitude_extraction_GPOD,...
                'input_path_L2_GPOD',input_path_L2_GPOD,...
                'Sig0AtmCorr_path',Sig0AtmCorr_path,...
                'filename_dist_to_coast',filename_dist_to_coast);
            if(log_flag)
                str_err{i_fileL1B_input}=exit_flag;
            end
            
        catch err
            if(log_flag)
                Logtime = datestr(now, 'yyyymmddTHHMMSS');
                str1=sprintf('%s -> Error processing file %s\n',Logtime, filesBulk.L1BFiles(i_fileL1B_input).name);
                str2=sprintf('%s\n', err.getReport('extended', 'hyperlinks','off'));
                str_err{i_fileL1B_input}=strcat(strcat('\n',repmat('-',[1,length(str1)]),'\n'),...
                                                str1,'\n',str2,...
                                                strcat('\n',repmat('-',[1,length(str1)]),'\n'))
               %disp(char(str_err{i_fileL1B_input}));
            end            
        end
    end
    
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end
else
    for i_fileL1B_input=1:filesBulk.nFilesL1B
        try
            exit_flag=L2_processing (filesBulk, i_fileL1B_input,...
                cnf_p,cst_p,chd_p,...
                'filename_mask_KML',filename_mask_KML,...
                'targz_option_active_L1B',targz_option_active_L1B,...
                'targz_option_active_L1BS',targz_option_active_L1BS,...
                'attitude_extraction_GPOD',attitude_extraction_GPOD,...
                'input_path_L2_GPOD',input_path_L2_GPOD,...
                'Sig0AtmCorr_path',Sig0AtmCorr_path,...
                'filename_dist_to_coast',filename_dist_to_coast);
            if(log_flag)
                str_err{i_fileL1B_input}=exit_flag;
            end
        catch err
            if(log_flag)
                Logtime = datestr(now, 'yyyymmddTHHMMSS');
                str1=sprintf('%s -> Error processing file %s\n',Logtime, filesBulk.L1BFiles(i_fileL1B_input).name);
                str2=sprintf('%s\n', err.getReport('extended', 'hyperlinks','off'));
                str_err{i_fileL1B_input}=strcat(strcat('\n',repmat('-',[1,length(str1)]),'\n'),...
                                                str1,'\n',str2,...
                                                strcat('\n',repmat('-',[1,length(str1)]),'\n'))
               %disp(char(str_err{i_fileL1B_input}));
            end 
            continue;
        end
    end
end

time = toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time: ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);

%close the error annotation file
if(log_flag)

    fid_log = fopen([output_path_L2 'data' filesep 'LogError.txt'],'w');

    if any(~cellfun(@isempty,str_err))
        fprintf(fid_log,strjoin(str_err,'\n \n'));
    else
        fprintf(fid_log,'\n No errors \n');
    end
    fprintf(fid_log,strcat('Total processing time: ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds \n \n'));
    fclose (fid_log);
end

%exit
end
    

