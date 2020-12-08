function exit_flag=L2_processing (filesBulk, i_fileL1B_input,cnf_p,...
                                cst_p,chd_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs the L2 processing
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1.2 21/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filesBulk    =   structure of input files within the folder where
%       to process the data (including the L1B as well as
%       configuration/characterization files):
%        filesBulk.inputPath        --> full path to the input L1B folder
%        filesBulk.resultPath       --> full path to results folder save L2 prod
%        filesBulk.inputPath_L1BS   --> full path to input L1BS folder
%        filesBulk.cnf_chd_cst_path --> full path to cnf/chd/cst/LUTs folder
%        filesBulk.CNF_file         --> full filename cnf config file
%        filesBulk.CHD_file         --> full filename chd charac file
%        filesBulk.CST_file         --> full filename cst const file
%        filesBulk.LUT_f0_file      --> full filename LUT for f0 function
%        filesBulk.LUT_f1_file      --> full filename LUT for f1 function
%        filesBulk.nFilesL1B        --> total number of L1B to be processed
%        filesBulk.nFilesL1BS       --> total number of L1BS to be processed
%        filesBulk.L1BFiles         --> information of the L1B files to be
%                                       processed: fields--> {name,date,bytes,isdir,datenum}
%        filesBulk.L1BSFiles        --> information of the L1BS files to be
%                                       processed: fields--> {name,date,bytes,isdir,datenum}
%       -i_fileL1B_input = index of the L1B within the filesBulk.LBFiles
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       Will write the corresponding L2 product in the specified resultPath
%       folder whenever the related flag cnf_p.write_output (n the cnf_file.m) is activated 
%       exit_flag = flag indicating whether the processing has been
%       successful or not o there is an error value to -1
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%   - preparing_inputFiles_L2processing: prepare and set the L1B, L1Bs and
%   KMl files when applicable from the input filesBulk data and cnf_p
%   - read_alt_data_EM: read input data from any combination of mission and
%   processor either in netcdf mat or DBL and save it in a common data
%   structure to be used in the L2 processing
%   - analytical_retracker: performs the analytical retracking based on the
%   original model developed by Chirs Ray et al. in IEEE TGRS "SAR
%   Altimeter Backscattered Waveform Model": DOI:
%   10.1109/TGRS.2014.23330423
%   - output_data_generation: acomodate the data of interest at the output
%   of the L2 processor in a given output structure which can be optionally 
%   saved (depending on the cnf_p.write_output) in .mat or  in netCDF 
%   depending on flag cnf_p.output_product_format
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
% 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 
% v1.1: Inclusion of option to run several retrackers: 
%   all the functions and processing related to the analytical or SAMOSA
%   retracker have been moved within a single function "analytical_retracker.m"
% v1.2: Accomodation and preparation of the input files into a separate
% auxiliar function "preparing_inputfiles_L2processing.m" in order to keep
% a more compact definition of L2_processing.m function
% v1.3: 28.03.2017: Include explicitly the option to run analytical
% retracker in a 2step using cnf_p.analytical_type_of_fitting (so we can run we the same cnf file the three options SWH fitting, MSS fitting and 2step for analytical)

try
time_init=tic;
% %% ------------ RUN global variables --------------------------------------
% %global cnf_p reduced_set_cnf
% % run(filesBulk.CNF_file); %CNF --> generate a cnf_p structure
% % run(filesBulk.CST_file); %CST --> provide the global constant variables
% % run(filesBulk.CHD_file); %CHD --> provide the global characterization variables
% 
% read_CNF_json(filesBulk.CNF_file); %CNF --> generate a cnf_p structure
% read_CST_json(filesBulk.CST_file); %CST --> provide the global constant variables
% read_CHD_json(filesBulk.CHD_file,cnf_p); %CHD --> provide the global characterization variables
% 
% if reduced_set_cnf
%     %force definition of cnf parameters not to be included in the CNF file
%     call_force_cnf_parameters_L2_GPP;
% end


%% ---------------- Handling input variables ------------------------------
if(nargin<5 || nargin>(5+7*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.addParamValue('input_path_L2_GPOD',{''}); 
p.addParamValue('attitude_extraction_GPOD',0);
p.addParamValue('Sig0AtmCorr_path',{''});
p.addParamValue('filename_dist_to_coast',{''});
p.parse(varargin{:});
filename_mask_KML=char(p.Results.filename_mask_KML);
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
input_path_L2_GPOD  = p.Results.input_path_L2_GPOD;
attitude_extraction_GPOD = p.Results.attitude_extraction_GPOD;
Sig0AtmCorr_path = p.Results.Sig0AtmCorr_path;
filename_dist_to_coast = char(p.Results.filename_dist_to_coast);
clear p;

exit_flag='';

%% ------------------ Ploting formating -----------------------------------
% set_default_plot;
% if cnf_p.visible_figures
%     set(0, 'DefaultFigureVisible', 'on');
% else
%     set(0, 'DefaultFigureVisible', 'off');
% end

%% ------------ Prepare/organize the input files definition ---------------
%--------------------------------------------------------------------------
[filename_L1B,filename_L1B_nopath,fileext_L1B,filename_L1BS,filename_L2,cnf_p,file_unique_id]=preparing_inputFiles_L2processing(filesBulk, i_fileL1B_input, cnf_p,...                                               
                                               'targz_option_active_L1B',targz_option_active_L1B,...
                                               'targz_option_active_L1BS',targz_option_active_L1BS);


disp(strcat('Processing L1B: ',filename_L1B_nopath,fileext_L1B))

%% --------------------- LOAD DATA FROM INPUT FILE L1B --------------------
%--------------------------------------------------------------------------
%reading the data from the specific L1B & load it in a specific "data"
%structure & specific filtering (geographical and/or number of looks is further applied)
[data,flag] = read_alt_data_EM (filename_L1B, cnf_p,cst_p,chd_p,...
    'filename_L1BS',filename_L1BS,'filename_L2',filename_L2,'filename_mask_KML',filename_mask_KML,...
                                'DEM_ref',cnf_p.wvfm_portion_selec_DEM_ref,'dir_DEM',cnf_p.wvfm_portion_selec_DEM_path,...
                                'path_Results',filesBulk.resultPath,...
                                'attitude_extraction_GPOD',attitude_extraction_GPOD,...
                                'input_path_L2_GPOD',input_path_L2_GPOD,...
                                'file_unique_id',file_unique_id,...
                                'Sig0AtmCorr_path',Sig0AtmCorr_path,...
                                'filename_dist_to_coast',filename_dist_to_coast);
%return;
if flag==-1
    %track outside the mask ROI
    Logtime = datestr(now, 'yyyymmddTHHMMSS');
    str1=sprintf('%s -> Error processing file %s\n',Logtime, filesBulk.L1BFiles(i_fileL1B_input).name);
    str2='Track outside the limits of the geographical mask';
    exit_flag=strcat(strcat('\n',repmat('-',[1,length(str1)]),'\n'),...
        str1,'\n',str2,...
        strcat('\n',repmat('-',[1,length(str1)]),'\n'));    
    return;
end

%% --------------------- RETRACKER RUN ------------------------------------
%--------------------------------------------------------------------------
time_init_retracking=tic;
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            disp(strcat('Processing with Analytical retracker: fitting',{' '},char(cnf_p.analytical_type_of_fitting(i_index_analytical)),'....'))                      
            
            %Patch to be able to run the analytical more than once (including the results in the same file) with
            %different options on the type of fitting : SWH, ROUGH or BOTH
            %Since we are using a single CNF file in this case for all of
            %them: other option is to create a different CNF for each
            %retracker            
                switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                    case 'SWH'
                        cnf_p.rou_flag=0;
                        cnf_p.two_step_fitting=0; %used to perform first a fitting on SWH and if
                    case 'MSS'
                        if cnf_p.SCOOP_flag==0
                            cnf_p.rou_flag=1;
                            cnf_p.two_step_fitting=0; %used to perform first a fitting on SWH and if
                        else
                            disp('This type of analytical retracker fitting with MSS not available for SCOOP');
                            continue;
                        end
                    case '2step'
                        %------------------ Two-step fitting ----------------------------------
                        if cnf_p.SCOOP_flag==0
                            cnf_p.rou_flag=0;
                            cnf_p.two_step_fitting=1; %used to perform first a fitting on SWH and if
                        else
                            disp('This type of analytical retracker with 2-step fitting not available for SCOOP');
                            continue;
                        end
                end
            [retracker_results.ANALYTICAL(i_index_analytical)]=analytical_retracker(data,cnf_p,...
                                                                chd_p,cst_p,...
                                                                'LUT_f0_file',filesBulk.LUT_f0_file,...
                                                                'LUT_f1_file',filesBulk.LUT_f1_file,...
                                                                'path_Results',filesBulk.resultPath,...
                                                                'L1B_filename',filename_L1B_nopath);
        case 'THRESHOLD'
            if cnf_p.SCOOP_flag==0                
                disp('Processing with Threshold retracker ....')
                [retracker_results.THRESHOLD]=threshold_retracker(data,cnf_p,chd_p,...
                    'path_Results',filesBulk.resultPath,...
                    'L1B_filename',filename_L1B_nopath);
            else
                disp('Threshold retracker not available for SCOOP');
                continue;
            end
        case 'OCOG'
            if cnf_p.SCOOP_flag==0
                disp('Processing with OCOG retracker ....')
                [retracker_results.OCOG]=OCOG_retracker(data,cnf_p,chd_p,...
                    'path_Results',filesBulk.resultPath,...
                    'L1B_filename',filename_L1B_nopath);
            else
                disp('OCOG retracker not available for SCOOP');
                continue;
            end
    end
end
time_end=toc(time_init_retracking);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Retracking time for ',filename_L1B_nopath,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);


%% -------------- ORGANIZE THE OUTPUT DATA --------------------------------
%--------------------------------------------------------------------------
%generating an output data structure and writing the corresponding product
%(either using .nc or .mat according to the configuration)
file.filename_L1B_nopath=filename_L1B_nopath;
file.fileext_L1B=fileext_L1B;
file.inputPath=filesBulk.inputPath;
file.resultPath=filesBulk.resultPath;
[out]=output_data_generation(file,retracker_results,data,cnf_p,chd_p,cst_p);


time_end=toc(time_init);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Processing time for ',file.filename_L1B_nopath,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);

catch err
    Logtime = datestr(now, 'yyyymmddTHHMMSS');
    str1=sprintf('%s -> Error processing file %s\n',Logtime, filesBulk.L1BFiles(i_fileL1B_input).name);
    str2=sprintf('%s\n', err.getReport('extended', 'hyperlinks','off'));
    str_err=strcat(strcat('\n',repmat('-',[1,length(str1)]),'\n'),...
        str1,'\n',str2,...
        strcat('\n',repmat('-',[1,length(str1)]),'\n'));
    %disp(sprintf(str_err));
    exit_flag=str_err;
    return;
    
end
    


end

