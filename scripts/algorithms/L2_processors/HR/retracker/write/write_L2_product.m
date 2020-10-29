function write_L2_product(file,out,cnf_p,chd_p,cst_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and geophysical correctiosn to be applied to the
% retracked range based on the info available in the L1B product and
% depending on the surface being observed
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 17/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -file    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -out = structure of the output data for the L2 product
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       
% RESTRICTIONS: 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: No tracking changes until v2.0
% v2.0: Include explicitly the option to run analytical
% retracker in a 2step using cnf_p.analytical_type_of_fitting
%

%% ------------- NAME DEFINITION ------------------------------------------
%-------------- Define the output name for the output file ----------------
date_creation = datestr(now, '_yyyymmddTHHMMSS_');
aux=strsplit(date_creation,'_');
out.date_creation=aux(2);
clear aux;
switch cnf_p.mission    
    case 'CS2'
    %--------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                string_ID=file.filename_L1B_nopath(20:20+30);
                name_L2_product=strcat('CR2_SR_2_',cnf_p.Product_type,'____',string_ID,date_creation,'isd');
            case 'ISD'
                %following a la Sentinel-3 format as per SEOMs projects
                aux=file.filename_L1B_nopath(1:47);
                aux(10:15)=strcat(cnf_p.Product_type,'___');
                name_L2_product=strcat(strrep(aux(1:47),'_1_','_2_'),date_creation,'isd');
            case 'GPOD'
                string_ID=file.filename_L1B_nopath(17+7:47+7);
                name_L2_product = strcat('CR2_SR_2_',cnf_p.Product_type,'___',string_ID,date_creation,'isd');
        end
    case 'S3'
    %--------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                % TBD
                %[data]=readL1B_S3_ESA(filename_L1B);
                if strfind(file.filename_L1B_nopath,'_SR_1')
                  name_L2_product=strrep(file.filename_L1B_nopath(1:47),'_1_SRA',strcat('_2_',cnf_p.Product_type));
                elseif strfind(file.filename_L1B_nopath,'_SR_2')
                  name_L2_product=file.filename_L1B_nopath(1:47);
                  name_L2_product(10:12)=cnf_p.Product_type;
                end
            case 'DeDop'
%                 if strfind(file.filename_L1B_nopath,'_SR_1')
%                   name_L2_product=strrep(file.filename_L1B_nopath,'_SR_1',strcat('_SR_2_',cnf_p.Product_type));
%                 elseif strfind(file.filename_L1B_nopath,'_SR_2')
%                   name_L2_product=strrep(file.filename_L1B_nopath,'_SR_2',strcat('_SR_2_',cnf_p.Product_type));
%                 else
                    name_L2_product=strrep(strcat(file.filename_L1B_nopath),'L1B','L2');
%                 end
            case 'ISD'
                  name_L2_product=strcat(strrep(file.filename_L1B_nopath,'_1_SRA',strcat('_2_',cnf_p.Product_type)));
%                 aux=file.filename_L1B_nopath(1:47);
%                 aux(10:15)=cnf_p.Product_type+'___';
%                 name_L2_product=strcat(ststrrep(aux(1:47),'_1_','_2_'),date_creation,'isd');
        end
    case 'S6'
    %--------------------- Sentinel-6 -----------------------------------        
        switch cnf_p.L1proc
            case 'ISD'
                %needs to be defined as using the final file naming convention
                %name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','L2'),date_creation,'isd');
                name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','L2'),'_isd');
        end 
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end

if cnf_p.optional_ext_file_flag
    %name_L2_product=strrep(name_L2_product,'isd',strcat(cnf_p.file_ext_string,'_isd'));
    name_L2_product=strcat(name_L2_product,'_',cnf_p.file_ext_string);
end

%% ---------------------- GENERATION OF OUTPUT PRODUCT --------------------
switch cnf_p.output_product_format
    case 'mat' % Matlab output
        % .mat file
        save([file.resultPath 'data' filesep name_L2_product '.mat'],'out');        
    case 'nc'  % NetCDF output product
        switch cnf_p.nc_format_type
            case 'default'                
                % for all projects except SS_CCI
                prepare_NetCDF_L2([file.resultPath 'data' filesep name_L2_product '.nc'],out,cnf_p,cst_p,chd_p); %can have a switch to define the variables if mission is S3 or S6
            case 'SS_CCI'
                %in order to use specific format for SS CCI
                prepare_NetCDF_L2_SS_CCI([file.resultPath 'data' filesep name_L2_product '.nc'],out,cnf_p,cst_p,chd_p); %can have a switch to define the variables if mission is S3 or S6           
        end
    otherwise
        error('The type of product output is not supported')        
end

%% Include the processing options product
fid_proc_opts=fopen([file.resultPath 'data' filesep name_L2_product,'_proc_opt.txt'],'w');
fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
fprintf(fid_proc_opts,'--------------- General Parameters & Options retrackers ----------------------------\n');
fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering Geographical mask',{': '},num2str(cnf_p.mask_ROI_flag))));
fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering land surfaces',{': '},num2str(cnf_p.filter_land))));
if cnf_p.filter_land
    fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering land surfaces type',{': '},(cnf_p.filter_land_type))));
end
fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering Number looks stack',{': '},num2str(cnf_p.mask_looks_flag))));
if cnf_p.mask_looks_flag
    fprintf(fid_proc_opts,'%s\n',char(strcat('Minimum Number looks stack (filtering)',{': '},num2str(cnf_p.Neff_thres))));
end
fprintf(fid_proc_opts,'%s\n',char(strcat('IFmask',{': '},num2str(cnf_p.IFmask_N))));

%--------------------------------------------------------------------------
fprintf(fid_proc_opts,'%s\n',char(strcat('Discard waveform samples',{': '},num2str(cnf_p.wvfm_discard_samples))));
fprintf(fid_proc_opts,'%s\n',char(strcat('Number of samples discarded from begining',{': '},num2str(cnf_p.wvfm_discard_samples_begin))));
fprintf(fid_proc_opts,'%s\n',char(strcat('Number of samples discarded from end',{': '},num2str(cnf_p.wvfm_discard_samples_end))));

%--------------------------------------------------------------------------
fprintf(fid_proc_opts,'%s\n',char(strcat('Waveform portion selection',{': '},num2str(cnf_p.wvfm_portion_selec))));
if cnf_p.wvfm_portion_selec
    fprintf(fid_proc_opts,'%s\n',char(strcat('Portion selection type',{': '},(cnf_p.wvfm_portion_selec_type))));
    switch lower(cnf_p.wvfm_portion_selec_type)
        case 'ref_height'
            fprintf(fid_proc_opts,'%s\n',char(strcat('Reference DEM',{': '},cnf_p.wvfm_portion_selec_DEM_ref)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Path Reference DEM',{': '},cnf_p.wvfm_portion_selec_DEM_path)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Number of left samples from ref. position (height)',{': '},num2str(cnf_p.wvfm_portion_selec_l_samples))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Number of right samples from ref. position (height)',{': '},num2str(cnf_p.wvfm_portion_selec_r_samples))));
        case 'peak_win'            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Number of left samples from ref. position (peak)',{': '},num2str(cnf_p.wvfm_portion_selec_l_samples))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Number of right samples from ref. position (peak)',{': '},num2str(cnf_p.wvfm_portion_selec_r_samples))));
        case 'peak_thresh'            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold value for left samples (w.r.t peak)',{': '},num2str(cnf_p.wvfm_portion_selec_l_thres))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold value for right samples (w.r.t peak)',{': '},num2str(cnf_p.wvfm_portion_selec_r_thres))));            
    end
end
%--------------------------------------------------------------------------
fprintf(fid_proc_opts,'%s\n',char(strcat('Reference sample of on-board tracker',{': '},num2str(cnf_p.ref_sample_wd))));
%--------------------------------------------------------------------------
fprintf(fid_proc_opts,'%s\n',char(strcat('Geophysical corrections applied',{': '},num2str(cnf_p.geo_corr_application_flag))));
fprintf(fid_proc_opts,'%s\n',char(strcat('force_geocorr_surf_type',{': '},num2str(cnf_p.force_geocorr_surf_type))));
if cnf_p.force_geocorr_surf_type
    fprintf(fid_proc_opts,'%s\n',char(strcat('product_type_surface',{': '},cnf_p.product_type_surface)));
end
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of Analytical retracker -----------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Type of analytical retracker',{': '},char(cnf_p.analytical_type_of_fitting(i_index_analytical)))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Use zeros in multilooking',{': '},num2str(cnf_p.use_zeros_cnf))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Antenna compensation along-track',{': '},num2str(cnf_p.antenna_compensation_al))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Doppler mask construction',{': '},cnf_p.Doppler_mask_cons_option)));
            if strcmpi(cnf_p.Doppler_mask_cons_option,'internal')
                fprintf(fid_proc_opts,'%s\n',char(strcat('Doppler mask construction internal type',{': '},cnf_p.Doppler_mask_cons_internal)));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('Zero Padding in range',{': '},num2str(cnf_p.ZP))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('window_type_a',{': '},cnf_p.window_type_a)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('window_type_r',{': '},cnf_p.window_type_r)));
            switch cnf_p.window_type_r
                case {'Adaptive','Adaptive_S3'}
                    %disp(strcat('Range PTR approx:',{''},'Adaptive'));
                otherwise
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*chd_p.alpha_gr_chd))))));
            end
            switch cnf_p.window_type_a
                case {'Adaptive','Adaptive_S3'}
                    %disp(strcat('Range PTR approx:',{''},'Adaptive'));
                otherwise
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*chd_p.alpha_ga_chd))))));
            end            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Noise estimation method:',{''},cnf_p.Thn_estimation_method)));
            switch lower(cnf_p.Thn_estimation_method)
                case 'external'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('External mean noise power:',{''},num2str(cnf_p.external_Thn_value))));
                case 'fixed_window'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('First noise sample',{': '},num2str(cnf_p.Thn_w_first))));
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Width noise window',{': '},num2str(cnf_p.Thn_w_width))));
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Noise estimation over',{': '},num2str(cnf_p.Thn_ML_SL_method))));
                case 'adaptive'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold noise',{': '},num2str(cnf_p.threshold_noise))));
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Max numer of iterations',{': '},num2str(cnf_p.max_iter_noise))));
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Mult. factor noise thresh.',{': '},num2str(cnf_p.factor_increase_noise_iter))));                    
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Noise estimation over',{': '},num2str(cnf_p.Thn_ML_SL_method))));
            end
            
%             disp(strcat('Sign pitch',{': '},num2str(cnf_p.sign_pitch)))
%             disp(strcat('Sign roll',{': '},num2str(cnf_p.sign_roll)))
            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Value of fixed roughness for SWH fitting',{': '},num2str(cnf_p.rou))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Value of fixed SWH for MSS fitting',{': '},num2str(cnf_p.Hs))));
                
            fprintf(fid_proc_opts,'%s\n',char(strcat('Indexation method',{': '},(cnf_p.looks_index_method))));
            switch cnf_p.looks_index_method
                case 'Look_angle'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Look angle method',{': '},cnf_p.look_ang_method)));
                case 'Doppler_freq'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Doppler freq. method',{': '},cnf_p.fd_method)));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('Power wvfm model',{': '},cnf_p.power_wfm_model)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('LUT flag',{': '},num2str(cnf_p.lut_flag))));
            if cnf_p.pre_processing
                fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold (leading edge pre-processing)',{': '},num2str(cnf_p.percent_leading_edge))));
            end
            switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                case '2step'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Two step fitting COR threshold',{': '},num2str(cnf_p.two_step_fitting_COR_threshold_rou))));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('initial_param_fit_feedback_flag',{': '},num2str(cnf_p.initial_param_fit_feedback_flag))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_opt',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_opt_discard_std_threshold',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_size',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_size))));
            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Range indexation method',{': '},cnf_p.range_index_method)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Fitting method',{': '},cnf_p.fitting_fun_type)));
            opti_cell=strcat(fieldnames(cnf_p.fitting_options),{': '},cellfun(@num2str,struct2cell(cnf_p.fitting_options),'UniformOutput',false));
            for i_field=1:length(opti_cell)
                fprintf(fid_proc_opts,'%s\n',char(opti_cell(i_field)));
            end
            if ~isempty(cnf_p.fitting_options_lb)
                fprintf(fid_proc_opts,'%s\n',char(strcat('Lower bounds fitting: ',sprintf('%d',cnf_p.fitting_options_lb))));
            end
            if ~isempty(cnf_p.fitting_options_ub)
                fprintf(fid_proc_opts,'%s\n',char(strcat('Upper bounds fitting: ',sprintf('%d',cnf_p.fitting_options_ub))));
            end 
            
        case 'THRESHOLD'
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of Threshold retracker ------------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold value',{': '},num2str(cnf_p.th_retracker.percentage_peak))));

        case 'OCOG'
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of OCOG retracker -----------------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Percentage OCOG Amplitude',{': '},num2str(cnf_p.OCOG_retracker.percentage_pow_OCOG))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('First ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n1))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Last ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n2))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Offset',{': '},num2str(cnf_p.OCOG_retracker.offset))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('OCOG method',{': '},num2str(cnf_p.OCOG_retracker.implementation_method))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('OCOG param computation method',{': '},num2str(cnf_p.OCOG_retracker.param_comp_method))));
    end
end
fclose(fid_proc_opts);


end

