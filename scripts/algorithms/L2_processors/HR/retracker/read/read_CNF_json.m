%function to read the configuration file in JSON and adapted to MAT struct
function [cnf_p]=read_CNF_json(cnf_file,reduced_set_cnf)
    %use global variable structure instead of return data
    %global cnf_p
    %global reduced_set_cnf
    struct=loadjson(cnf_file);
    
    version_matlab=version;
    
    %% ------------------------- MISSION & MODE -------------------------------
    % -------------------------------------------------------------------------
    if reduced_set_cnf.SCOOP==0
        cnf_p.mission           =   struct.mission.value;  % {'CS2','S3','S6'}
        cnf_p.mode              =   struct.mode.value;      % up to now only SAR mode operates: this shall be SAR, SARin --> as in CS2 or HRM_raw(2) and HRM_rmc(3) for Jason-CS/S6
    end
    
    
    %% ------------------------- L1 PROCESSOR ---------------------------------
    % -------------------------------------------------------------------------
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        cnf_p.L1proc            =   struct.L1proc.value;   %{'ESA','ISD'}
        % The L1 processor shall be
        %   ESA --> ESA processing baseline for the different type of missions
        %   ISD --> isardSAT processing baseline based on Sentinel-6 GPP
    end
    
    %--------------------------------------------------------------------------
    %--------------------- L1B processing characteristics ---------------------
    %--------------------------------------------------------------------------
    %zero-padding in range
    cnf_p.ZP                =   struct.ZP.value;              % 1 if no ZP, 2 if ZP of 2, 4 if ZP of 4, etc...
    %zeros in ML
    cnf_p.use_zeros_cnf     =   struct.use_zeros_cnf.value;
    %antenna pattern compensation along-track
    if reduced_set_cnf.SHAPE==0
        cnf_p.antenna_compensation_al = struct.antenna_compensation_along.value;
    end
    
    cnf_p.window_type_a     =   struct.window_type_a.value;       % ALONG TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)
    cnf_p.window_type_r     =   struct.window_type_r.value;       % ACROSS TRACK window type previous to AFFT can be Hamming, Hanning or Boxcar (Note: use capital letters)
    
    %read the forced values of the PTR from the configuration file itself
    switch cnf_p.window_type_a
        case 'Forced'
            cnf_p.A_s2Ga_chd=1;
            cnf_p.alpha_ga_chd=1.0/(2.0*(struct.PTR_forced_al.value).^2);
    end
    
    switch cnf_p.window_type_r
        case 'Forced'
            cnf_p.A_s2Gr_chd=1;
            cnf_p.alpha_gr_chd=1.0/(2.0*(struct.PTR_forced_ac.value).^2);
    end

    
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        %--------------------------------------------------------------------------
        %----------------------- SEED INFORMATION ---------------------------------
        %--------------------------------------------------------------------------
        cnf_p.seed              =   struct.seed_info.value;              % 0 - no seed introduced; 1 - seed provided by phase information of SARin data
        
        
        %--------------------------------------------------------------------------
        % -------------------- IF MASK --------------------------------------------
        %--------------------------------------------------------------------------
        cnf_p.IFmask_N          =   struct.IFmask_N.value;              % First/last N samples of the waveform affected by the IF Mask filter which should be excluded from processing. Specify if known, 0 if not known
        
        
        %------------------------------------------------------------------
        %------------------- FILTER LAND ----------------------------------
        %------------------------------------------------------------------
        cnf_p.filter_land       =   struct.filter_land.value;
        cnf_p.filter_land_type  =   struct.filter_land_type.value;
        
    end
    
    
    %% ------------------------ L2 PROCESSOR ----------------------------------
    if str2double(version_matlab(end-5:end-2))>2015
        cnf_p.retracker_name=cellstr(char(string(struct.retracker_name.value)));%{'THRESHOLD','OCOG','ANALYTICAL','ANALYTICAL','ANALYTICAL'}; %vector to run the different retrackers on the same data:
    else
        cnf_p.retracker_name=(((struct.retracker_name.value)));%{'THRESHOLD','OCOG','ANALYTICAL','ANALYTICAL','ANALYTICAL'}; %vector to run the different retrackers on the same data:
    end
    %'SAMOSA' or 'ANALYTICAL'--> analytical retracker
    %'THRESHOLD' --> Simple percentage of maximum of peak
    %'OCOG' --> Offset Center of Gravity
    %Rough approach when calling the ANALYTICAL more than once to differentiate
    %them whether fitting is on SWH, ROUGHNESS (MSS) or both of them (still to be implemneted)
    if str2double(version_matlab(end-5:end-2))>2015
        cnf_p.analytical_type_of_fitting=cellstr(char(string(struct.analytical_type_of_fitting.value))); %common structure
    else
        cnf_p.analytical_type_of_fitting=(((struct.analytical_type_of_fitting.value))); %common structure
    end
    
    %-------------------------------------------------------------------------
    %----------------------- REF. SAMPLES ------------------------------------
    %-------------------------------------------------------------------------
    %this is the reference sample for the window delay in L1B
    cnf_p.ref_sample_wd=struct.ref_sample_wd.value; % 64 in CS-2 no ZP, with ZP is 64*ZP factor and 43/42 in S-3
    
    
    %--------------------------------------------------------------------------
    % -------------------- FILTER BY GEOGRAPHICAL LOCATION --------------------
    %--------------------------------------------------------------------------
    cnf_p.mask_ROI_flag               =   struct.mask_ROI_flag.value;              % 0 processe all track; Otherwise specify a KML file
    
    %--------------------------------------------------------------------------
    % -------------------- FILTER BY # LOOKS/BEAMS ----------------------------
    %--------------------------------------------------------------------------
    % define the minimum number of looks available per waveform to consider it
    % for the processing --> we exclude those waveforms not generated with 95% of the looks available in the mode theoretically
    cnf_p.mask_looks_flag=struct.mask_looks_flag.value;
    cnf_p.Neff_thres=struct.Neff_thres.value;
    
    
    %--------------------------------------------------------------------------
    %--------------------- WAVEFORM PORTION SELECTION -------------------------
    %--------------------------------------------------------------------------
    if reduced_set_cnf.SCOOP==0
        %activate the filtering to select the part of waveform of interest
        cnf_p.wvfm_portion_selec=struct.wvfm_portion_selec.value;
        cnf_p.wvfm_portion_selec_type=struct.wvfm_portion_selec_type.value; % -'ref_height': using a reference height from DEM
        % -'peak_win': using the
        % maximum peak of waveform
        % with a window
        % around it (wvfm_portion_selec_l_samples, wvfm_portion_selec_r_samples)
        %-'peak_valley': using the
        %maximum peak selecting the
        %samples around it based on
        %two closest valleys srounding
        %the peak plus a margin
        %defined by the
        %wvfm_portion_selec_l_samples
        %and
        %wvfm_portion_selec_r_samples
        cnf_p.wvfm_portion_selec_DEM_ref=struct.wvfm_portion_selec_DEM_ref.value;
        cnf_p.wvfm_portion_selec_DEM_path=struct.wvfm_portion_selec_DEM_path.value; %path to the DEM to be used or saved in case of SRTM
        cnf_p.wvfm_portion_selec_l_samples=struct.wvfm_portion_selec_l_samples.value; % left samples w.r.t to position of ref. position within window (either peak or valley)
        cnf_p.wvfm_portion_selec_r_samples=struct.wvfm_portion_selec_r_samples.value; % right samples w.r.t to position of ref. position within window (either peak or valley)
        cnf_p.peak_prominence_norm=struct.peak_prominence_norm.value; % define the minimum prominence of the peaks to be sorted out in a multipeak scenario
        % currently used only with ref_height
        % option
    end
    %----------------------------------------------------------------------
    %------------------ Discarding samples begin/end waveform -------------
    %----------------------------------------------------------------------
    if reduced_set_cnf.SHAPE==0
        cnf_p.wvfm_discard_samples       = struct.wvfm_discard_samples.value;
        cnf_p.wvfm_discard_samples_begin = struct.wvfm_discard_samples_begin.value;
        cnf_p.wvfm_discard_samples_end   = struct.wvfm_discard_samples_end.value;
    end
    
    
    
    
    %% ----------------- THRESHOLD-BASED RETRACKER ----------------------------
    if any(strcmp(cnf_p.retracker_name,'THRESHOLD')) && reduced_set_cnf.SCOOP==0
        cnf_p.th_retracker.percentage_peak=struct.th_retracker_percentage_peak.value;
    end
    
    %% ----------------- OCOG RETRACKER ----------------------------
    if any(strcmp(cnf_p.retracker_name,'OCOG')) && reduced_set_cnf.SCOOP==0
        cnf_p.OCOG_retracker.percentage_pow_OCOG=struct.OCOG_retracker_percentage_pow_OCOG.value;
        cnf_p.OCOG_retracker.n1=struct.OCOG_retracker_n1.value; %first zero-padded sample to be used in leading edge estimation
        cnf_p.OCOG_retracker.n2=struct.OCOG_retracker_n2.value; %last zero-padded sample to be used in leading edge estimation
        cnf_p.OCOG_retracker.offset=struct.OCOG_retracker_offset.value; %additional offset
        cnf_p.OCOG_retracker.param_comp_method=str2double(struct.OCOG_retracker_param_comp_method.value); % 0: Using squares of power samples as per Frappart
        % 1: Using the power samples as per Wingham
        cnf_p.OCOG_retracker.implementation_method=str2double(struct.OCOG_retracker_implementation_method.value); % 0: Using threshold-based approach computed within n1 and n2
        % 1: Computing Amplitude reference
        % within limits defined by the window
        % centered at the COG
        % 2: Using epoch=offset+(COG-W/2)
        
    end
    
    
    %% ----------------- ANALYTICAL RETRACKER ---------------------------------
    % SAMOSA based retracker
    if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
        if reduced_set_cnf.SHAPE==0
            %--------------------------------------------------------------------------
            % ----------------------NOISE ESTIMATION ----------------------------------
            %--------------------------------------------------------------------------
            cnf_p.Thn_flag              =   struct.analytical_retracker_Thn_flag.value;              % 1 - To account for ThN in estimation; 0 - Otherwise
            cnf_p.Thn_estimation_method =   struct.analytical_retracker_Thn_estimation_method.value; % String:
            % -'external': using an external
            % fixed value (to be used in inland water)
            % -'fixed_window': using a set of
            % range bins in a given window to
            % estimate the noise given in
            % cnf_p.Thn_w_first and
            % cnf_p.Thn_w_width
            % -'adaptive': use an adaptive
            % window computed based on the
            % derivative of the waveform and a
            % given threshold cnf_p.threshold_noise, which is
            % accordingly increased with a
            % given maximum number of
            % iterations max_iter_noise
            
            % ------------- external noise ----------------------------------------
            cnf_p.external_Thn_value =   struct.analytical_retracker_external_Thn_value.value;              %in case of using an external source when THn_flag set to
            
            %--------------- fixed_window noise estimation ------------------------
            cnf_p.Thn_w_first       =   struct.analytical_retracker_Thn_w_first.value;%50;             % Gate number after IF Mask to start Thermal noise windowing; this is a subscript indice thus must be > 0
            cnf_p.Thn_w_width       =   struct.analytical_retracker_Thn_w_width.value;             % Thermal noise window width in range bins
            
            % --------------- adaptive noise esti. --------------------------------
            cnf_p.threshold_noise   = struct.analytical_retracker_threshold_noise.value; %threshold used to estimate the samples used in the noise window estimation
            % based on the derivative of the window
            cnf_p.factor_increase_noise_iter=struct.analytical_retracker_factor_increase_noise_iter.value; %factor to multiply the threshold of noise per iteration
            cnf_p.max_iter_noise    = struct.analytical_retracker_max_iter_noise.value; %number of iterations
            
            % ---------------- ML or SL estimation --------------------------------
            if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
                cnf_p.Thn_ML_SL_method        =  struct.analytical_retracker_Thn_ML_SL_method.value; % String:
                % -'ML': using the multilloked waveform
                % -'SL': estimate the noise per look from
                %  the real stack (require to pass the whole stack)
            end
        end
        
        
        %--------------------------------------------------------------------------
        % ------------------- MODEL CONFIGURATION ---------------------------------
        %--------------------------------------------------------------------------
        % ------------------ Surface parameters -----------------------------------
        
        cnf_p.rou           =   struct.analytical_retracker_rou.value; % Roughness or MSS when fitting the SWH or HS: Gaussian PDF show a roughness value in the order of 10-2
        % For Ocean is 1e-2 and for inland
        % waters 1e-4 (to be tested)
        if reduced_set_cnf.SCOOP==0
            cnf_p.Hs            =   struct.analytical_retracker_Hs.value; % value of the Hs when fitting the MSS or roughness
        end
        
        
        %------------------- Power waveform model ---------------------------------
        cnf_p.power_wfm_model=struct.analytical_retracker_power_wfm_model.value; %Define the model approximation of power wfm whether to compute
        % 'simple': Pkl=Bkl*sqrt(gl)*func_f0
        % 'complete':
        % Pkl=Bkl*sqrt(gl)*(func_f0+Tkl*gl*simga_s^2*func_f1);
        
        if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
            %-------------------- Multilooking option ---------------------------------
            cnf_p.multilook_option=struct.analytical_retracker_multilook_option.value; %String:
            % -Cris_old: thermal noise is added to
            % each SL
            % -NormML_plus_noise: ML is normalized
            % and the thermal noise is added
            % outside
            
            
            %-------------------- Masking construction --------------------------------
            cnf_p.Doppler_mask_cons_option=struct.analytical_retracker_Doppler_mask_cons_option.value; %-'external': loaded from L1B product
            %-'internal': computed in
            %the construction of the
            %model (two options)
            cnf_p.Doppler_mask_cons_internal=struct.analytical_retracker_Doppler_mask_cons_internal.value; % -'l1b_like': as done for L1B "good samples" approach
            % -'epoch_like': epoch
            % fixes the minimum point
            % of closest approach and
            % the range history is
            % computed on top of it (to be implemented)
            
            
            %--------------------------------------------------------------------------
            %----------------------- LOOK INDEXATION METHOD ---------------------------
            %--------------------------------------------------------------------------
            cnf_p.looks_index_method=   struct.analytical_retracker_looks_index_method.value;   % String:
            % - Cris_old (old version),
            % - Norm_index (norm. #pulses),
            % - Doppler_freq (exploiting beam angle, velocity and PRI info)
            % - Look_angle angle method (for flat surface: equivalent to pitch + pointing angle)
            cnf_p.look_ang_method=struct.analytical_retracker_look_ang_method.value; % String:
            % - approximate: exploiting the min and max
            % values of the Dopp., and Look angles and
            % velocity satellite's over the surface with the given PRI information per stack.
            % - exact : using the beam angles,
            % satellite's velocity, PRI
            % information for each beam in stack
            cnf_p.fd_method=struct.analytical_retracker_fd_method.value; % String:
            % - approximate: exploiting the min and max
            % values of the Dopp., and Look angles and
            % velocity satellite's over the surface with the given PRI information per stack.
            % - exact : using the beam angles,
            % satellite's velocity, PRI
            % information for each beam in stack
            
            %------------------ Range indexation ----------------------------------
            cnf_p.range_index_method=struct.analytical_retracker_range_index_method.value; % 'conventional': when creating the range vector used in fitting use directly 1:1:N_samples (even with ZP)
            % 'resampling': creating the
            % range vector as a resampling
            % of the original one (if there is a resampling)
        end
        
        
        % -------------------------------------------------------------------------
        % ---------------------- LOOK UP TABLES CONFIGURATION ---------------------
        % -------------------------------------------------------------------------
        % LUT tables in .mat should be in the processor folder
        cnf_p.lut_flag          =   struct.analytical_retracker_lut_flag.value;              % 0 or 1;if null do not use Look up tables
        if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
            cnf_p.LUT_ximax         =   struct.analytical_retracker_LUT_ximax.value; % maximum value of the gl*k parameter LUT
            cnf_p.LUT_ximin         =   struct.analytical_retracker_LUT_ximin.value; % minimum value of the gl*k parameter LUT
            cnf_p.LUT_step          =   struct.analytical_retracker_LUT_step.value; % step of gl*k parameter in the LUT
        end
        
        %--------------------------------------------------------------------------
        %----------------- PRE-PROCESSING STAGE------------------------------------
        %--------------------------------------------------------------------------
        %provide better initial estimate of epoch--
        %include a pre_processing stage to compute the initial epoch if desired
        %based on a mid-leading edge estimation as a percentatge of the peak power
        %in the actual waveform and eventually include the estimation of the noise
        %floor adjusting automatically
        cnf_p.pre_processing=struct.analytical_retracker_pre_processing.value;
        cnf_p.percent_leading_edge=struct.analytical_retracker_pre_processing_percent_leading_edge.value; % percentage of peak detect to establish the mid-point leading edge
        
        %------------------------------------------------------------------
        % ----------- Sliding window for normalization --------------------
        %based on SAMOSA v_2_5 
        if isfield(struct,'analytical_retracker_wf_Norm_AW')
            cnf_p.wf_Norm_AW         =   struct.analytical_retracker_wf_Norm_AW.value;
        end
        
        
        % -------------------------------------------------------------------------
        % ------------- FITTING CONFIGURATION -------------------------------------
        %---- Initial values-------------------------------------------------------
        cnf_p.ini_Epoch         =   struct.analytical_retracker_ini_Epoch.value;             % 45 for BB, 145 for MEAS18
        cnf_p.ini_Hs            =   struct.analytical_retracker_ini_Hs.value;
        cnf_p.ini_Pu            =   struct.analytical_retracker_ini_Pu.value;
        if reduced_set_cnf.SCOOP==0
            cnf_p.init_rou          =   struct.analytical_retracker_ini_MSS.value;
        end
        
        if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
            %------------- Fitting feedback -------------------------------------------
            cnf_p.initial_param_fit_feedback_flag = struct.analytical_retracker_initial_param_fit_feedback_flag.value; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
            cnf_p.ini_Hs_rou_sliding_win_opt = struct.analytical_retracker_ini_Hs_rou_sliding_win_opt.value;  % activate the averaging of Hs/roughness for next initial seed
            cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold = struct.analytical_retracker_ini_Hs_rou_sliding_win_opt_discard_std_thr.value; %number of times std for removal (it is also used for comparison with the accumulated one)
            cnf_p.ini_Hs_rou_sliding_win_size = struct.analytical_retracker_ini_Hs_rou_sliding_win_size.value;
        end
        %------------------ Fitting function --------------------------------------
        if reduced_set_cnf.SHAPE==0
            cnf_p.fitting_fun_type=struct.fitting_fun_type.value; % Define the fitting procedure used
            
            cnf_p.lsq_algorithm      = struct.lsq_algorithm.value;
            cnf_p.fitting_options_lb = struct.fitting_options_lb.value; %lower boundaries
            if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
                cnf_p.fitting_options_ub = struct.fitting_options_ub.value; %upper boundaries
            end
            %algorithm
            switch cnf_p.fitting_fun_type
                case 'lsq'
                    %switch cnf_p.lsq_algorithm
                    %case 'trust-region-reflective'
                    cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'Display','off'); %'UseParallel','always' %'TolFun',1e-8,'TolX',1e-8
                    %case 'levenberg-marquardt'
                    %    cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'Display','off'); %'TolFun',1e-8,'TolX',1e-8
                    %end
                    %cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'TolFun',1e-8,'TolX',1e-8,'Display','off','UseParallel','always');
                case 'fmin'
                    cnf_p.fitting_options     =   optimset('Display','off','TolFun',1e-7,'TolX',1e-7);
            end
        end
        %------------------ Two-step fitting ----------------------------------
        if reduced_set_cnf.SCOOP==0
            cnf_p.two_step_fitting_COR_threshold_rou=struct.analytical_retracker_two_step_fitting_COR_threshold_rou.value; %threshold below which second fitting on roughness is used
        end
        
    end
    
    %% -------------------- MISFIT QUALITY  -------------------------------
    %used for Sea State CCI to set quality of the retrieval:
    %above the MISfit value a bad retrieval is considered
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        cnf_p.quality_flag_misfit_th=struct.quality_flag_misfit_th.value;
    end
    
    %% -------------------- GEOPHYSICAL CORRECTIONS ----------------------------
    % -------------------------------------------------------------------------
    % Corrections on SSH
    cnf_p.geo_corr_application_flag=struct.geo_corr_application_flag.value;
    if reduced_set_cnf.SCOOP==0
        cnf_p.force_geocorr_surf_type=struct.force_geocorr_surf_type.value; % 0: use the geophysical corrections depending on the flag surface type included in the L1B product
        % 1: use the same group of geophysical corrections to all the records independently of the flag of the surface type: forcing the type surface
        cnf_p.product_type_surface=struct.product_type_surface.value; % open_ocean, sea_ice, land_ice ...
    end
    %----------------------------------------------------------------------
    % Atmospheric attenuation correction
    if reduced_set_cnf.SHAPE==0
        cnf_p.atm_att_correction_flag = struct.atm_att_correction_flag.value;
    end
    
    %% -------------------- OUTPUT PRODUCT -------------------------------------
    % -------------------------------------------------------------------------
    %Inherete from Sentinel-3 file formating
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        cnf_p.Product_type=struct.Product_type.value; % String: 'LAN' for land products      
        %'WAT' for water products
    end
    
    cnf_p.write_output=struct.write_output.value;
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        cnf_p.output_product_format=struct.output_product_format.value; % nc: NetCDF, mat: Matlab
        cnf_p.nc_format_type = struct.nc_format_type.value; %use default or SS CCI
        cnf_p.nc_name_surface_height=struct.nc_name_surface_height.value; %string indicating which name to be used in the netcdf for the surface height: SSH (sea surface height) or SH (surface height)
    end
    
    if reduced_set_cnf.SCOOP==0 && reduced_set_cnf.SHAPE==0
        cnf_p.optional_ext_file_flag =struct.optional_ext_file_flag.value;
        %cnf_p.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
        cnf_p.file_ext_string=struct.file_ext_string.value;
    end
    
    
    %% -------------------- PLOTING OPTIONS -----------------------------------
    cnf_p.plot_fits_flag    =   struct.plot_fits_flag.value;              % 0 - not plot; 1 - plot fit results
    cnf_p.plot_fits_lat_range = struct.plot_fits_lat_range.value;%[-91,91]; %latitude range of waveforms fit to be plotted (by default between -91 and 91 to force all of them to be plotted)
    cnf_p.plot_fits_downsampling =struct.plot_fits_downsampling.value;
    cnf_p.visible_figures=struct.visible_figures.value;
    if ~exist('struct.text_interpreter')
        struct.text_interpreter = 'latex';
    end
    cnf_p.text_interpreter=struct.text_interpreter;

    clear struct;

end
