%to avoid redefining cnf file and all the implemented functions
%--------------------------------------------------------------------------
%----------------------- SEED INFORMATION ---------------------------------
%--------------------------------------------------------------------------
cnf_p.seed              =   0;              % 0 - no seed introduced; 1 - seed provided by phase information of SARin data


%--------------------------------------------------------------------------
% -------------------- IF MASK --------------------------------------------
%--------------------------------------------------------------------------
cnf_p.IFmask_N          =   0;              % First/last N samples of the waveform affected by the IF Mask filter which should be excluded from processing. Specify if known, 0 if not known


%--------------------------------------------------------------------------
% ----------------------NOISE ESTIMATION ----------------------------------
%--------------------------------------------------------------------------
cnf_p.Thn_flag              =   1;              % 1 - To account for ThN in estimation; 0 - Otherwise
cnf_p.Thn_estimation_method =   'external'; % String:
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
cnf_p.external_Thn_value =   0;              %in case of using an external source when THn_flag set to

%--------------- fixed_window noise estimation ------------------------
cnf_p.Thn_w_first       =   1;%50;             % Gate number after IF Mask to start Thermal noise windowing; this is a subscript indice thus must be > 0
cnf_p.Thn_w_width       =   20;%20;             % Thermal noise window width in range bins

% --------------- adaptive noise esti. --------------------------------
cnf_p.threshold_noise   = 1e-3; %threshold used to estimate the samples used in the noise window estimation
% based on the derivative of the window
cnf_p.factor_increase_noise_iter=1.5; %factor to multiply the threshold of noise per iteration
cnf_p.max_iter_noise    = 100; %number of iterations

% ---------------- ML or SL estimation --------------------------------
cnf_p.Thn_ML_SL_method        =  'ML'; % String:
% -'ML': using the multilloked waveform
% -'SL': estimate the noise per look from
%  the real stack (require to pass the whole stack)



%--------------------------------------------------------------------------
%----------------------- LOOK INDEXATION METHOD ---------------------------
%--------------------------------------------------------------------------
cnf_p.looks_index_method=   'Look_angle';   % String:
% - Cris_old (old version),
% - Norm_index (norm. #pulses),
% - Doppler_freq (exploiting beam angle, velocity and PRI info)
% - Look_angle angle method (for flat surface: equivalent to pitch + pointing angle)
switch cnf_p.looks_index_method
    case 'Look_angle'
        cnf_p.look_ang_method='approximate'; % String:
        % - approximate: exploiting the min and max
        % values of the Dopp., and Look angles and
        % velocity satellite's over the surface with the given PRI information per stack.
        % - exact : using the beam angles,
        % satellite's velocity, PRI
        % information for each beam in stack
    case 'Doppler_freq'
        cnf_p.fd_method='approximate'; % String:
        % - approximate: exploiting the min and max
        % values of the Dopp., and Look angles and
        % velocity satellite's over the surface with the given PRI information per stack.
        % - exact : using the beam angles,
        % satellite's velocity, PRI
        % information for each beam in stack
end

%------------------- Model parameters -------------------------------------
%-------------------- Multilooking option ---------------------------------
cnf_p.multilook_option='Cris_old'; %String:
% -Cris_old: thermal noise is added to
% each SL
% -NormML_plus_noise: ML is normalized
% and the thermal noise is added
% outside

%-------------------- Masking construction --------------------------------
cnf_p.Doppler_mask_cons_option='external'; %-'external': loaded from L1B product
%-'internal': computed in
%the construction of the
%model (two options)
cnf_p.Doppler_mask_cons_internal='l1b_like'; % -'l1b_like': as done for L1B "good samples" approach
% -'epoch_like': epoch
% fixes the minimum point
% of closest approach and
% the range history is
% computed on top of it (to be implemented)

%------------------ LUT cnf -----------------------------------------------
cnf_p.LUT_ximax         =   50; % maximum value of the gl*k parameter LUT
cnf_p.LUT_ximin         =   -10; % minimum value of the gl*k parameter LUT
cnf_p.LUT_step          =   1e-5; % step of gl*k parameter in the LUT

%------------- Fitting feedback -------------------------------------------
cnf_p.initial_param_fit_feedback_flag = 0; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
cnf_p.ini_Hs_rou_sliding_win_opt = 0;  % activate the averaging of Hs/roughness for next initial seed
cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold = 3; %number of times std for removal (it is also used for comparison with the accumulated one)
cnf_p.ini_Hs_rou_sliding_win_size = 50;

%------------------ Range indexation ----------------------------------
cnf_p.range_index_method='conventional'; % 'conventional': when creating the range vector used in fitting use directly 1:1:N_samples (even with ZP)
% 'resampling': creating the
% range vector as a resampling
% of the original one (if there is a resampling)


cnf_p.optional_ext_file_flag =0;
%cnf_p.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
cnf_p.file_ext_string='adaptive_noise_whole_win';


%------------------ Fitting function --------------------------------------
cnf_p.fitting_fun_type='lsq'; % Define the fitting procedure used

%String:
%-'lsq': Using the lsqcurvfit based on levegender-Marquard
%-'fmin': Using the fminsearch
%algorithm
switch cnf_p.fitting_fun_type
    case 'lsq'
        %trust-region-reflective
        %levenberg-marquardt
        %'MaxIter',400
        cnf_p.lsq_algorithm      = 'trust-region-reflective';
        cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'Display','off','UseParallel','always'); %'TolFun',1e-8,'TolX',1e-8
        %cnf_p.fitting_options    = optimset('Algorithm',cnf_p.lsq_algorithm,'TolFun',1e-8,'TolX',1e-8,'Display','off','UseParallel','always');
        cnf_p.fitting_options_lb = [0,0,0]; %lower boundaries
        cnf_p.fitting_options_ub = []; %upper boundaries
    case 'fmin'
        cnf_p.fitting_options     =   optimset('Display','off','TolFun',1e-7,'TolX',1e-7);
        cnf_p.fitting_options_lb = []; %lower boundaries
        cnf_p.fitting_options_ub = []; %upper boundaries
end
