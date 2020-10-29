%to avoid redefining cnf file and all the implemented functions
%--------------------------------------------------------------------------
%----------------------- MISSION AND MODE ---------------------------------
%--------------------------------------------------------------------------
cnf_p.mission           =   'CS2';  % {'CS2','S3','S6'}
cnf_p.mode              =   'SAR';


%--------------------------------------------------------------------------
%----------------------- L1B PROCESSOR ------------------------------------
%--------------------------------------------------------------------------
cnf_p.L1proc            =   'ISD';

%--------------------------------------------------------------------------
%----------------------- SEED INFORMATION ---------------------------------
%--------------------------------------------------------------------------
cnf_p.seed              =   0;              % 0 - no seed introduced; 1 - seed provided by phase information of SARin data


%--------------------------------------------------------------------------
% -------------------- IF MASK --------------------------------------------
%--------------------------------------------------------------------------
cnf_p.IFmask_N          =   0;              % First/last N samples of the waveform affected by the IF Mask filter which should be excluded from processing. Specify if known, 0 if not known


%--------------------------------------------------------------------------
%--------------------- WAVEFORM PORTION SELECTION -------------------------
%--------------------------------------------------------------------------
%activate the filtering to select the part of waveform of interest
cnf_p.wvfm_portion_selec=0;
cnf_p.wvfm_portion_selec_type=''; % -'ref_height': using a reference height from DEM
cnf_p.wvfm_portion_selec_DEM_ref='';
cnf_p.wvfm_portion_selec_DEM_path=''; %path to the DEM to be used or saved in case of SRTM
cnf_p.wvfm_portion_selec_l_samples=[]; % left samples w.r.t to position of ref. position within window (either peak or valley)
cnf_p.wvfm_portion_selec_r_samples=[]; % right samples w.r.t to position of ref. position within window (either peak or valley)
cnf_p.peak_prominence_norm=[]; % define the minimum prominence of the peaks to be sorted out in a multipeak scenario


%% THRESHOLD AND OCOG RETRACKERS
%--------------------------------------------------------------------------
%--------------------- THRESHOLD RETRACKER --------------------------------
%--------------------------------------------------------------------------
cnf_p.th_retracker.percentage_peak=[];

%--------------------------------------------------------------------------
%--------------------- OCOG RETRACKER --------------------------------
%--------------------------------------------------------------------------
cnf_p.OCOG_retracker.percentage_pow_OCOG=[];
cnf_p.OCOG_retracker.n1=[]; %first zero-padded sample to be used in leading edge estimation
cnf_p.OCOG_retracker.n2=[]; %last zero-padded sample to be used in leading edge estimation
cnf_p.OCOG_retracker.offset=[]; %additional offset
cnf_p.OCOG_retracker.param_comp_method=[];
% 1: Using the power samples as per Wingham
cnf_p.OCOG_retracker.implementation_method=[];


%% ---------------------- ANALYTICAL RETRACKER OPTIONS ----------------------
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
%---------------------- NOISE ESTIMATION ML O SL --------------------------
% -------------------------------------------------------------------------
cnf_p.Thn_ML_SL_method        =  'ML'; % String:
% -'ML': using the multilloked waveform
% -'SL': estimate the noise per look from
%  the real stack (require to pass the whole stack)

%--------------------------------------------------------------------------
%---------------------- SWH VALUE WHEN FITTING MSS ------------------------
% -------------------------------------------------------------------------
cnf_p.Hs = [];

%--------------------------------------------------------------------------
%-------------------- Multilooking option ---------------------------------
%--------------------------------------------------------------------------
cnf_p.multilook_option='Cris_old'; 

%--------------------------------------------------------------------------
%-------------------- Masking construction --------------------------------
%--------------------------------------------------------------------------
cnf_p.Doppler_mask_cons_option='external';
cnf_p.Doppler_mask_cons_internal=[];


%--------------------------------------------------------------------------
%----------------------- LOOK INDEXATION METHOD ---------------------------
%--------------------------------------------------------------------------
cnf_p.looks_index_method=   'Look_angle';   
cnf_p.look_ang_method='approximate';
cnf_p.fd_method=[];

%------------------ Range indexation ----------------------------------
cnf_p.range_index_method='conventional';



%------------------ LUT cnf -----------------------------------------------
cnf_p.LUT_ximax         =   50; % maximum value of the gl*k parameter LUT
cnf_p.LUT_ximin         =   -10; % minimum value of the gl*k parameter LUT
cnf_p.LUT_step          =   1e-5; % step of gl*k parameter in the LUT

% -------------------FITTING OPTIONS --------------------------------------
%-------------- Initial values---------------------------------------------
cnf_p.init_rou = [];

%------------- Fitting feedback -------------------------------------------
cnf_p.initial_param_fit_feedback_flag = 1; %activate the updating of the initial parameters per fitted waveform either using the sliding window or with acumulation of previous info
cnf_p.ini_Hs_rou_sliding_win_opt = 1;  % activate the averaging of Hs/roughness for next initial seed
cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold = 3; %number of times std for removal (it is also used for comparison with the accumulated one)
cnf_p.ini_Hs_rou_sliding_win_size = 10;

% ---------------- Upper bounds on fitting --------------------------------
cnf_p.fitting_options_ub =[];

%------------------ Two-step fitting threshold ----------------------------
cnf_p.two_step_fitting_COR_threshold_rou=[];


cnf_p.optional_ext_file_flag =0;
%cnf_p.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
cnf_p.file_ext_string='adaptive_noise_whole_win';

%--------------------------------------------------------------------------
%-------------------------- MISFIT QUALITY --------------------------------
%--------------------------------------------------------------------------
cnf_p.quality_flag_misfit_th=100; % so no filtering applies
%only used in SS-CCI


%--------------------------------------------------------------------------
%------------------------GEOPHYSICAL CORRECTIONS --------------------------
%--------------------------------------------------------------------------
cnf_p.force_geocorr_surf_type=1;
cnf_p.product_type_surface='open_ocean';

%--------------------------------------------------------------------------
%------------------------OUTPUT PRODUCT -----------------------------------
%--------------------------------------------------------------------------
cnf_p.Product_type='WAT';
cnf_p.output_product_format='nc';
cnf_p.nc_format_type = 'default';
cnf_p.nc_name_surface_height='ssh';

cnf_p.optional_ext_file_flag =0;
%cnf_p.file_ext_string='PTR_sigmaa_new_sigmar_cris_Lz_Ly_ZP_rg_idx_conv';
cnf_p.file_ext_string='';




