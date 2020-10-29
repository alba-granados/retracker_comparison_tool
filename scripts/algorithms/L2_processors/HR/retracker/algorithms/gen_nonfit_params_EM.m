function nf_p = gen_nonfit_params_EM (data,cnf_p,chd_p,cst_p) 
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This function initializates the non fitting parameters' structure
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT v1 16/06/2016
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       data    =   input data structure for the L2 processor
%       cnf_p           =   configuration parameters
% OUTPUT:
%       ini_p           =   Fitting parameters seed

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
% Versions control:
% v1.0:

%% ------------------- Height statistical distribution --------------------
%if ~cnf_p.rou_flag
    nf_p.waveskew   =   0;
    nf_p.EMbias     =   0;
    nf_p.rou    =   cnf_p.rou; % from Yuguan Liu
%end

%% ------------------ Array values ----------------------------------------
%% ------------------ Geometry/Projections --------------------------------
% Height
nf_p.h       =   data.GEO.H;
% Spherical earth term
R_Earth =sqrt((cst_p.semi_major_axis_cst^4*cos(data.GEO.LAT.*pi/180).^2+cst_p.semi_minor_axis_cst^4*sin(data.GEO.LAT*pi/180).^2)./((cst_p.semi_major_axis_cst*cos(data.GEO.LAT.*pi/180)).^2+(cst_p.semi_minor_axis_cst*sin(data.GEO.LAT.*pi/180)).^2));
nf_p.alpha        =   1 + (nf_p.h)./(R_Earth);

% Pitch and roll angles projection as distances
% assuming no change of sign is performed in the lecture
%take into account that the reference system for Model (x,y,z) is having z
%normal to the surface and x in the plane formed by satellite velocity and
%z : positive y component on the left_side; x positive forward movement
%satellite
switch cnf_p.mission
    case {'CS2','CR2'}
        % Satellite-frame coordinates: x-along-track, z-opposite direction
        % to nadir and positive y on left-side of satellite
        %pitch is defined as positive/negative nose down/up respectively        
        nf_p.xp      =   data.GEO.H.*(-1.0*data.GEO.pitch);
        %nf_p.xp      =   data.GEO.H.*(cnf_p.sign_pitch*1.0*data.GEO.pitch);
        %roll-angle is defined as positive/negative antenna right down/up        
        nf_p.yp      =   data.GEO.H.*(1.0*data.GEO.roll);
        %nf_p.yp      =   data.GEO.H.*(cnf_p.sign_roll*data.GEO.roll);
    case 'S3' % alva
        %TBD
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   data.GEO.H.*(1.0*data.GEO.pitch);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   data.GEO.H.*(1.0*data.GEO.roll);
    case 'S6'
        % Satellite-frame coordinates: x-along-track, z direction
        % to nadir and positive y on right-side of satellite
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   data.GEO.H.*(1.0*data.GEO.pitch);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   data.GEO.H.*(data.GEO.roll);        
end


%Equivalent antenna beamwidth projection on ground
nf_p.alphax  =   (8*log(2))./(data.GEO.H.*chd_p.antenna_beamwidth_alt_ku_chd).^2;
nf_p.alphay  =   (8*log(2))./(data.GEO.H.*chd_p.antenna_beamwidth_act_ku_chd).^2; %8 comes from fact that using the two way antenna pattern

%% --------------- Sampling spacing ---------------------------------------
% along-track
nf_p.Lx      =   cst_p.c_cst*data.GEO.H./(2*data.GEO.V.*chd_p.freq_ku_chd*chd_p.N_total_pulses_b_chd.*data.HRM.pri_surf);

% across-track & vertical 
switch cnf_p.range_index_method
        case 'conventional'
            nf_p.Ly      =   sqrt(cst_p.c_cst*data.GEO.H./(nf_p.alpha.*data.HRM.fs_clock_ku_surf*cnf_p.ZP));
            nf_p.Lz      =   0.5*cst_p.c_cst./(data.HRM.fs_clock_ku_surf*cnf_p.ZP);
            
        case 'resampling'           
%             nf_p.Ly      =   sqrt(cst_p.c_cst*data.GEO.H./(nf_p.alpha*chd_p.bw_rx_ku_chd));
%             nf_p.Lz      =   0.5*cst_p.c_cst/(chd_p.bw_rx_ku_chd);
            nf_p.Ly      =   sqrt(cst_p.c_cst*data.GEO.H./(nf_p.alpha.*data.HRM.fs_clock_ku_surf));
            nf_p.Lz      =   0.5*cst_p.c_cst./(data.HRM.fs_clock_ku_surf);
end

% not clear from where this term is coming
nf_p.Lgamma  =   nf_p.alpha./(2*data.GEO.H.*nf_p.alphay);

%% --------------- Number of looks that form the stack --------------------
% it would be probably more convenient to define a unique structure value of Neff
% no need to use data.SIN.Neff (will be kept as originally proposed just in case there are some future elaborations on the SARin data)
switch  cnf_p.mode
%     case 'SARin'
%         nf_p.Neff    =   data.SIN.Neff;
    case {'SAR','SARin','RAW','HR','RMC','LR-RMC','FF-RAW','FF-RMC'}        
        nf_p.Neff    =   data.HRM.Neff;
    otherwise
        error('No valid operational mode')
end

%% ------------------- Constant Values ------------------------------------
nf_p.Npulses    =   chd_p.N_total_pulses_b_chd;
nf_p.alphag_a   =   chd_p.alpha_ga_chd;
nf_p.alphag_r   =   chd_p.alpha_gr_chd;
nf_p.A_s2Ga     =   chd_p.A_s2Ga_chd;
nf_p.A_s2Gr     =   chd_p.A_s2Gr_chd;
nf_p.Nbcycle    =   chd_p.N_bursts_cycle_chd;
nf_p.bw_Rx      =   chd_p.bw_rx_ku_chd;
nf_p.Nsamples   =   chd_p.N_samples_chd;


end