function fit_res = fitting_noise_method(data, cnf_p, chd_p , nf_p, ini_p,varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to perform the fitting of the waveforms
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%                   
%                   
%
% Reviewer:        M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul /
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% INPUT:
%   MANDATORY:
%       data        =   SAR data structure 
%       nf_p        =   Non-fit parameters of the fitting
%       cnf_p       =   Configuration parameters of the L2 processor
%       ini_p       =   initialization parameters of need for the fitting
%   OPTIONAL:
%       LUT_f0_file =   file containing the look up table for the f0
%                       function parametrization (full path)
%       LUT_f1_file =   file containing the look up table for the f1
%                       function parametrization (full path)
%       path_Results=   path where to save the fitting plot results
%       L1B_filename=   Filename of the reference L1B processed file to be
%                       used as extension for the plot file name
%       
%
% OUTPUT:
%       fit_res     =   structure containing fitting results
%           .Epoch  =   Epoch not compensated for IF_Mask
%           .Hs     =   significant wave height
%           .Pu     =   Pu if to be fitted
%           .GOF    =   Goodness of fitting in %
%           .flag   =   fitting flag output = the same as for lsqcurve fit
%                       flags
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - ml_wave_gen_EM: genrates the multilook waveform based on Chris Ray
% analytical model in IEEE TGRS "SAR Altimeter Backscattered Waveform Model" 
% DOI:10.1109/TGRS.2014.23330423
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - The routine will not work for option of fixing the size of the window
% to estimate the noise (not estimating it from the waveform) when fitting
% only a portion of the waveform to be fitted
% - the defintiion of the two structures nf_p and nf_params is inhereted
% from the code of cristina and should be probably refined in the
% definition
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% V7.1: Upgrading the look indexation using different methods based on
% normalized indexation (first rough approach), Doppler freq info, look
% angle information. The latter two methods have the exact and approximate approaches
% with the possibility to use the whole angular information (beam angles
% and look angles respectively) from the stack itself or using the start
% and stop angles. Vector of look indexation is passed to the ml waveform
% function.
% V7.2: Option to perform the estimation of the noise per look provided the
% whole stack is accessible to the L2 processing (to date such option will be 
% considered as no big differences have been observed)
% V7.3: It has been incorportaed the option to construct the Power waveform
% model as simple (using only f0 function: equation (43) in paper Chris
% Ray) or the complete option (including f0 and f1, as in equation (42)).
% It has been considered the option to implement sucgh f0 and f1 using
% their expansions in termns of bessel functions or using Look Up Tables to
% accelerate the fitting.
% V7.4: It has been included the option to perform the LSE minimization
% either using the lsqcurvefit or the fminsearch (minimizing the MSE)
% V7.5: A pre-processing optional stage has been incorporated to perform a
% rough threshold-based retracker to have a first initial guess of the epoch
%, which can be used instead of the one provided by the user 
%(some issues encountered in the fitting procedure depending on 
% the initial epoch guess provided by the user --> local minimums). In such
% pre-processing is it also possible to perform an adaptive estimation of
% the noise thermal noise based on the derivative of the real waveform
% without using a fixed set of initial sample and window size for the
% estimation of the noise floor.
% V7.6: Construction of the stack using matrix notation instead of for
% loops to accelerate the formation of the multilook waveform. Computation
% and generation of the stack mask based on the stack_mask vector provided
% as input in the L1B.
% V7.7: Provide an average estimation of the initial parameters for the
% next surface based on a sliding window or as an accumulated average value
% from the first surface for both the Hs and roughness. Specific filtering
% in the average computation is made based on the mean and std.
% v7.8: Incorporating the option to perform a two-step fitting first on a
% ocean-like waveform using SWH and the second depending on the correlation
% threshold based on roughness fitting leads-like surface. Keep the best
% fitting among the two of them
% v7.9: Changes on the computation of the Beam angles for the approximate
% case using the Doppler angle and the look angle considering
% linspace(pi/2+Dopp_start-look_start,pi/2+Dopp_stop-look_stop);
% conditioning of the statics for feedback initialization of the SWH and
% roughness in case of NaN values present
% v7.10: Include the possibility to perform a pre-processing to select the
% portion of the waveform to be fitted.
% v7.11: 29.03.2017 For noise inclusion three different options with a single parameter
% has been considered to avoid misunderstanding with the previous version
% using different configuration flags for external, window and adaptive
% v7.12: 21/11/2017: Possibility to discard samples at the beginning and
% end of the waveform



%% ---------------- Handling input variables ------------------------------
if(nargin<5 || nargin>(5+4*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('LUT_f0_file',{''},@(x)ischar(x));
p.addParamValue('LUT_f1_file',{''},@(x)ischar(x));
p.addParamValue('path_Results',{''},@(x)ischar(x));
p.addParamValue('L1B_filename',{''},@(x)ischar(x));
p.parse(varargin{:});
LUT_f0_file=char(p.Results.LUT_f0_file);
LUT_f1_file=char(p.Results.LUT_f1_file);
path_Results=char(p.Results.path_Results);
L1B_filename=char(p.Results.L1B_filename);
clear p;

%file name identification
switch cnf_p.mission
    case {'CS2','CR2'}
        % --------------------- CroySAT-2 ---------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                file_id=L1B_filename(20:50);
            case 'GPOD'
                file_id=L1B_filename(24:54);
            case 'ISD'
                file_id=L1B_filename(17:47);
        end
    case {'S6','JCS'}
        % --------------------- Sentinel-6 -----------------------------
        file_id=L1B_filename;
    case {'S3A','S3B','S3'}
        switch cnf_p.L1proc
            case 'ESA'
                file_id=L1B_filename(17:47);
            case 'ISD'
                file_id=L1B_filename(17:47);
            case 'DeDop'
                file_id=L1B_filename;
        end
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end


%% -------------- LOADING THE LUTS -----------------------------------------
%--------------------------------------------------------------------------
% added by EM 01.03.2016 include the LUT load
% modified by EM 15.03.2016 for func_f1
if (cnf_p.lut_flag)
    [~,~,LUT_f0_file_ext]=fileparts(LUT_f0_file);
    switch LUT_f0_file_ext
        case '.mat'
            load(LUT_f0_file,'func_f0');
        otherwise
            error('No valid LUT f0 file format');
    end
    switch cnf_p.power_wfm_model
        case 'complete'
            [~,~,LUT_f1_file_ext]=fileparts(LUT_f1_file);
            switch LUT_f1_file_ext
                case '.mat'
                    load(LUT_f1_file,'func_f1');
                otherwise
                    error('No valid LUT f1 file format');
            end
        case 'simple'    
            %
        otherwise
            error('No valid power waveform model')
    end
end



%% ---------- NON-FIT MODEL PARAMETERS SINGLE VARIABLE --------------------
% read non-fit parameters which are single value
% -------------------------------------------------------------------------
    nf_params.waveskew  =   nf_p.waveskew;
    nf_params.EMbias    =   nf_p.EMbias;
    nf_params.rou       =   nf_p.rou;
    nf_params.Npulses   =   nf_p.Npulses;
    nf_params.alphag_a  =   nf_p.alphag_a;
    nf_params.alphag_r  =   nf_p.alphag_r;
    nf_params.Nbcycle   =   nf_p.Nbcycle;
    nf_params.bw_Rx     =   nf_p.bw_Rx;
    nf_params.Nsamples  =   nf_p.Nsamples;
    
    nf_params.A_s2Ga= nf_p.A_s2Ga;
    nf_params.A_s2Gr= nf_p.A_s2Gr;
    
    
%% ------------------ FITTING SEEDS ---------------------------------------
% Initialize parameters to be fitting according to flags
% -------------------------------------------------------------------------
    %EM: 13.09.2016
    if cnf_p.two_step_fitting && cnf_p.rou_flag
        cnf_p.rou_flag=0; %ensure first the SWH fitting is performed
    end
    
    if cnf_p.rou_flag
        if strcmp(cnf_p.multilook_option,'Cris_old')
            fit_params_ini      =   [ini_p.Epoch ini_p.rou ini_p.Pu];
        else
            fit_params_ini      =   [ini_p.Epoch ini_p.rou];
        end
        
    else
        if strcmp(cnf_p.multilook_option,'Cris_old')
            fit_params_ini      =   [ini_p.Epoch ini_p.sigmaz ini_p.Pu];
        else
            fit_params_ini      =   [ini_p.Epoch ini_p.sigmaz];
        end
        
    end
%     %added by Em 16.03.2016 to include fitting optimization of thermal
%     %noise floor
%     if cnf_p.fit_noise
%         fit_params_ini=[fit_params_ini,ini_p.Thn_1st,ini_p.Thn_wdth];
%     end
    

% from the above it is clear that maximum three parameters are to be
% fitted. 
% (1)- Epoch --- is always fitted
% (2)- mss/Hs   --- roughtness or significant wave height are to be fitted.
% Either one or the other
% (3)- Pu   --- amplitude is to be fitted

% -------------------------------------------------------------------------
% IF MASK - nOT TO BE USED IN THIS VERSION OF THE RETRACKER
% -------------------------------------------------------------------------
switch cnf_p.mode
    case {'SAR','SARin','RAW','LR-RMC','FF-RAW','FF-RMC'}
        data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    case {'RMC'} %RMC
        data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end,:)';
%     case {'SARin'}
%         data.SIN.power_wav_filtered    =   data.SIN.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    otherwise
        error('not a valid mode');
end

%% --------------- Define output structure --------------------------------
%--------------------------------------------------------------------------
fit_res.Epoch=zeros(1,data.N_records);
fit_res.rou=zeros(1,data.N_records);
fit_res.Hs=zeros(1,data.N_records);
fit_res.Pu=zeros(1,data.N_records);                     
fit_res.flag=-1.0*ones(1,data.N_records);
fit_res.COR=zeros(1,data.N_records);
fit_res.misfit=zeros(1,data.N_records);

accumulated=0; %keep an accumulated value of the Hs for 
%parameter initialization for next surface
m_rou=0; %keep control on how many surfaces enter the second step-fitting on roughness
m_fitted=0; %keep number of fitted waveforms (avoiding those where L1B waveform is not valid)

%due to missing initialization in L1B processor for check continuty
idx_pri_non_zero=find(data.HRM.pri_surf,1);


% modelled_stacks = NaN(max(data.HRM.Neff),data.N_samples,data.N_records);
% look_angle_stacks = NaN(max(data.HRM.Neff),data.N_records);

%% ----------------- INITIALIZATION OF INDEX RANGE VECTOR ------------------
%--------------------------------------------------------------------------
switch cnf_p.range_index_method
        case 'conventional'
            range_index_comm=1:1:data.N_samples;
        case 'resampling'            
            switch cnf_p.mission
                case {'S6'}
                  %there is an on-board resampling (fs~=BW): 
                  %range indexation original shall be accordingly resampled
                  %for on-board and/or any on-ground resampling
%                   
%                   %compute the N_samples original without any resampling
%                   N_samples_orig = ceil(data.N_samples/(nf_p.oversampling*cnf_p.ZP));
%                   range_index_comm=1:1/(nf_p.oversampling*cnf_p.ZP):N_samples_orig;
                  
%                   range_index_comm=1:1/(cnf_p.ZP):chd_p.N_samples_chd;
                  
                  range_index_comm=[1:1/cnf_p.ZP:(data.N_samples/cnf_p.ZP),...
                                    (data.N_samples/cnf_p.ZP+1/cnf_p.ZP):1/cnf_p.ZP:(data.N_samples/cnf_p.ZP+1/cnf_p.ZP*(cnf_p.ZP-1))];
                  
                otherwise
%                   %range_index_comm=interp((1:data.N_samples/cnf_p.ZP),cnf_p.ZP);
                  range_index_comm=[1:1/cnf_p.ZP:(data.N_samples/cnf_p.ZP),...
                                    (data.N_samples/cnf_p.ZP+1/cnf_p.ZP):1/cnf_p.ZP:(data.N_samples/cnf_p.ZP+1/cnf_p.ZP*(cnf_p.ZP-1))];
                                
%                  range_index_comm=1:1/(cnf_p.ZP):chd_p.N_samples_chd;
            end
end

%ThN_estimate_accumulated=0.0;

% % EM 05.09.2017 Save the modelled stacks and the corresponding look angles
% %modelled_stacks = NaN(data.N_records,max(nf_p.Neff),data.N_samples);
% modelled_stacks = NaN(max(nf_p.Neff),data.N_samples,data.N_records);
% look_angle_stacks = NaN(data.N_records,max(nf_p.Neff));

%normalized_noise_estimation = NaN(1,data.N_records);

%% ---------------------- Run fitting per waveform ------------------------
% -------------------------------------------------------------------------
% Least Square Fitting
% -------------------------------------------------------------------------

disp(strcat('Total # wvfms ',num2str(data.N_records)))
    for m = 1:data.N_records
        if mod(m,10)==0
            disp(strcat('Surface #',num2str(m)));
        end
        %t_fitting_init=tic;
        %disp(m)        
%% ---------------------- WAVEFORM SELECTION PORTION ----------------------
%--------------------------------------------------------------------------            
        %added by EM: 21.03.2017
        if cnf_p.wvfm_portion_selec
            if strcmpi(cnf_p.wvfm_portion_selec_type,'ref_height')
                [start_sample,stop_sample,flag_nadir_return_in_win,epoch_ref_DEM] = wvfm_portion_selec(data.HRM.power_wav_filtered(m,:),data.MEA.win_delay(m),cnf_p,chd_p,m,...
                                                                    'wd_ref_DEM',data.GEO.wd_ref_DEM(m),'path_Results',path_Results,'L1B_filename',L1B_filename);
            elseif strcmpi(cnf_p.wvfm_portion_selec_type,'CP4O')
                [start_sample,stop_sample,flag_nadir_return_in_win,epoch_ref_DEM] = wvfm_portion_selec(data.HRM.power_wav_filtered(m,:),data.MEA.win_delay(m),cnf_p,chd_p,m,...
                                                                    'wd_ref_DEM',data.CP4O.seedpw(m),'path_Results',path_Results,'L1B_filename',L1B_filename);
            else
                [start_sample,stop_sample,flag_nadir_return_in_win,epoch_ref_DEM] = wvfm_portion_selec(data.HRM.power_wav_filtered(m,:),data.MEA.win_delay(m),cnf_p,chd_p,m,...
                                                                      'path_Results',path_Results,'L1B_filename',L1B_filename);
            end
        else
            start_sample=1; %first sample to be used
            stop_sample=data.N_samples; %last sample to be used
            flag_nadir_return_in_win=NaN;
        end
       
        
%% ----------------- Discarding samples beginning and end -----------------
if cnf_p.wvfm_discard_samples
    start_sample=max(start_sample,1+cnf_p.wvfm_discard_samples_begin);
    stop_sample=min(stop_sample,data.N_samples-cnf_p.wvfm_discard_samples_end);    
end

range_index=range_index_comm(start_sample:stop_sample);

        
%% ------------------- Checking validity of L1B waveform ------------------
%--------------------------------------------------------------------------
        if ~any(~isnan(data.HRM.power_wav_filtered(m,start_sample:stop_sample))) || ~any(data.HRM.power_wav_filtered(m,start_sample:stop_sample))
            fit_res.Epoch(m)    =   NaN;
            if cnf_p.rou_flag
                fit_res.rou(m)  =   NaN; 
            else
                fit_res.Hs(m)   =   NaN;
                if cnf_p.two_step_fitting
                    fit_res.rou(m)  =   NaN;
                end
            end
            fit_res.Pu(m)       =   NaN;
            fit_res.flag(m)     =   NaN;           
            fit_res.COR(m)      =   NaN; 
            fit_res.misfit(m)   =   NaN;
            fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
            fit_res.Flag_quality(m) = 0; % waveform with no information
            continue;
        end
        m_fitted=m_fitted+1;

%             fprintf('-------\n');
%             fprintf('- Fitting wave num = %d\n',m-1);

%% --------- NON-FIT MODEL PARAMETERS DEFINED AS ARRAY
% read non-fit parameters which are single value
% -------------------------------------------------------------------------
            nf_params.h         =   nf_p.h(m);
            nf_params.xp        =   nf_p.xp(m);
            nf_params.yp        =   nf_p.yp(m);
            nf_params.alphax    =   nf_p.alphax(m);
            nf_params.alphay 	=   nf_p.alphay(m);
            nf_params.Lx        =   nf_p.Lx(m);
            nf_params.Ly        =   nf_p.Ly(m);
            nf_params.Lz        =   nf_p.Lz(m);
            nf_params.Lgamma    =   nf_p.Lgamma(m);
            nf_params.Neff      =   nf_p.Neff(m); 
            nf_params.fs_clock  =   data.HRM.fs_clock_ku_surf(m);

%% ---------------- DEFINE data to fit  -----------------------------------
% -------------------------------------------------------------------------
            switch cnf_p.mode
                case {'SAR','SARin','RAW','RMC','LR-RMC','FF-RAW','FF-RMC'}
                    max_power_wav       =   max(data.HRM.power_wav_filtered(m,start_sample:stop_sample));
                    fit_data            =   data.HRM.power_wav_filtered(m,start_sample:stop_sample)/max_power_wav;
%                 case {'SARin'}
%                     max_power_wav       =   max(data.SIN.power_wav_filtered(m,start_sample:stop_sample));
%                     fit_data            =   data.SIN.power_wav_filtered(m,start_sample:stop_sample)/max_power_wav;
            end            
            nf_params.max_power_wav=max_power_wav;
% Added EM 21.01.2016            
%% -------------- Define the look indexation ------------------------------
%--------------------------------------------------------------------------
       switch cnf_p.looks_index_method
            case 'Cris_old'
                N_Tlooks    =   nf_p.Neff(m);
                looks       =   (-floor((N_Tlooks-1)/2):floor(N_Tlooks/2)).';                
           case 'Norm_index'
                N_Tlooks    =   nf_p.Neff(m);
                looks       =   -floor((N_Tlooks-1)/2):floor(N_Tlooks/2);
                looks       =   (looks)/(nf_p.Npulses); % changed by EM: Model Chris l was per burst; not for stack (l should be normalized accordingly) 
           case 'Doppler_freq'
           %Generate the look index according to the relationship between
           %beam angles and Doppler
                switch cnf_p.fd_method
                    case 'exact'
                        %as a future option not available
                        looks=2.0/chd_p.wv_length_ku.*data.GEO.V_stack(m,1:nf_p.Neff(m)).*cos(data.HRM.beam_ang_stack(m,1:nf_p.Neff(m))).*...
                            data.HRM.pri_stack(m,1:nf_p.Neff(m)).*nf_p.Npulses;
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        %EM: 31.08.2016: Sentinel-6 Doppler angle
                        %corresponds to the complementary of Beam angle for
                        %SEOMs case Doppler angle still corresponds to
                        %angle between velocity vector & tangential
                        %direction                        
                        switch cnf_p.mission
                            case {'S6'}
                                beam_angles=...
                                    linspace(pi/2-data.HRM.doppler_angle_start_surf(m),...
                                    pi/2-data.HRM.doppler_angle_stop_surf(m),nf_p.Neff(m)).';                                                            
                            otherwise
                                beam_angles=...
                                    linspace(pi/2+data.HRM.doppler_angle_start_surf(m)-data.HRM.look_ang_start_surf(m),...
                                    pi/2+data.HRM.doppler_angle_stop_surf(m)-data.HRM.look_ang_stop_surf(m),nf_p.Neff(m)).';                                
                        end     
                        looks=2.0/chd_p.wv_length_ku.*data.GEO.V(m).*cos(beam_angles).*data.HRM.pri_surf(m).*nf_p.Npulses;                                
                        fd=2.0/chd_p.wv_length_ku.*data.GEO.V(m).*cos(beam_angles);
                end 
           case 'Look_angle'
               if data.HRM.pri_surf(m)==0
                   data.HRM.pri_surf(m)=data.HRM.pri_surf(idx_pri_non_zero);
               end
               delta_look_angle=asin(chd_p.wv_length_ku./(data.HRM.pri_surf(m).*2.0*...
                                    nf_p.Npulses*data.GEO.V(m))); 
               nf_params.delta_look_angle=delta_look_angle;
               switch cnf_p.look_ang_method
                    case 'exact'
                        %as a future option not available
                        looks=data.HRM.look_ang_stack(m,1:nf_p.Neff(m)).'/delta_look_angle;
                        look_angles = looks*delta_look_angle;
%                         save('Look_angles.mat','look_angles');
                    case 'approximate'
                        %compute the max and min beam angles from 
                        % Look and Doppler Angles information
                        looks=(linspace(data.HRM.look_ang_start_surf(m),...
                            data.HRM.look_ang_stop_surf(m),nf_p.Neff(m))/delta_look_angle).';
                        look_angles=(looks.*delta_look_angle).';
%                         save('looks.mat','look_angles');
%                         clear look_angles;
                        
                end 
               
       end 
       %!! Prova
       %looks=looks(1:2:nf_p.Neff(m));
       
%% ---------- Create the mask (as a matrix)--------------------------------
%--------------------------------------------------------------------------
nf_params.Doppler_mask=ones(nf_params.Neff,data.N_samples);
switch lower(cnf_p.Doppler_mask_cons_option)
%generated once avoiding to create it within the ML waveform
%function for computational load issues    
    case 'external'
        % Based on the information provided in the 
       if isfield(data.HRM,'Doppler_mask')           
           %to be in-line with L1B-processing
           if cnf_p.use_zeros_cnf == 1
               value_mask=1*0;
           else
               value_mask=NaN;
           end
           for i_look=1:nf_params.Neff
               if data.HRM.Doppler_mask(i_look,m)>0
                   %geometry corrections mask
                   %if data.HRM.Doppler_mask(i_look,m)<data.N_samples
                       nf_params.Doppler_mask(i_look,data.HRM.Doppler_mask(i_look,m):end)=value_mask;
                   %end
               else
                   %applying the removal of noise beams or gaps in-between
                   %all forced to NaN as we want to discard the beams
                   nf_params.Doppler_mask(i_look,:)=NaN;
               end
           end
       end
       % added by EM 21.03.2017 filter only the part of the waveform of
       % interest
       nf_params.Doppler_mask=nf_params.Doppler_mask(:,start_sample:stop_sample);
    case 'internal'
        % Constructed using information from Beam angles Doppler to compute
        % the slant range correction
        nf_params.range_history=compute_range_history_approx(data,nf_p,cnf_p,chd_p,m); 
        if isfield(data.HRM,'Doppler_mask')
            nf_params.Doppler_mask_beams_all_zeros_original= data.HRM.Doppler_mask(1:nf_params.Neff,m)==1; %mask always located in the first nf_params.Neff valid looks            
        else
            nf_params.Doppler_mask_beams_all_zeros_original=logical(zeros(1,nf_params.Neff));
        end
        switch lower(cnf_p.Doppler_mask_cons_internal)
            case 'l1b_like'
                %as done for L1B no consideration regarding where the epoch is
                %located
                if cnf_p.use_zeros_cnf == 1
                    value_mask=1*0;
                else
                    value_mask=NaN;
                end
                %construct the geometry corrections impact
                for i_look=1:nf_params.Neff
                    if nf_params.range_history(i_look)<data.N_samples
                        nf_params.Doppler_mask(i_look,(data.N_samples-nf_params.range_history(i_look)):end)=value_mask;
                    else
                        nf_params.Doppler_mask(i_look,1:end)=value_mask;
                    end
                end
                %add any possible masking induced by noisy beams removal or others
                nf_params.Doppler_mask(nf_params.Doppler_mask_beams_all_zeros_original,:)=value_mask;
                
                % Potential RMC mask for FF operation in RMC S-6
                switch cnf_p.mission
                    case {'S6'}
                        switch cnf_p.mode
                            case {'FF-RMC'}
                                %Values shall be included in chd (defined in the processor)
                                i_sample_start_chd = 1;
                                rmc_margin = 6;
                                rmc_mask = ones(nf_params.Neff,data.N_samples);
                                
                                rmc_mask(:,(data.N_samples/cnf_p.ZP/ 2 + i_sample_start_chd - 1 - rmc_margin)...
                                        * cnf_p.ZP + 1 : data.N_samples/cnf_p.ZP * cnf_p.ZP) = 0; 
            
                                nf_params.Doppler_mask=nf_params.Doppler_mask.*rmc_mask;
                        end
                end
                
                nf_params.Doppler_mask=nf_params.Doppler_mask(:,start_sample:stop_sample);
        end % switch type of internal mask construction     
        
end %switch type of Doppler mask construction

%nf_params.Doppler_mask=nf_params.Doppler_mask(1:2:nf_p.Neff(m),:);

%% ------------------ Pre-processing --------------------------------------
% -------------------------------------------------------------------------
% added by EM 29.03.2016: pre-processing stage to estimate the initial
% epoch as mid-point leading edge (some percent of peak power) & for
% estimation of the floor noise using derivative of the data itself
% detecting the leading edge (not fixing the window size initially)
if cnf_p.pre_processing
    %used to estimate the initial value of the epoch 
    %& eventually if desired the noise floor
    [peak_pow,idx_max_peak]=max(fit_data);
    idx_leading_edge=find(find(fit_data<=cnf_p.percent_leading_edge/100.0*peak_pow)<idx_max_peak, 1, 'last' );
    if isempty(idx_leading_edge)
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge=idx_max_peak;
    end
    fit_params_ini(1)=range_index(idx_leading_edge);
    
  
    %clear idx_leading_edge peak_pow idx_max_peak;
else
    %in case the noise adaptive processing estimation is activated specific
    %indication of the leading edge is required: force a value by default
    %87% of the peak
    [peak_pow,idx_max_peak]=max(fit_data);
    idx_leading_edge=find(find(fit_data<=0.87*peak_pow)<idx_max_peak, 1, 'last' );
    if isempty(idx_leading_edge)
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge=idx_max_peak;
    end
    %added by EM 16.05.2017
    if cnf_p.wvfm_portion_selec 
        if strcmpi((cnf_p.wvfm_portion_selec_type),'ref_height')
            fit_params_ini(1)=range_index(epoch_ref_DEM);
        end
        if strcmpi((cnf_p.wvfm_portion_selec_type),'CP4O')
            fit_params_ini(1)=data.CP4O.seedpw(m);
        end

    end
    
    
end

%Modified by EM 23.03.2017
%% --------------------------- ThN ----------------------------------------
% -------------------------------------------------------------------------
if cnf_p.Thn_flag
    switch lower(cnf_p.Thn_estimation_method)
        case 'external'
            %using a forced value provided in the configuration file
            nf_params.ThN        =ones(1,nf_p.Neff(m)).*cnf_p.external_Thn_value; %normalized as the waveform is normalized
        case 'fixed_window'
            %using a fixed window size providing the inital and final range
            %bins in the configuration file
            switch cnf_p.Thn_ML_SL_method
                case 'ML'
                    nf_params.ThN        =   mean(fit_data(cnf_p.Thn_w_first:cnf_p.Thn_w_first+cnf_p.Thn_w_width-1))*ones(1,nf_p.Neff(m));
                case 'SL'
                    %only possible when stack is available;
                    dumm=squeeze(data.HRM.beams_rng_cmpr_mask(m,1:nf_p.Neff(m),:));
                    dumm(~isfinite(dumm))=NaN;
                    dumm(dumm==0)=NaN;
                    max_mean_dumm=max(mean(dumm,1,'omitnan'));
                    for i_beam=1:nf_p.Neff(m)
                        nf_params.ThN(i_beam)        =   mean(dumm(i_beam,cnf_p.Thn_w_first:cnf_p.Thn_w_first+cnf_p.Thn_w_width-1)/max_mean_dumm,'omitnan');
                    end
                    clear dumm;
            end
            
      %      normalized_noise_estimation(m)=nf_params.ThN(1);
      %      if m==data.N_records
      %          mkdir([path_Results,'data/noise/']);
      %          save([path_Results,'data/noise/',file_id,'_noise_estimate.mat'],'normalized_noise_estimation');
      %          return
      %      else
      %          continue;
      %      end
            
        case 'adaptive'
            %not using a fixed window size for all the waveforms, based on
            %the derivative
            switch cnf_p.Thn_ML_SL_method
                case 'ML'
                    idx_noise=[];
                    temp_noise_thr=cnf_p.threshold_noise;
                    iter_noise=1;
                    while isempty(idx_noise) && iter_noise<cnf_p.max_iter_noise
                        idx_noise=find(abs(diff(fit_data))<=temp_noise_thr);
                        idx_noise=idx_noise(idx_noise<idx_leading_edge);
                        temp_noise_thr=temp_noise_thr*cnf_p.factor_increase_noise_iter;
                        iter_noise=iter_noise+1;
                    end
                    if iter_noise<cnf_p.max_iter_noise
                        nf_params.ThN        =   mean(fit_data(idx_noise))*ones(1,length(looks));
                        clear idx_noise tempo_noise_thr;
                    else
                        %take an average mean value of the previous Thermal
                        %noise estimation and keep it as reference for this
                        %record
                        if m~=1
                            nf_params.ThN         = mean(nf_params.ThN)*ones(1,length(looks));
                        else
                            %!!!need to be reviewed
                            %peak at the beginning of the window and so no
                            %space left for estimation of the noise floor
                            nf_params.ThN         = min(fit_data)*ones(1,length(looks));
                        end
                    end
                case 'SL'
                    %To be Implemented
                    nf_params.ThN        =zeros(1,nf_p.Neff(m));                    
            end
    end
else
    nf_params.ThN        =zeros(1,nf_p.Neff(m));
end

%nf_params.ThN=nf_params.ThN(1:2:nf_p.Neff(m));



       


%% ----------- DEFINE FITTING MODEL & CALL FITTING ROUTINE  ---------------
% -------------------------------------------------------------------------                        
            % -------------------------------------------------------------
            % --------- DEFINE FITTING FUNCTION MODEL ---------------------
            % -------------------------------------------------------------            
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p, chd_p, looks,func_f0);
                            case 'fmin'
                                fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p,chd_p,looks,func_f0)-fit_data).^2);
                        end
                        
                    case 'complete'
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p,chd_p,looks,func_f0,func_f1);
                            case 'fmin'
                                fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p,chd_p,looks,func_f0,func_f1)-fit_data).^2);
                        end
                end
            else
                switch cnf_p.fitting_fun_type
                    case 'lsq'
                        mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p, chd_p, looks);
                    case 'fmin'
                        fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p, chd_p, looks)-fit_data).^2);
                end
            end
            
            % -------------------------------------------------------------
            % --------- DEFINE TYPE OF MINIMIZATION PROCEDURE -------------
            % -------------------------------------------------------------
            switch cnf_p.mode
                case {'SAR','SARin','SARIN','HR','RMC','RAW','LR-RMC','FF-RMC','FF-RAW'}
                    switch cnf_p.fitting_fun_type
                        case 'lsq'
                            [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini,range_index,fit_data,cnf_p.fitting_options_lb,cnf_p.fitting_options_ub,cnf_p.fitting_options);
                        case 'fmin'
                            [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini,zeros(1,length(fit_params_ini)),cnf_p.fitting_options_lb,cnf_p.fitting_options);
                    end
                otherwise
                    error('No valid mode')
                    
            end
            
            %update the quality flag as bad if exit flag -1 or -2:            
            if flag==-1 
                %-1 algorithm terminated by the output function
                fit_res.Flag_quality(m) = 0;
            elseif flag==-2
                %-2 Problem is infeasible: the bounds lb and ub are
                %inconsistent (only for lsqcurvefit)
                fit_res.Flag_quality(m) = 0;
            else
                fit_res.Flag_quality(m) = 1;
            end
            
%% --------------------- DEFINE OUTPUT STRUCTURE --------------------------
% -------------------------------------------------------------------------        
            if cnf_p.lut_flag
                switch cnf_p.power_wfm_model
                    case 'simple'
                        [ml_wav,~,~]          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks,func_f0);
                        
                    case 'complete'
                        [ml_wav,~,~]          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks,func_f0,func_f1);
                end
            else
                [ml_wav,~,~]          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks);
            end
            
%             %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%             %save the modelled stack and the look angles
%             if m>=1 & m<=100
%                 modelled_stacks(1:nf_p.Neff(m),:,m) = stack_real;
%                 look_angle_stacks(1:nf_p.Neff(m),m) = look_angles;
%                 N_beams(m) = nf_p.Neff(m);
%                 latitude(m)=data.GEO.LAT(m);
%                 longitude(m)=data.GEO.LON(m);
%                 clear stack_real look_angle;
%                 %if (m==data.N_records) || mod(m,100)
%                 save([path_Results,'data/',file_id,'_Stacks_info.mat'],'modelled_stacks','look_angle_stacks','latitude','longitude','N_beams');
%                 %end
%             end
%             %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            correlation_fit         =   corrcoef(fit_data,ml_wav);
            fit_res.COR(m)          =   correlation_fit(1,2);
            fit_res.misfit(m)       =   100*sqrt(nanmean((abs(fit_data-ml_wav)).^2));
            
            % For Sea State CCI
            %Set to a bad quality retrieval if misfit is above a given
            %threshold
            if fit_res.misfit(m)> cnf_p.quality_flag_misfit_th
                fit_res.Flag_quality(m) = 0;          
            end
            
            fit_res.Epoch(m)    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP - 1; % include again the IF mask contribution
            %fit_res.Epoch(m)    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP; % include again the IF mask contribution
            
            fit_params(2)=abs(fit_params(2));
            
            % !!!!!!!!!!!!!!!!!! THE -1 VALUE MIGHT CAUSE SOME BIAS
            % CHECK IT !!!
            if cnf_p.rou_flag
                fit_res.rou(m)  =   fit_params(2); 
            else
                fit_res.Hs(m)   =   (fit_params(2)*4); % sigma_z = Hs/4
                if cnf_p.two_step_fitting
                    fit_res.rou(m)  =   NaN;
                end
            end
            if strcmp(cnf_p.multilook_option,'Cris_old')                
                fit_res.Pu(m)       =   10*log10(fit_params(3)*max_power_wav);
            else
                fit_res.Pu(m)       =   10*log10(max_power_wav);
            end            
            fit_res.flag(m)     =   flag;            
            if cnf_p.two_step_fitting
                ml_wav_1st_iter = ml_wav;
                COR_1st_iter=fit_res.COR(m);
                misfit_1st_iter = fit_res.misfit(m);
                Epoch_1st_iter=fit_res.Epoch(m);
                SWH_1st_iter=fit_res.Hs(m);
                Pu_1st_iter=fit_res.Pu(m);
            end            
            fit_res.flag_validity_L1B(m) =1; %fiting valid 
            
%% ----------- SET NEXT SEEDS / INITIAL VALUES OF FITTING PARAMS ----------
% -------------------------------------------------------------------------
            if cnf_p.initial_param_fit_feedback_flag
                fit_params_ini(1) = fit_params(1);
                fit_params_ini(3) = fit_params(3);
                
                                
                % for SWH estimation
                if cnf_p.ini_Hs_rou_sliding_win_opt==1 && m~=1
                    
                    % Use sliding window to estimate
                    %Compute the statistics for the whole window up to
                    %first estimated surface (for filtering purposes)
                    if cnf_p.rou_flag
                        average_total=nanmean(fit_res.rou(1:m-1));
                        std_total=nanstd(fit_res.rou(1:m-1));
                    else
                        average_total=nanmean(fit_res.Hs(1:m-1)/4);
                        std_total=nanstd(fit_res.Hs(1:m-1)/4);
                    end
                    if isnan(average_total)
                        average_total=fit_params(2);
                    end
                    if isnan(std_total)
                        std_total=0.0;
                    end
                    
                    %use sliding window to estimate the next initial seed
                    init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
                    if cnf_p.rou_flag
                        average_sliding=nanmean(fit_res.rou(init_sample:m-1));
                        std_sliding=nanstd(fit_res.rou(init_sample:m-1));
                    else
                        average_sliding=nanmean(fit_res.Hs(init_sample:m-1)/4);
                        std_sliding=nanstd(fit_res.Hs(init_sample:m-1)/4);
                    end
                    
                    if isnan(average_sliding)
                        average_sliding=fit_params(2);
                    end
                    if isnan(std_sliding)
                        std_sliding=0.0;
                    end
                    
                    %to discard effects of having a set of different
                    %consecutive non-valid SWH/rou estimates
                    if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
                        %take the average from the very beginning
                        fit_params_ini(2)=average_total;
                    else
                        current_value=fit_params(2);
                        if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
                            %if within the limits +- std of sliding include
                            %them in average
                            if cnf_p.rou_flag
                                fit_params_ini(2)=nanmean(fit_res.rou(init_sample:m));
                            else
                                fit_params_ini(2)=nanmean(fit_res.Hs(init_sample:m))/4.0;
                            end
                        else
                            %if not wihtin limits +- std sliding take the
                            %sliding average
                            fit_params_ini(2)=average_sliding;
                        end
                    end
                    
                else
                    % Use the actual estimation to feed subsquent surface
                    fit_params_ini(2)=abs(fit_params(2));
                    
                end
            end % end flag feedback
            % If not activated take the initial values set from CNF
                                             
            

%% ---------------------- TWO STEP FITTING --------------------------------
% -------------------------------------------------------------------------
            if cnf_p.two_step_fitting && fit_res.COR(m)<cnf_p.two_step_fitting_COR_threshold_rou/100
                %run the roughness fitting
                cnf_p.rou_flag=1;
                m_rou=m_rou+1;
                
                if m_rou==1
                    if strcmp(cnf_p.multilook_option,'Cris_old')
                        fit_params_ini_rou      =   [fit_params_ini(1) ini_p.rou fit_params_ini(3)];
                    else
                        fit_params_ini_rou      =   [fit_params_ini(1) ini_p.rou];
                    end
                else
                    if cnf_p.pre_processing
                        fit_params_ini_rou(1)  = range_index(idx_leading_edge);
                    else
                        fit_params_ini_rou(1)  = fit_params(1);
                    end
                end
                
                                
                % -------------------------------------------------------------
                % --------- DEFINE FITTING FUNCTION MODEL ---------------------
                % -------------------------------------------------------------
                if cnf_p.lut_flag
                    switch cnf_p.power_wfm_model
                        case 'simple'
                            switch cnf_p.fitting_fun_type
                                case 'lsq'
                                    mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p, chd_p,looks,func_f0);
                                case 'fmin'
                                    fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p, chd_p,looks,func_f0)-fit_data).^2);
                            end
                            
                        case 'complete'
                            switch cnf_p.fitting_fun_type
                                case 'lsq'
                                    mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p, chd_p,looks,func_f0,func_f1);
                                case 'fmin'
                                    fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p, chd_p,looks,func_f0,func_f1)-fit_data).^2);
                            end
                    end
                else
                    switch cnf_p.fitting_fun_type
                        case 'lsq'
                            mpfun               =   @(fit_params,x)ml_wave_gen(x, fit_params, nf_params, cnf_p, chd_p,looks);
                        case 'fmin'
                            fminfun               =   @(fit_params)sum((ml_wave_gen(range_index, fit_params, nf_params, cnf_p, chd_p,looks)-fit_data).^2);
                    end
                end
                
                % -------------------------------------------------------------
                % --------- DEFINE TYPE OF MINIMIZATION PROCEDURE -------------
                % -------------------------------------------------------------
                switch cnf_p.mode
                    case {'SAR','SARin','SARIN','HR','RMC','RAW','LR-RMC','FF-RAW','FF-RMC'}
                        switch cnf_p.fitting_fun_type
                            case 'lsq'
                                [fit_params,~,~,flag]     =   lsqcurvefit (mpfun,fit_params_ini_rou,range_index,fit_data,cnf_p.fitting_options_lb,cnf_p.fitting_options_ub,cnf_p.fitting_options);
                            case 'fmin'
                                [fit_params,~,flag]     =   fminsearchbnd (fminfun,fit_params_ini_rou,zeros(1,length(fit_params_ini_rou)),cnf_p.fitting_options_lb,cnf_p.fitting_options);
                        end
                    otherwise
                        error('No valid mode')
                        
                end
                
                
                
                % -------------------------------------------------------------------------
                % DEFINE OUTPUT STRUCTURE
                % -------------------------------------------------------------------------
                
                if cnf_p.lut_flag
                    switch cnf_p.power_wfm_model
                        case 'simple'
                            ml_wav          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks,func_f0);
                            
                        case 'complete'
                            ml_wav          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks,func_f0,func_f1);
                    end
                else
                    ml_wav          =   ml_wave_gen(range_index,fit_params,nf_params, cnf_p, chd_p,looks);
                end
                ml_wav_2nd_iter=ml_wav;
                correlation_fit         =   corrcoef(fit_data,ml_wav);
            
                misfit_2nd_iter       = 100*sqrt(nanmean((abs(fit_data-ml_wav)).^2));
                COR_2nd_iter          =   correlation_fit(1,2);
                Epoch_2nd_iter    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP - 1; % include again the IF mask contribution
                %Epoch_2nd_iter    =   fit_params(1) + cnf_p.IFmask_N*cnf_p.ZP; % include again the IF mask contribution
                rou_2nd_iter  =   fit_params(2);               
                if strcmp(cnf_p.multilook_option,'Cris_old')
                    Pu_2nd_iter       =   10*log10(fit_params(3)*max_power_wav);
                else
                    Pu_2nd_iter       =   10*log10(max_power_wav);
                end
                
                if correlation_fit(1,2)>COR_1st_iter                    
                    fit_res.COR(m)          =   COR_2nd_iter;
                    fit_res.misfit(m)  =    misfit_2nd_iter;
                    fit_res.Epoch(m)    =   Epoch_2nd_iter; % include again the IF mask contribution
                    fit_res.rou(m)  =   rou_2nd_iter;                    
                    fit_res.Hs(m) = NaN;                    
                    fit_res.Pu(m)       =   Pu_2nd_iter;                    
                    fit_res.flag(m)     =   flag;
                    
                    % For Sea State CCI
                    %Set to a bad quality retrieval if misfit is above a
                    %given threshold
                    if fit_res.misfit(m)> cnf_p.quality_flag_misfit_th
                        fit_res.Flag_quality(m) = 0;
                    end
                    %update the quality flag as bad if exit flag -1 or -2:
                    if flag==-1
                        %-1 algorithm terminated by the output function
                        fit_res.Flag_quality(m) = 0;
                    elseif flag==-2
                        %-2 Problem is infeasible: the bounds lb and ub are
                        %inconsistent (only for lsqcurvefit)
                        fit_res.Flag_quality(m) = 0;
                    else
                        fit_res.Flag_quality(m) = 1;
                    end  
                    
                    
                    
                else
                    ml_wav=ml_wav_1st_iter;
                end
                
                % -------------------------------------------------------------------------
                % SET NEXT SEED IF NEEDED
                % -------------------------------------------------------------------------

                if cnf_p.initial_param_fit_feedback_flag
                    
                    fit_params_ini_rou(1)  = fit_params(1);
                    fit_params_ini_rou(3)  = fit_params(3);
                    
                    
                    if cnf_p.ini_Hs_rou_sliding_win_opt==1 && m~=1
                        % Use a sliding window
                        
                        fit_params(2)=abs(fit_params(2));
                        
                        %Compute the statistics for the whole window from
                        %very beginning for outliers filtering
                        average_total=nanmean(fit_res.rou(1:m-1));
                        std_total=nanstd(fit_res.rou(1:m-1));
                        if isnan(average_total)
                            average_total=fit_params(2);
                        end
                        if isnan(std_total)
                            std_total=0.0;
                        end
                        
                        %use sliding window to estimate the next initial seed
                        init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
                        
                        average_sliding=nanmean(fit_res.rou(init_sample:m-1));
                        std_sliding=nanstd(fit_res.rou(init_sample:m-1));
                        
                        if isnan(average_sliding)
                            average_sliding=fit_params(2);
                        end
                        if isnan(std_sliding)
                            std_sliding=0.0;
                        end
                        
                        %to discard effects of having a set of different
                        %consecutive non-valid SWH/rou estimates
                        if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
                            %take the average from the very beginning
                            fit_params_ini_rou(2)=average_total;
                        else
                            current_value=fit_params(2);
                            if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
                                %if within the limits +- std of sliding include
                                %them in average
                                fit_params_ini_rou(2)=nanmean(fit_res.rou(init_sample:m));
                            else
                                %if not wihtin limits +- std sliding take the
                                %sliding average
                                fit_params_ini_rou(2)=average_sliding;
                            end
                        end
                        
                    else
                        fit_params_ini_rou(2)  = abs(fit_params(2));
                    end % End sliding window operation over SWH
                    
                end % End feedback feeding for initial values subsequent surface

            end % End 2-step fitting 


%% ---------------------- PLOT FITTING RESULTS ----------------------------
% -------------------------------------------------------------------------       
            if cnf_p.plot_fits_flag && ((mod(m,cnf_p.plot_fits_downsampling)==0) || (m==1)) && (data.GEO.LAT(m)>=cnf_p.plot_fits_lat_range(1) && data.GEO.LAT(m)<=cnf_p.plot_fits_lat_range(2))

                text_interpreter=get(0, 'defaultAxesTickLabelInterpreter'); %cnf_p.text_interpreter;
%                 colors=[86/255, 180/255, 233/255; 230/255, 159/255, 0/255; 0/255, 158/255, 115/255; 213/255, 94/255, 0/255; 0/255, 114/255, 178/255; 204/255, 121/255, 167/255; 240/255, 228/255, 66/255; 0/255, 0/255, 0/255]; % colorblind palette for plots
colors=[20/255, 115/255, 175/255; 209/255, 96/255, 47/255;  229/255, 158/255, 55/255; 25/255, 158/255, 116/255; 204/255, 121/255, 167/255; 240/255, 228/255, 66/255; 0/255, 0/255, 0/255];
                f1=figure;
                plt=plot(1:data.N_samples,data.HRM.power_wav_filtered(m,:)/max(data.HRM.power_wav_filtered(m,:)),'Color', colors(1,:), 'LineStyle', '-');
                plt.Color(4) = 0.5; % transparency
                hold on                
                if (cnf_p.two_step_fitting) && COR_1st_iter<cnf_p.two_step_fitting_COR_threshold_rou/100
                    plot(start_sample:stop_sample,ml_wav_1st_iter,'--r'); %ocean-like
                    plot(start_sample:stop_sample,ml_wav_2nd_iter,'.-g'); %ice-like 
                    h_leg=legend('L1b-Waveform', 'Analytical fit (SWH)','Analytical fit (Rou)','Location','northeastoutside');
                    pos_leg=get(h_leg,'Position'); 
       
                    textbox_string = [sprintf('Epoch (1st)= %.4g [r.b.]\nEpoch (2nd)=%.4g [r.b.]\nSWH=%.4g [m]\n$\rho$=%.4g [-]\nPu (1st)=.4g [-]\nPu (2nd)=%.4g [-]\nr (1st)=%.4g [\%]\nr (2nd)=%.4g [\%]\n', ...
                        Epoch_1st_iter, Epoch_2nd_iter, abs(SWH_1st_iter), abs(rou_2nd_iter), Pu_1st_iter, Pu_2nd_iter, COR_1st_iter*100, COR_2nd_iter*100)];
                    annotation('textbox', [pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],...
                        'String',textbox_string,...
                        'FitBoxToText','on',  'Interpreter',text_interpreter);
                    type_fit='2step';  

                else
                    plot(start_sample:stop_sample,ml_wav, 'Color', colors(1,:), 'LineStyle', '-')
                    if cnf_p.rou_flag
                        h_leg=legend('L1b-Waveform', 'Analytical fit','Location','northeastoutside');
                        pos_leg=get(h_leg,'Position');
                        annotation('textbox', [pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],...
                            'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['Rou = ' ,...
                            num2str(abs(fit_res.rou(m)),'%10.5e'), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [ - ]'],...
                            ['r = ', num2str(fit_res.COR(m)*100,5), '[\%]'], ['LSFlag = ', num2str(fit_res.flag(m))],...
                            ['$\sigma^0$ = ', num2str(fit_res.Pu(m)+data.HRM.s0_sf(m)), ' [ dB ]']},...
                            'FitBoxToText','on',  'Interpreter',text_interpreter);
                        type_fit='MSS';
                    else   
                        h_leg=legend('L1b-Waveform', 'Analytical fit','Location','northeastoutside');
                        pos_leg=get(h_leg,'Position');

                        if strcmp(text_interpreter, 'latex')
                            textbox_string = {['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['SWH = ' ,...
                                num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [dB]'], ...
                                ['$\rho$ = ', num2str(fit_res.COR(m)*100,5), ' [\%]'], ['LSFlag = ', num2str(fit_res.flag(m))],...
                                ['$\sigma^2$ = ', num2str(fit_res.Pu(m)+data.HRM.s0_sf(m)), ' [dB]'] };
                        else
                            textbox_string = {['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],['SWH = ' ,...
                                num2str(abs(fit_res.Hs(m)),4), ' [m]'],['Pu = ', num2str(fit_res.Pu(m),4), ' [dB]'], ...
                                ['\rho = ', num2str(fit_res.COR(m)*100,5), ' [%]'], ['LSFlag = ', num2str(fit_res.flag(m))],...
                                ['\sigma^0 = ', num2str(fit_res.Pu(m)+data.HRM.s0_sf(m)), ' [dB]'] };
                        end
                        annotation('textbox', [pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],...
                            'String',textbox_string,...
                            'FitBoxToText','on',  'Interpreter',text_interpreter, 'LineWidth', 1);
                        type_fit='SWH';
                    end
                    
                end
                if cnf_p.wvfm_portion_selec && strcmpi((cnf_p.wvfm_portion_selec_type),'ref_height')
                    plot(epoch_ref_DEM.*ones(1,2),[0 1],'-.k')
                end
                grid on
                xlabel('range bin','Interpreter',text_interpreter);
                if strcmp(text_interpreter, 'latex')
                    title(sprintf('wav. \\# %d (LAT: %.4g deg)', m, data.GEO.LAT(m)), 'Interpreter',text_interpreter); 
                else
                    title(sprintf('wav. # %d (LAT: %.4g deg)', m, data.GEO.LAT(m)), 'Interpreter',text_interpreter); 
                end
                axis([1 data.N_samples 0 1.0]);                
                
               
                if cnf_p.optional_ext_file_flag                    
                    print('-dpng ','-r0',[path_Results,'plots',filesep,'fitted_waveforms',filesep,file_id,'_',cnf_p.file_ext_string,'_W',num2str(m,'%04.0f'),'.png']);
                else
                    print('-dpng ','-r0',[path_Results,'plots',filesep,'fitted_waveforms',filesep,file_id,'_W',num2str(m,'%04.0f'),'.png']);
                end
                
%                 savefig([path_Results,'plots',filesep,'fitted_waveforms',filesep,file_id,'_W',num2str(m,'%04.0f'), '.fig']) 
                
                close(f1);
                
            end % end ploting option
            %disp(strcat('Time wvfm ',num2str(m),' [sec]:',num2str(toc(t_fitting_init))));
            if cnf_p.two_step_fitting
                cnf_p.rou_flag=0;
            end
    end % end for over surfaces    
end % end function

% Next seeds code from 26.03.2018
% %% ----------- SET NEXT SEEDS / INITIAL VALUES OF FITTING PARAMS ----------
% % -------------------------------------------------------------------------
% %             if fit_params(3) <0
% %                 %Pu
% %                 %force setting this parameter to 0
% %                 fit_params(3)=0;
% %             end
%             % forcing the amplitude fitting to 1
%             if cnf_p.initial_param_fit_feedback_flag
%                 fit_params_ini(3) = 1;
%                 
%                 % initial epoch fitting
% %                 % it has much no sense as right now take one provided by the
% %                 % peak retracker
% %                 if cnf_p.seed && m<length(data.GEO.LAT)-1
% %                     if (mean(data.MEA.seed)*0.75 <= data.MEA.seed(m+1)) && (data.MEA.seed(m+1)<  mean(data.MEA.seed)*1.25) % rank +-25% over the mean value
% %                         fit_params_ini(1)  =   data.MEA.seed(m+1);
% %                     end
% %                 else
% %                     fit_params_ini(1)  = fit_params(1);
% %                 end
%                 fit_params_ini(1)=fit_params(1);
% 
%                 % initial SWH/roughness seed
%                 fit_params(2)=abs(fit_params(2));
%                 if m~=1
%                     if cnf_p.rou_flag
%                         average_total=nanmean(fit_res.rou(1:m-1));
%                         std_total=nanstd(fit_res.rou(1:m-1));                        
%                     else
%                         average_total=nanmean(fit_res.Hs(1:m-1)/4);
%                         std_total=nanstd(fit_res.Hs(1:m-1)/4);                                              
%                     end
%                     
%                     if isnan(average_total)
%                         average_total=fit_params(2);
%                     end
%                     if isnan(std_total)
%                         std_total=0.0;
%                     end
%                     
%                     if cnf_p.ini_Hs_rou_sliding_win_opt==1
%                         %use sliding window to estimate the next initial seed
%                         init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
%                         if cnf_p.rou_flag
%                             average_sliding=nanmean(fit_res.rou(init_sample:m-1));
%                             std_sliding=nanstd(fit_res.rou(init_sample:m-1));                                                        
%                         else
%                             average_sliding=nanmean(fit_res.Hs(init_sample:m-1)/4);
%                             std_sliding=nanstd(fit_res.Hs(init_sample:m-1)/4);
%                         end
%                         
%                         if isnan(average_sliding)
%                             average_sliding=fit_params(2);
%                         end
%                         if isnan(std_sliding)
%                             std_sliding=0.0;
%                         end
%                         
%                         %to discard effects of having a set of different
%                         %consecutive non-valid SWH/rou estimates
%                         if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
%                             %take the average from the very beginning
%                             fit_params_ini(2)=average_total;
%                         else
%                             current_value=fit_params(2);
%                             if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
%                                 %if within the limits +- std of sliding include
%                                 %them in average
%                                 if cnf_p.rou_flag
%                                     fit_params_ini(2)=nanmean(fit_res.rou(init_sample:m));
%                                 else
%                                     fit_params_ini(2)=nanmean(fit_res.Hs(init_sample:m))/4.0;
%                                 end
%                             else
%                                 %if not wihtin limits +- std sliding take the
%                                 %sliding average
%                                 fit_params_ini(2)=average_sliding;
%                             end
%                         end
%                     else
%                         accumulated=accumulated+fit_params(2);
%                         fit_params_ini(2) = (accumulated/m_fitted);
%                     end
%                 else
%                     fit_params_ini(2)=fit_params(2);
%                 end
%             else
%                 
%                 fit_params_ini(1)=fit_params(1);
%                 fit_params_ini(2)=abs(fit_params(2));
%                 fit_params_ini(3)=fit_params(3);
%                 
%             end
                                             

%                 % -------------------------------------------------------------------------
%                 % SET NEXT SEED IF NEEDED
%                 % -------------------------------------------------------------------------
%                 %             if fit_params(3) <0
%                 %                 %Pu
%                 %                 %force setting this parameter to 0
%                 %                 fit_params(3)=0;
%                 %             end
%                 if cnf_p.initial_param_fit_feedback_flag
%                     % forcing the amplitude fitting to 1
%                     fit_params_ini(3) = 1;
%                     
%                     % initial epoch fitting
%                     % it has much no sense as right now take one provided by the
%                     % peak retracker
%                     if cnf_p.seed && m<length(data.GEO.LAT)-1
%                         if (mean(data.MEA.seed)*0.75 <= data.MEA.seed(m+1)) && (data.MEA.seed(m+1)<  mean(data.MEA.seed)*1.25) % rank +-25% over the mean value
%                             fit_params_ini(1)  =   data.MEA.seed(m+1);
%                         end
%                     else
%                         fit_params_ini_rou(1)  = fit_params(1);
%                     end
%                     
%                     % initial roughness seed
%                     fit_params(2)=abs(fit_params(2));
%                     if m~=1
%                         average_total=nanmean(fit_res.rou(1:m-1));
%                         std_total=nanstd(fit_res.rou(1:m-1));
%                         
%                         if isnan(average_total)
%                             average_total=fit_params(2);
%                         end
%                         if isnan(std_total)
%                             std_total=0.0;
%                         end
%                         
%                         if cnf_p.ini_Hs_rou_sliding_win_opt==1
%                             %use sliding window to estimate the next initial seed
%                             init_sample=max(m-(cnf_p.ini_Hs_rou_sliding_win_size-1),1);
%                             
%                             average_sliding=nanmean(fit_res.rou(init_sample:m-1));
%                             std_sliding=nanstd(fit_res.rou(init_sample:m-1));
%                             
%                             if isnan(average_sliding)
%                                 average_sliding=fit_params(2);
%                             end
%                             if isnan(std_sliding)
%                                 std_sliding=0.0;
%                             end
%                             
%                             %to discard effects of having a set of different
%                             %consecutive non-valid SWH/rou estimates
%                             if (average_sliding < average_total-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total) || (average_sliding > average_total+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_total)
%                                 %take the average from the very beginning
%                                 fit_params_ini_rou(2)=average_total;
%                             else
%                                 current_value=fit_params(2);
%                                 if ((current_value >= average_sliding-cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding) && (current_value <= average_sliding+cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold*std_sliding))
%                                     %if within the limits +- std of sliding include
%                                     %them in average
%                                     fit_params_ini_rou(2)=nanmean(fit_res.rou(init_sample:m));
%                                 else
%                                     %if not wihtin limits +- std sliding take the
%                                     %sliding average
%                                     fit_params_ini_rou(2)=average_sliding;
%                                 end
%                             end
%                         else
%                             accumulated_rou=accumulated_rou+fit_params(2);
%                             fit_params_ini_rou(2) = (accumulated_rou/m_rou);
%                         end
%                     else
%                         fit_params_ini_rou(2)=fit_params(2);
%                     end
%                 else
%                     fit_params_ini_rou(1)=fit_params(1);
%                     fit_params_ini_rou(2)=fit_params(2);
%                     fit_params_ini_rou(3)=fit_params(3);                    
%                 end


% L1BS_folder = 'C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\L1B_L1BS_GPOD_SARVATORE\agulhas\2013\';
% L2_filename = strcat('C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\inputs\L2_GPOD_SARVATORE\Agulhas\2013\',L1B_filename,'.nc');
% index_stack = m;
% [stack,beam_angles,beam_angles_tangent,lon,lat,time,index_stack,Stack_Beam_Index,Stack_Doppler_Beam_Angle,NLook_20Hz]=...
%         return_stack_GPOD(data.GEO.TAI.total,L2_filename,L1BS_folder,'index_stack_in',index_stack);
%     
% size_stack=size(stack(Stack_Beam_Index(1:NLook_20Hz),:));
% mask_from_L1BS = ones(size_stack(1),size_stack(2));
% mask_from_L1BS(isnan(stack(Stack_Beam_Index(1:NLook_20Hz),:))) = 0;
% figure; 
% subplot(2,2,1); imagesc(mask_from_L1BS); colormap('jet'); colorbar; title('Mask L1BS- GPOD');
% subplot(2,2,2); imagesc(nf_params.Doppler_mask); colormap('jet'); colorbar; title('Mask constructed at L2- GPOD');
% subplot(2,2,3.5); imagesc(-nf_params.Doppler_mask+mask_from_L1BS); colormap('jet'); colorbar; title('Mask L1BS - Mask constructed at L2 [GPOD]');
