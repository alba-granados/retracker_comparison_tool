function [range_history] = compute_range_history_approx(data,nf_p,cnf_p,chd_p,m)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to compute the approximate slant range history variation as a function of Doppler for the
% given look/beam angle information of the corresponding stack (using theoretical approximation)
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
%       m           =   Surface of interest
%   OPTIONAL:
%       
%
% OUTPUT:
%       range_history = range history for all the beams within the stack of
%       interest (value in number of samples)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% V1: From the start and stop look and Doppler angle satellite (original definition in the DPM)
% the beam angles is computed and the Doppler frequency as well
%


%% ------------ Compute the Doppler frequency -----------------------------

%-------------- Compute the Beam angles -----------------------------------
switch cnf_p.mission
    case {'S6'}
        %in Sentinel-6 the Doppler angle is defined as the complementary of
        %the Beam angle
        beam_angles=linspace(pi/2-data.HRM.doppler_angle_start_surf(m),...
            pi/2-data.HRM.doppler_angle_stop_surf(m),nf_p.Neff(m)).';    
    case {'CS2','S3'}
        switch cnf_p.L1proc
            case {'GPOD'}
                beam_angles=pi/2-data.HRM.look_ang_stack(m,:);
            case {'ESA'}
                beam_angles=linspace(pi/2-data.HRM.look_ang_start_surf(m),...
                    pi/2-data.HRM.look_ang_stop_surf(m),nf_p.Neff(m)).';
            otherwise                
                beam_angles=linspace(pi/2+data.HRM.doppler_angle_start_surf(m)-data.HRM.look_ang_start_surf(m),...
                    pi/2+data.HRM.doppler_angle_stop_surf(m)-data.HRM.look_ang_stop_surf(m),nf_p.Neff(m)).';
        end
end
        
fd=2.0/chd_p.wv_length_ku.*data.GEO.V(m).*cos(beam_angles);      

%-------------- Azimuth chirp rate ----------------------------------------
v_g=data.GEO.V(m)/nf_p.alpha(m); %ground velocity projection
v_e=sqrt(data.GEO.V(m)*v_g); %effective velocity
R=data.MEA.win_delay(m)*chd_p.wv_length_ku/2.0;
Ka=2*v_e.^2/(chd_p.wv_length_ku*R);

%--------------- Slant range history variation ----------------------------
switch cnf_p.mission
    case {'S3A','S3B','S3'}
        switch cnf_p.L1proc
            case {'ESA','DeDop'}
                range_history=(chd_p.wv_length_ku.*fd.^2)./(4*Ka);
            otherwise
                range_history=data.GEO.H_rate(m).*fd/Ka+(chd_p.wv_length_ku.*fd.^2)./(4*Ka); %
        end
    case {'CS2','CR2'}
        switch cnf_p.L1proc
            case {'GPOD','ESA'}
                range_history=(chd_p.wv_length_ku.*fd.^2)./(4*Ka);
            otherwise
                range_history=data.GEO.H_rate(m).*fd/Ka+(chd_p.wv_length_ku.*fd.^2)./(4*Ka); %
        end
    case {'S6'}
        range_history=data.GEO.H_rate(m).*fd/Ka+(chd_p.wv_length_ku.*fd.^2)./(4*Ka); %
end
delta_range_spacing=chd_p.wv_length_ku/(2*cnf_p.ZP*chd_p.fs_clock_ku_chd); %spacing in meters
range_history=floor(range_history/delta_range_spacing);


end % end function