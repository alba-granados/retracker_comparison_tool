function [sl_lrm_wav] = sl_lrm_wave_gen_LA (t, fit_p, nf_p)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to generate Laiba's LRM waveform model. The code is optimized for the
% fitting of epoch, sigmaz(Hs) and Pu
% -------------------------------------------------------------------------
% 
% Author:           Cristina Martin-Puig / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT v1 07/01/2014
%
% This software is built under the Jason-CS contract
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       t           =   Range bin index in time 
%       fit_p       =   Parameters to be fitted 
%       nf_p        =   parameters of need for the waveform generation, but not to be
%                       fitted
%       cnf_p       =   configuration parameters of L2 processor
% OUTPUT:
%       sl_LRM_wav  =   LRM waveform 
% -------------------------------------------------------------------------  

    if length(fit_p)==3
        % we are working in MLE3 mode
        epoch   =   fit_p(1);
        Hs      =   fit_p(2);
        Pu      =   fit_p(3);
        xi      =   nf_p.xi;
        
    elseif length(fit_p)==4
        % we are working in MLE4 mode
        epoch   =   fit_p(1);
        Hs      =   fit_p(2);
        Pu      =   fit_p(3);
        xi      =   fit_p(4);
    else
        error('error in the number of fitted params');
    end

    sigma_c     =   sqrt((0.513*nf_p.rt)^2+(0.5*Hs/nf_p.c)^2); %Hs = 4*sigma_s
%     sigma_c     =   sqrt((1/.513*nf_p.rt)^2+(0.5*Hs/nf_p.c)^2); %Hs = 4*sigma_s
    %     sigma_c     =   sqrt((0.425*nf_p.rt)^2+(0.5*Hs/nf_p.c)^2);
    
    %for 0.425
    gamma      =   (sin(nf_p.beamwidth).^2)   /(2*log(2));          % <--- from Amarouche paper
%     gamma      =   (sin(nf_p.beamwidth/2).^2) *2/log(2);            % <--- from S3 L2 ADAS <--- SAME AS THE PREVIOUS ONE
%     gamma      =   (sin(nf_p.beamwidth).^2)   /(sqrt(2)*log(2));    % <--- THIS IS THE BEST SOLUTION, Amarouche * sqrt(2)
%     gamma      =   (sin(nf_p.beamwidth).^2)   /(sqrt(2)*2*log(2));    % <--- THIS IS THE BEST SOLUTION, Amarouche / sqrt(2)
%     gamma      =   (sin(nf_p.beamwidth).^2) / (1/0.513);               % using the approximation in the paper instead of these logs
    
    %for 0.513
%     gamma      =   (sin(nf_p.beamwidth).^2)*1.052676;             % <--- from Amarouche paper
%     gamma      =   (sin(nf_p.beamwidth).^2)*1.052676/sqrt(2);     %     <--- THIS IS THE BEST SOLUTION, Amarouche * sqrt(2)
%     gamma      =   (sin(nf_p.beamwidth).^2)*1.052676/(2);         %     <--- THIS IS THE BEST SOLUTION, Amarouche * 2

    
    a           =   exp(-4*sin(xi).^2./gamma);
    
%     b           =   cos(2*nf_p.xi)-(sin(2*nf_p.xi).^2./gamma);
%     
%     alpha1       =   b.*(4*nf_p.c/(gamma*nf_p.h));

    
    
    %----------------------------------------------------------------------
    %Resizing the x-axis. The reference must be the reality, not the sampling.
% 	t           =   t*(320e6/nf_p.BWClock);
%     t           =   t/(320e6*nf_p.rt);
%     t           =   t*nf_p.rt;
    %----------------------------------------------------------------------
    
%     delta       =   (4*nf_p.c*cos(2*xi))/(gamma*nf_p.h); % <--- THIS HAS BEEN CHANGED
    delta       =   (4*nf_p.c*cos(2*xi))/(gamma*(nf_p.h*(1+nf_p.h/nf_p.semi_major_axis))); % <--- THIS IS THE NEW VALUE, FOR NOW [18.02.2015]
    
%     beta        =   (4/gamma)*sqrt(nf_p.c/nf_p.h)*sin(2*xi); % <--- THIS HAS BEEN CHANGED
    beta        =   (4/gamma)*sqrt(nf_p.c/(nf_p.h*(1+nf_p.h/nf_p.semi_major_axis)))*sin(2*xi); % <--- THIS IS THE NEW VALUE, FOR NOW [18.02.2015]
    
%     alpha       =   delta - beta.^2/4;
    alpha1       =   delta - (beta.^2)/8;
    
%     v1          =   alpha.*((t-epoch)*nf_p.rt-0.5*alpha*sigma_c);
    v1          =   alpha1.*((t-epoch)*nf_p.rt-0.5*alpha1*sigma_c.^2);
    
%     u1          =   ((t-epoch)*nf_p.rt-alpha*sigma_c^2)/(sqrt(2)*sqrt(sigma_c));
    u1          =   ((t-epoch)*nf_p.rt-alpha1*sigma_c.^2)/(sqrt(2)*sigma_c);
    
    alpha2       =   delta;
    
%     v2          =   alpha.*((t-epoch)*nf_p.rt-0.5*alpha*sigma_c);
    v2          =   alpha2.*((t-epoch)*nf_p.rt-0.5*alpha2*sigma_c.^2);
    
%     u2          =   ((t-epoch)*nf_p.rt-alpha*sigma_c^2)/(sqrt(2)*sqrt(sigma_c));
    u2          =   ((t-epoch)*nf_p.rt-alpha2*sigma_c.^2)/(sqrt(2)*sigma_c);
    
    sl_lrm_wav  =   a*Pu*(exp(-v1).*(1+erf(u1))-0.5*exp(-v2).*(1+erf(u2)))+nf_p.TN;
    
    

end
