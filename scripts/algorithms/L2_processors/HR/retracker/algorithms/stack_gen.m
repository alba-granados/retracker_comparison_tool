function stack = stack_gen(x, l, fit_p, nf_p, cnf_p, chd_p, func_f0, func_f1)%varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to the modelled stack based on the Chris et al SAR waveform model
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V9 07/07/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       x           =   Ranges bin index matrix: varying in columns cte in
%                       rows
%       l           =   Looks index matrix: vrying in rows cte in columnts 
%       fit_p       =   Parameters to be fitted 
%       nf_p        =   parameters of need for the waveform generation, but not to be
%                       fitted
%       cnf_p       =   configuration parameters of L2 processor
% OUTPUT:
%       stack =   single look waveform for a given look and all range
%                       bins
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - sl_wave_gen: in charge of genrating the single look waveform based on
% the model proposed by Chris in IEEE TGRS "SAR Altimeter Backscattered Waveform Model"
% DOI:10.1109/TGRS.2014.23330423
% -
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Versions control:
% Based on the V9 version defined by Cristina Martin-Puig and updated
% according to revisit of the different expressions for antenna and
% radiation patterns and the f0 and f1 functions (first and second order power waveform model), 
% inclusion of the option to use external values based on LUTs for the f0 and f1: 
% It is extended in matrix formulation to generate the whole stack avoiding
% the external loop call
% 20.06.2017: Defintiion of the Bkl and Tkl for k<0

%delta_rb_original=cst_p.c_cst/(2.0*chd_p.bw_rx_ku_chd);%sampling spacing in the range dimension with no ZP

    %**********************************************************************
    %*************** Input optional parameters ****************************
    %**********************************************************************
%     p = inputParser;
%     defaultff0=0.0;
%     addOptional(p,'func_f0',defaultff0,@isnumeric);
%     addOptional(p,'func_f1',defaultff0,@isnumeric);
%     parse(p,varargin{:});
%     func_f0=p.Results.func_f0;
%     func_f1=p.Results.func_f1;
%     clear p

    if nargin < 8
        func_f1=0; % initialize altitude to 0;
    end
    if nargin <7
        func_f0=0; % initialize altitude to 0;
    end
    
    epoch=fit_p(1);
 
    if cnf_p.rou_flag
        sigmaz=cnf_p.Hs/4; rou=fit_p(2);
    else        
        rou=nf_p.rou; sigmaz=fit_p(2);
    end
    
      
    %**********************************************************************
    %% ******************** Samosa3 ***************************************
    %**********************************************************************   
    
    %------------dilation term---------------------------------------------
    switch cnf_p.window_type_a 
        case 'Adaptive'
            %PTR definition as function of SWH according to Dinardo
            nf_p.alphag_a = 1./(2*(0.4178+sqrt(0.0019+((sigmaz*4.0-0.9689)/30.6673).^2)).^2);
            
            %nf_p.alphag_a = 1./(2*(0.59676+sqrt(0.0019+((sigmaz*4.0-0.9689)/30.6673).^2)).^2);
            
        case 'Adaptive_S3'
            %Use the polynomial fitting of bias between ISR and GPOD to get
            %the SWH GPOD and use the PTR adequately
            SWH_ISR = sigmaz*4;
            bias_SWH = -0.36 + 0.08924.*(SWH_ISR) - 0.007074.*(SWH_ISR).^2;
            SWH_GPOD =  SWH_ISR + bias_SWH;
            
            nf_p.alphag_a = 1./(2*(0.4178+sqrt(0.0019+((SWH_GPOD-0.9689)/30.6673).^2)).^2);            
    end
    
    switch cnf_p.window_type_r
        case 'Adaptive'
            %PTR definition as function of SWH according to Dinardo
            nf_p.alphag_r = 1./(2*(0.4178+sqrt(0.0019+((sigmaz*4.0-0.9689)/30.6673).^2)).^2);
            
        case 'Adaptive_S3'
            %Use the polynomial fitting of bias between ISR and GPOD to get
            %the SWH GPOD and use the PTR adequately
            SWH_ISR = sigmaz*4;
            bias_SWH = -0.36 + 0.08924.*(SWH_ISR) - 0.007074.*(SWH_ISR).^2;
            SWH_GPOD =  SWH_ISR + bias_SWH;
            
            nf_p.alphag_r = 1./(2*(0.4178+sqrt(0.0019+((SWH_GPOD-0.9689)/30.6673).^2)).^2);             
    end
    
    
    g       =   sqrt(2*nf_p.alphag_a*nf_p.alphag_r./(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * l.^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));        
    
    %--------------Varibales  ---------------------------------------------
    k=(x - epoch); xl=nf_p.Lx.*l; 
    if rou==-1
        alpha_sigma=0.0; 
    else
        alpha_sigma=1/(nf_p.h^2*rou); 
    end    
    yk=nf_p.Ly.*abs(sqrt(k)); gk=g.*k;
    yk(k<0)=0;
    %            
    %**********************************************************************
    %************************ Antenna & Surface ***************************
    %**********************************************************************
    %Constant Term
    %Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2).*exp(-alpha_sigma * xl.^2).*exp(-nf_p.alphay * nf_p.yp^2).*exp(-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
    
    %consider the impact of compensating the antenna along-track
    if cnf_p.antenna_compensation_al
        Bkl=2.0*exp(-nf_p.alphay * nf_p.yp^2-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
    else
        Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-nf_p.alphay * nf_p.yp^2-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
    end    
    %Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-(alpha_sigma).*(yk).^2);
    
    switch cnf_p.power_wfm_model
        case 'complete'
            % Linear Term
            Tkl=(nf_p.Ly./abs(sqrt(k))).*(nf_p.alphay*nf_p.yp).*tanh(2*nf_p.alphay*nf_p.yp*yk)-...
                (nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
%             Tkl(k==0)=nf_p.Ly^2.*(nf_p.alphay*nf_p.yp)^2*2-(nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
%             Tkl(k<0)=0;
            Tkl(k<=0)=nf_p.Ly^2.*(nf_p.alphay*nf_p.yp)^2*2-(nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
            %small angle approximation for 0
            %Tkl=nf_p.Ly^2.*(nf_p.alphay*nf_p.yp)^2*2-(nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
    end
    
    %**********************************************************************
    %************************ Dilation functions **************************
    %**********************************************************************
    funcf0=0.*gk;
    if cnf_p.lut_flag  
        funcf0(gk==0)=1.077900274770464;
        indexes_1=gk>=cnf_p.LUT_ximin & gk<=cnf_p.LUT_ximax;
        indexes_2=floor((gk(indexes_1)-cnf_p.LUT_ximin)./cnf_p.LUT_step)+1;
        indexes_3=gk>cnf_p.LUT_ximax;
        funcf0(indexes_1)=func_f0(indexes_2); 
        funcf0(indexes_3)=sqrt(pi./(2.0*gk(indexes_3))).*(1+3./(8*(gk(indexes_3)).^2)+105./(16*8*(gk(indexes_3)).^4)); 
%         
        switch cnf_p.power_wfm_model
            case 'complete'
                funcf1=0.*gk; 
                funcf1(gk==0)=0.515224256147498; 
                funcf1(indexes_1)=func_f1(indexes_2);  
                funcf1(indexes_3)=-1.0*sqrt(pi.*gk(indexes_3)/8).*(1./((gk(indexes_3)).^2)+15./(8.*(gk(indexes_3)).^4));                
                
        end
    else   
        funcf0(gk~=0)=pi/4.0*sqrt(abs(gk(gk~=0))).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*besseli(1/4,1/4*(gk(gk~=0)).^2,1)); funcf0(gk==0)=1.077900274770464;%2^(1/4)*gamma(5/4);
        switch cnf_p.power_wfm_model
            case 'complete'
              funcf1=zeros(1,length(x)); 
              funcf1(gk~=0)=-1.0*pi/8.0*(abs(gk(gk~=0)).^(3/2)).*(besseli(1/4,1/4*(gk(gk~=0)).^2,1)-besseli(-3/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)-besseli(3/4,1/4*(gk(gk~=0)).^2,1))); 
              funcf1(gk==0)=0.515224256147498;%gamma(3/4)/(2.0*(2)^(1/4));
        end
        
    end
    
    %**********************************************************************
    %************************ Common Cte param ****************************
    %**********************************************************************
    %K_cte=sqrt(2*pi)*sqrt(1/(2*alphag_r))*sqrt(1/(2*alphag_a))*;
   
    
    
    %**********************************************************************
    %************************ Power Waveform ******************************
    %**********************************************************************
    switch cnf_p.power_wfm_model
        case 'simple'
            %stack         =   sqrt(nf_p.fs_clock/chd_p.bw_rx_ku_chd)*sqrt(g).*Bkl.*funcf0;
            stack         =   (nf_p.fs_clock/chd_p.bw_rx_ku_chd)*sqrt(g).*Bkl.*funcf0;
        case 'complete'
            %stack         =   *sqrt(g).*Bkl.*(funcf0+Tkl.*g.*((sigmaz/nf_p.Lz)^2).*funcf1);            
            stack         =   (nf_p.fs_clock/chd_p.bw_rx_ku_chd)*sqrt(g).*Bkl.*(funcf0+Tkl.*g.*((sigmaz/nf_p.Lz)^2).*funcf1);            
    end
    
    % **************** Zero-padding Impact ********************************
    switch cnf_p.range_index_method
        case 'conventional'
%             switch cnf_p.mission
%                 case {'S6'}
%                     stack=sqrt(cnf_p.ZP).*stack;    
%             end
%             %stack=stack;    
            stack=sqrt(cnf_p.ZP).*stack;    
    end
       
end
