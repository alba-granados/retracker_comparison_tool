function ini_p = gen_fit_params (data,cnf_p) 
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This function initializates fitting parameters structure
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT v1 15/06/2016
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
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
%

    % initialitzate fitting parameters    
    if ~cnf_p.seed
        ini_p.Epoch     =   cnf_p.ini_Epoch;
    else
        ini_p.Epoch     =   data.MEA.seed(1);    
    end
    
    ini_p.Hs        =   cnf_p.ini_Hs;
    ini_p.sigmaz    =   ini_p.Hs/4;
    if (strcmp(cnf_p.multilook_option,'Cris_old'))
        ini_p.Pu        =   cnf_p.ini_Pu;  
    end
    ini_p.rou       =   cnf_p.init_rou;
%     %added by EM 16.03.2016 include optimization of the noise thermal
%     %estimation window
%     if cnf_p.fit_noise
%        ini_p.Thn_1st=cnf_p.ini_Thn_1st;
%        ini_p.Thn_wdth=cnf_p.ini_Thn_wdth;
%     end
    
%     % if LUT initializate LUT tables
%     if (cnf_p.lut_flag)
%         % Look up tables initialization 
%         % Generate LUT configuration file
%         switch cnf_p.pdf_surf
%             case 'Gaussian'
%                 order = 1; % this generates up to f_order(xi)
%                 gen_LUT(order,[cnf_p.LUT_ximin:((cnf_p.LUT_ximax-cnf_p.LUT_ximin)/(cnf_p.LUT_NLUT)):cnf_p.LUT_ximax],cnf_p);
%             case 'non-Gaussian'
%                 print('WARNING: Code not optimized for non-Gaussian')
%                 order = 4; % this generates up to f_order(xi)
%                 gen_LUT(order,[cnf_p.LUT_ximin:((cnf_p.LUT_ximax-cnf_p.LUT_ximin)/(cnf_p.LUT_NLUT)):cnf_p.LUT_ximax],cnf_p);
%             otherwise
%                 error('Not a valid surface pdf')
%         end
%     end

        
end