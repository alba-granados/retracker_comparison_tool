function [fit_res]=analytical_retracker(data,cnf_p,chd_p,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs analytical retracker based on the original model developed 
% by Chirs Ray et al. in IEEE TGRS "SAR Altimeter Backscattered Waveform Model" 
% DOI:10.1109/TGRS.2014.23330423
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 17/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -data    =  data input structure for L2 processing
%       -cnf_p = configuration parameters structure for L2 processing
%      OPTIONAL
%       -LUT_f0_file =   file containing the look up table for the f0
%                       function parametrization (full path)
%       -LUT_f1_file =   file containing the look up table for the f1
%                       function parametrization (full path)
%       -path_Results=  full path to the folder results
%       
% OUTPUT:
%       -fit_res        =   structure of fitting results {Hs,Pu,SWH,sigma_rou}
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - gen_fit_params: creates a structure ini_p with the intial guess values of the
%                   fitting parameters
% - gen_nonfit_params_EM: generates and initializes the structure of
%                         non-fitting data parameters
% - fitting: implements the fitting of the data with the model analytical retracker 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - Need to optimize the code avoiding so many different data structures
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: Include optional input variables: LUTs files for f0 and f1
% functions
%% ---------------- Handling input variables ------------------------------
if(nargin<4 || nargin>(4+4*2))
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

%% --------------------- SEED GENERATION based on Pablo's Approach --------
% TBC and TBD

%% --------------------- INITIALIZATION OF FITTING PARAMETERS -------------
%--------------------------------------------------------------------------
% following structure defined by Cristina
[ini_p]=gen_fit_params(data,cnf_p);

%% -------------- INITIALIZATION OF THE NON-FITTING PARAMETERS ------------
%--------------------------------------------------------------------------
%create a structure with the necessary information to be used by the
%fitting algorithm for all the waveforms 
[nf_p] = gen_nonfit_params_EM (data,cnf_p,chd_p,cst_p); 

%% -------------- RUN THE FITTING -----------------------------------------
%--------------------------------------------------------------------------
[fit_res]       =   fitting_noise_method(data, cnf_p, chd_p , nf_p, ini_p,'LUT_f0_file',LUT_f0_file,'LUT_f1_file',LUT_f1_file,'path_Results',path_Results,'L1B_filename',L1B_filename);

end

