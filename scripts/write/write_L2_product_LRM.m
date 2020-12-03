function write_L2_product_LR(file,out,cnf_p,cst_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and geophysical correctiosn to be applied to the
% retracked range based on the info available in the L1B product and
% depending on the surface being observed
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
%       -file    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -out = structure of the output data for the L2 product
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       
% RESTRICTIONS: 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: No tracking changes until v2.0
% v2.0: Include explicitly the option to run analytical
% retracker in a 2step using cnf_p.analytical_type_of_fitting
%

%% ------------- NAME DEFINITION ------------------------------------------
%-------------- Define the output name for the output file ----------------
date_creation = datestr(now, '_yyyymmddTHHMMSS_');
aux=strsplit(date_creation,'_');
out.date_creation=aux(2);
clear aux;

%needs to be defined as using the final file naming convention
%name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','L2'),date_creation,'isd');
name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','L2'),'_isd');

if cnf_p.optional_ext_file_flag
    %name_L2_product=strrep(name_L2_product,'isd',strcat(cnf_p.file_ext_string,'_isd'));
    name_L2_product=strcat(name_L2_product,'_',cnf_p.file_ext_string);
end

%% ---------------------- GENERATION OF OUTPUT PRODUCT --------------------
switch cnf_p.output_product_format
    case 'mat' % Matlab output
        % .mat file
        save([file.resultPath 'data' filesep name_L2_product '.mat'],'out');        
    case 'nc'  % NetCDF output product
        prepare_NetCDF_L2_LRM([file.resultPath 'data' filesep name_L2_product '.nc'],out,cnf_p,cst_p); %can have a switch to define the variables if mission is S3 or S6
end


end

