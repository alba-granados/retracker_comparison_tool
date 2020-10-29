function [out_data,indx]=remove_outliers(y, varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This functions removes outliers of a given input vector, e.g., estimates
% of SSH, SWH, etc. The default outlier definition is based on Tukey's
% fences. Nans are placed instead.
%
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 17/09/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -y    =   input vector
% 
%      OPTIONAL
%       - outlier_percentil_low, _high: lower and upper quartiles
%       - IQR_times = nonnegative constant. If 1.5 "outlier", if 3 "far out"
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% Code based on /retracker/algorithms/outliers_by_percentiles.m by EM / isardSAT 
% 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:


p = inputParser;
p.addParamValue('type_outliers_removal','',@(x)ischar(x));
p.addParamValue('outlier_percentil_low', 25);
p.addParamValue('outlier_percentil_high', 75);
p.addParamValue('IQR_times', 1.5);

p.parse(varargin{:});
type_outliers_removal=p.Results.type_outliers_removal;
outlier_percentil_low=p.Results.outlier_percentil_low;
outlier_percentil_high=p.Results.outlier_percentil_high;
IQR_times=p.Results.IQR_times;
clear p;

if ~strcmp(type_outliers_removal, 'tukey_fence')
    fprintf('No such outlier definition in type_outliers_removal flag. Tukey s fences definition considered instead.');
    type_outliers_removal = 'tukey_fence';
end

if strcmp(type_outliers_removal, 'tukey_fence')
    IQR=iqr(y);
    percentiles = prctile(y,[outlier_percentil_low outlier_percentil_high]);
    indx=(y<(percentiles(1)-IQR_times*IQR) | y>(percentiles(2)+IQR_times*IQR));
    out_data=y;     %force the outliers as missing data or NaN
    out_data(indx)=NaN;
end
    
    
end