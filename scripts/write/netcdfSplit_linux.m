function [file_p2]=netcdfSplit_linux(file_origin_path, file_origin_name, interval_to_select, time_vector_name, varargin)	
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code creates a netcdf from input interval samples [initial final] of
% original netcdf track using ncks (sudo apt-get install nco) linux command
%
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 15/12/2020
% This software is built within the Sentinel-6 P4 L1 GPP project 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       - interval_to_select in initial index sample - end index sample
%       - time dimension name (time_vector_name)
% OUTPUT:
%       - 
% COMMENTS/RESTRICTIONS
% This script requires linux OS and nco software for linux installed.
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% 15/12/2020: Alba Granados, script creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 4
   file_ext = varargin{1};
else
    file_ext = '.nc';
end

file_p = [file_origin_path, file_origin_name, file_ext];

%% Prepare new netcdf file
ref_time_vector = ncread(file_p,time_vector_name);

name_time_dimension = split(time_vector_name, '/');
if any(strcmp(name_time_dimension(end), 'time_20_ku')) % reading L2 output product from retracker -> time TAI, not UTC
    ref_time_vector = ref_time_vector - 37; % convert to UTC to build filename with new sensing time
end
start_sens_new = datestr(datenum([2000, 1, 1, 0, 0, ref_time_vector(interval_to_select(1))]), 'yyyymmddTHHMMSS');
stop_sens_new = datestr(datenum([2000, 1, 1, 0, 0, ref_time_vector(interval_to_select(2))]), 'yyyymmddTHHMMSS');

new_file_name = [file_origin_name file_ext];
new_file_name(21:21+30) = [start_sens_new '_' stop_sens_new]; % new sensing time

file_p2 = [file_origin_path, new_file_name];

system(sprintf('ncks -d %s,%d,%d %s %s\n', name_time_dimension(end), interval_to_select(1), interval_to_select(2), file_p, file_p2));

end
