function [start_sample,stop_sample,flag_nadir_return_in_win] = pre_processing(wvfm,wd_surf,cnf_p,chd_p,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code is in charge of limiting the samples around which the
% retracking shall be considered (waveforms)
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 20/03/2017
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -wvfm    =   input waveform in the L2 processing   
%       -wd_surf =   measured window delay for the surface under analysis
%       -cnf_p = configuration parameters structure for L2 processing
%      OPTIONAL 
%       -wd_ref_DEM = window delay of the nadir return extracted from a DEM
% OUTPUT:
%       start_sample             =   first sample to be used
%       stop_sample              =   last sample to be used
%       flag_nadir_return_in_win =   flag indicating whether the nadir
%       return is within the range window or not
% RESTRICTIONS: 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: First option is to use a threshold peak retracker to locate the
% epoch: from this point the selection of the waveform samples is purely
% based on a noise threshold
% In a first approach lets consider the waveform extraction based on the
% primary peak retracker filtering
if(nargin<3 || nargin>(5+1*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('wd_ref_DEM',[]);
p.parse(varargin{:});
wd_ref_DEM=p.Results.wd_ref_DEM;
clear p;

flag_nadir_return_in_win=1;

N_samples=length(wvfm);

delta_tau=1.0/(chd_p.fs_clock_ku_chd*cnf_p.ZP); %sampling spacing in window delay

if cnf_p.pre_processing_ref_height
    %% --------------- Search for the portion of waveform from nadir ----------
    %Based on the DEM and height of the satellite compute the expected window
    %delay
    vec_wd=((1:N_samples)-cnf_p.ref_sample_wd).*delta_tau+wd_surf;
    
    if wd_ref_DEM<=max(vec_wd) && wd_ref_DEM>=min(vec_wd)
        epoch_ref_DEM=find(min(abs(vec_wd-wd_ref_DEM)),1,'first');
    else   
        %take the peak position
        [~,dumm]=max(wvfm);
        epoch_ref_DEM=dumm(1);
        flag_nadir_return_in_win=0;
    end    
    start_sample=max(epoch_ref_DEM-cnf_p.pre_processing_l_samples,1);
    stop_sample=min(epoch_ref_DEM+cnf_p.pre_processing_r_samples,N_samples);        
end



% %% --------------- Threshold retracker isolation --------------------------
% %Take the maximum peak in waveform: isolate the samples left and right to
% %the peak above a very low percentage of the peak
% [peak_pow,idx_max_peak]=max(wvfm);
% samples_above_threshold_left=find(wvfm>peak_pow*cnf_p.pre_processing_l_thres/100);
% samples_above_threshold_right=find(wvfm<=peak_pow*cnf_p.pre_processing_r_thres/100);
% start_sample=samples_above_threshold_left(find(samples_above_threshold_left<idx_max_peak(1),1,'last'));
% stop_sample=samples_above_threshold_right(find(samples_above_threshold_right>idx_max_peak(1),1,'first'));







% %% --------------- Primary peak retracker filtering -----------------------
% s=size(wvfm);
% if s(1)>1
%     N_samples=s(1);
% else
%     N_samples=s(2);
% end
% 
% %consecutive differences
% d_consecutive=diff(wvfm);
% 
% if s(1)>1
%     dumm=circshift(wvfm,-2,1)-wvfm;
% else
%     dumm=circshift(wvfm,-2,2)-wvfm;
% end
% %alternate differences
% d_alternate=dumm(1:end-2);
% clear dumm;
% 
% %---------- Compute the start & stop thresholds ---------------------------
% Th_start=sqrt(((N_samples-2)*sum(d_alternate.^2)-(sum(d_alternate)).^2)./((N_samples-2)*(N_samples-3)));
% Th_stop=sqrt(((N_samples-1)*sum(d_consecutive.^2)-(sum(d_consecutive)).^2)./((N_samples-1)*(N_samples-2)));
% 
% start_sample=find(d_consecutive>Th_start,1,'first');
% if isempty(start_sample)
%     start_sample=1;
% end
% stop_sample=find(d_consecutive<Th_stop,1,'last');
% if isempty(stop_sample)
%     stop_sample=1;
% end










end

