function [fit_res]=threshold_retracker(data,cnf_p,chd_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs a simple threshold retracker: based on determining the
% leading edge as the % of the peak
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 01/07/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -data    =  data input structure for L2 processing
%       -cnf_p = configuration parameters structure for L2 processing
%      OPTIONAL
%       
% OUTPUT:
%       -fit_res        =   structure of fitting results {Epoch,sigma0}
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
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
% V1.1: Include the filtering of the portion of the waveform of interest
%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(3+2*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('path_Results',{''},@(x)ischar(x));
p.addParamValue('L1B_filename',{''},@(x)ischar(x));
p.parse(varargin{:});
path_Results=char(p.Results.path_Results);
L1B_filename=char(p.Results.L1B_filename);
clear p;
%% --------------------- RUN DETECTOR -------------------------------------
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% IF MASK - Small aliasing neglection
% Due to the fact that the IFMask is not perfect and it does not go to zero
% in the stop band some aliasing is introduced in the first few samples.
% It' been agreed among the community that these samples ~12 first and last
% shall not be considered for the fitting of the waveform
% -------------------------------------------------------------------------
switch cnf_p.mode
    case {'SAR','SARin','RAW','LR-RMC','RMC','FF-RAW','FF-RMC'}
        data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
%     case {'RMC'}
%         data.HRM.power_wav_filtered    =   data.HRM.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:cnf_p.Nrmc,:)';
%     case {'SARin'}
%         data.SIN.power_wav_filtered    =   data.SIN.power_wav(cnf_p.IFmask_N*cnf_p.ZP+1:end-(cnf_p.IFmask_N*cnf_p.ZP),:)';
    otherwise
        error('not a valid mode');
end
%--------------------------------------------------------------------------
%-------------------DETECTOR-BASED ----------------------------------------
%--------------------------------------------------------------------------
for m = 1:data.N_records
    
    
    %--------------------------------------------------------------------------
    %---------------------- waveform portion selection ------------------------
    %--------------------------------------------------------------------------
    
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
    
    
    %----------------- Discarding samples beginning and end -----------------
    if cnf_p.wvfm_discard_samples
        start_sample=max(start_sample,1+cnf_p.wvfm_discard_samples_begin);
        stop_sample=min(stop_sample,data.N_samples-cnf_p.wvfm_discard_samples_end);
    end
    
    
    % added EM: 20.09.2016
    %checking waveform valid due to L1B processing related to noise
    %filtering outer beams
    if ~any(~isnan(data.HRM.power_wav_filtered(m,start_sample:stop_sample))) || ~any(data.HRM.power_wav_filtered(m,start_sample:stop_sample))
        fit_res.Epoch(m)    =   NaN;
        fit_res.Pu(m)       =   NaN;
        fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
        continue;
    end
    
    [peak_pow,idx_max_peak]=max(data.HRM.power_wav_filtered(m,start_sample:stop_sample),[],2);
    idx_max_peak=start_sample-1+idx_max_peak(1);
    dumm=start_sample-1+find(data.HRM.power_wav_filtered(m,start_sample:stop_sample)>cnf_p.th_retracker.percentage_peak/100.0*peak_pow, 1);
    if ~isempty(dumm)
        %position by simple linear interpolation
        if dumm>1
            idx_leading_edge_ISD=dumm-1+((cnf_p.th_retracker.percentage_peak/100.0*peak_pow-data.HRM.power_wav_filtered(m,dumm-1))/(data.HRM.power_wav_filtered(m,dumm)-data.HRM.power_wav_filtered(m,dumm-1)));
            if idx_leading_edge_ISD<0 || idx_leading_edge_ISD>data.N_samples
                idx_leading_edge_ISD=dumm; %it could happen that power at dumm and dumm-1 is almost the same (linear interpolation does not work properly)
                %leading to a really high value of epoch.
            end
        else
            idx_leading_edge_ISD=dumm;
        end
    else
        %if there is no leading edge or the waveform has displaced that
        %much from the window to the left select the peak as leading
        %edge
        idx_leading_edge_ISD=idx_max_peak;
    end
    fit_res.Epoch(m)    =   idx_leading_edge_ISD + cnf_p.IFmask_N*cnf_p.ZP - 1;
    fit_res.Pu(m)       =   10*log10(peak_pow); %dB
    fit_res.flag_validity_L1B(m) = 1;
    fit_res.flag_nadir_return_in_win(m)=flag_nadir_return_in_win;
end

end

