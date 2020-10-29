function [start_sample,stop_sample,flag_nadir_return_in_win,epoch_ref_DEM] = wvfm_portion_selec(wvfm,wd_surf,cnf_p,chd_p,m,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code is in charge of limiting the samples around which the
% retracking shall be considered (waveforms)
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
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
%       - m = waveform index
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
% v1.0: First option is to use the difference between window delay from DEM
% w.r.t measured window delay
% v1.1: Add two more options to select the portion of the waveform on the
% basis of the peak using 1) fixed value of samples on left and right and
% 2) using the samples on the left and right based on a threshold w.r.t
% peak
if(nargin<4 || nargin>(5+3*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('wd_ref_DEM',[]);
p.addParamValue('path_Results',{''},@(x)ischar(x));
p.addParamValue('L1B_filename',{''},@(x)ischar(x));
p.parse(varargin{:});
wd_ref_DEM=p.Results.wd_ref_DEM;
path_Results=char(p.Results.path_Results);
L1B_filename=char(p.Results.L1B_filename);
clear p;


flag_nadir_return_in_win=1;

N_samples=length(wvfm);

delta_tau=1.0/(chd_p.fs_clock_ku_chd*cnf_p.ZP); %sampling spacing in window delay

switch (cnf_p.wvfm_portion_selec_type)
    %% --------------- Search for the portion of waveform from nadir ----------
    %Based on the DEM (corrected for difference geoid ellipsoid) and height of the satellite compute the expected window
    %delay
    case 'CP4O'
      
        epoch_ref_DEM=wd_ref_DEM;
        start_sample=max(wd_ref_DEM-cnf_p.wvfm_portion_selec_l_samples,1);
        stop_sample=min(wd_ref_DEM+cnf_p.wvfm_portion_selec_r_samples,N_samples);
        

    case 'ref_height'
        
        vec_wd=((1:N_samples)-cnf_p.ref_sample_wd).*delta_tau+wd_surf;
        
        if wd_ref_DEM<=max(vec_wd) && wd_ref_DEM>=min(vec_wd)
            %epoch_ref_DEM=find(min(abs(vec_wd-wd_ref_DEM)),1,'first');
            [~,epoch_ref_DEM]=min(abs(vec_wd-wd_ref_DEM));
        else
            %take the peak position
            [~,dumm]=max(wvfm);
            epoch_ref_DEM=dumm(1);
            flag_nadir_return_in_win=0;
            %disp(strcat({'Nadir return for Waveform #: '},num2str(m),{' outside the window, taking portion around max. peak'}))
        end
        
        
        %combine the information with the closest peak & valley info to
        %select the samples of interest around it
        [pks,loc_pks,w_pks,p_pks]=findpeaks(wvfm./max(wvfm),'MinPeakProminence',cnf_p.peak_prominence_norm); %peaks        
        [vll,loc_vll]=findpeaks((max(wvfm)-wvfm)); %valleys
        [~,idx_peak]=min(abs(loc_pks-epoch_ref_DEM));

        % ----------- Compute the start sample ----------------------------
        %take the 1st peak to the left & valley to the left of it
        idx_peak_left=find(loc_pks<epoch_ref_DEM,1,'last');
        if ~isempty(idx_peak_left)
            %find the first valley to left & set first left sample
            idx_valley_left=find(loc_vll<loc_pks(idx_peak_left),1,'last');
            if ~isempty(idx_valley_left)
                start_sample=max(loc_vll(idx_valley_left)-cnf_p.wvfm_portion_selec_l_samples,1);
            else
                start_sample=max(loc_pks(idx_peak_left)-cnf_p.wvfm_portion_selec_l_samples,1);
            end
            
        else
            %take the reference epoch as "peak"
            idx_valley_left=find(loc_vll<epoch_ref_DEM,1,'last');
            if ~isempty(idx_valley_left)
                start_sample=max(loc_vll(idx_valley_left)-cnf_p.wvfm_portion_selec_l_samples,1);
            else
                start_sample=max(loc_pks(idx_peak_left)-cnf_p.wvfm_portion_selec_l_samples,1);
            end
        end
        
        %-------------- Compute the stop sample ---------------------------
        %tale the first right peak
        idx_peak_right=find(loc_pks>epoch_ref_DEM,1,'first');
        if ~isempty(idx_peak_right)
            %find the first valley to the right of left peakpeak and below the second peak
            if ~isempty(idx_peak_left)
                idx_valley_right=find(loc_vll>=loc_pks(idx_peak_left) & loc_vll<loc_pks(idx_peak_right));
            else
                %as the epoch_ref_DEM as "peak"
                idx_valley_right=find(loc_vll>=epoch_ref_DEM & loc_vll<loc_pks(idx_peak_right));
            end
            if ~isempty(idx_valley_right)
                [~,min_vll_right_idx]=min(wvfm(loc_vll(idx_valley_right)));
                idx_valley_right=idx_valley_right(min_vll_right_idx);
                stop_sample=min(loc_vll(idx_valley_right)+cnf_p.wvfm_portion_selec_r_samples,round(loc_pks(idx_peak_right)-w_pks(idx_peak_right)/2));
            else
                %stop_sample=min(loc_pks(idx_peak)+cnf_p.wvfm_portion_selec_r_samples,N_samples);
                stop_sample=min(epoch_ref_DEM+cnf_p.wvfm_portion_selec_r_samples,round(loc_pks(idx_peak_right)-w_pks(idx_peak_right)/2));
            end
        else
            loc_valley_right=find(loc_vll>=epoch_ref_DEM,1,'first');
            if ~isempty(loc_valley_right)
                stop_sample=min([loc_vll(loc_valley_right)+cnf_p.wvfm_portion_selec_r_samples,N_samples]);
            else
                stop_sample=min([epoch_ref_DEM+cnf_p.wvfm_portion_selec_r_samples,N_samples]);
            end
        end
        
        
%         % ALTERNATIVE:
%         %take the closest peak to ref_epoch and define the halfbeamwidth
%         %peak left and right + left and right margin samples
%         [~,idx_peak_closest_ref_epoch]=min(abs(loc_pks-epoch_ref_DEM));
%         start_sample=max(loc_pks(idx_peak_closest_ref_epoch)-ceil(w_pks(idx_peak_closest_ref_epoch)/2)-cnf_p.wvfm_portion_selec_l_samples,1);
%         stop_sample=min(loc_pks(idx_peak_closest_ref_epoch)+ceil(w_pks(idx_peak_closest_ref_epoch)/2)+cnf_p.wvfm_portion_selec_r_samples,N_samples);
        
        
    %% --------------- Maximum peak and window around ---------------------
    %Taking a window around left and right samples
    case 'peak_win'  
        [peak_pow,idx_max_peak]=max(wvfm);
        start_sample=max(idx_max_peak(1)-cnf_p.wvfm_portion_selec_l_samples,1);
        stop_sample=min(idx_max_peak(1)+cnf_p.wvfm_portion_selec_r_samples,N_samples);
    %% ---------- Maximum peak and window around base on threshold --------
    %Taking a window around left and right samples
%     case 'peak_thresh'  
%         [peak_pow,idx_max_peak]=max(wvfm);
% %         samples_above_threshold_left=find(wvfm>peak_pow*cnf_p.wvfm_portion_selec_l_thres/100);
% %         samples_above_threshold_right=find(wvfm>peak_pow*cnf_p.wvfm_portion_selec_r_thres/100);
% %         start_sample=samples_above_threshold_left(find(samples_above_threshold_left<idx_max_peak(1),1,'last'));
% %         stop_sample=samples_above_threshold_right(find(samples_above_threshold_right>idx_max_peak(1),1,'first'));
%         samples_above_threshold_left=find(wvfm>peak_pow*cnf_p.wvfm_portion_selec_l_thres/100);
%         samples_above_threshold_right=find(wvfm>peak_pow*cnf_p.wvfm_portion_selec_r_thres/100);
%         start_sample=samples_above_threshold_left(find(samples_above_threshold_left<idx_max_peak(1),1,'last'));
%         stop_sample=samples_above_threshold_right(find(samples_above_threshold_right>idx_max_peak(1),1,'first'));
    case 'peak_valley'
        %take the maximum peak and the first left and right valleys 
        [pks,loc_pks,w_pks,p_pks]=findpeaks(wvfm); %peaks
        [vll,loc_vll]=findpeaks(max(wvfm)-wvfm); %valleys
        [~,idx_max_peak]=max(pks);
        idx_int=find(loc_vll<loc_pks(idx_max_peak),1,'last');
        if ~isempty(idx_int)
            start_sample=max(loc_vll(idx_int)-cnf_p.wvfm_portion_selec_l_samples,1);
        else
            start_sample=max(loc_vll(idx_max_peak)-cnf_p.wvfm_portion_selec_l_samples,1);
        end
        idx_int=find(loc_vll>loc_pks(idx_max_peak),1,'first');
        if ~isempty(idx_int)
            stop_sample=min(loc_vll(idx_int)+cnf_p.wvfm_portion_selec_r_samples,N_samples);
        else
            stop_sample=min(loc_vll(idx_max_peak)+cnf_p.wvfm_portion_selec_r_samples,N_samples);
        end
%         start_sample=max(loc_vll(find(loc_vll<loc_pks(idx_max_peak),1,'last'))-cnf_p.wvfm_portion_selec_l_samples,1);
%         stop_sample=min(loc_vll(find(loc_vll>loc_pks(idx_max_peak),1,'first'))+cnf_p.wvfm_portion_selec_r_samples,N_samples);

    case 'rmc_range'  
        start_sample=1;
        stop_sample=0.5*N_sample;
                
    otherwise
        start_sample=1;
        stop_sample=N_samples;
end

% %checking purposes
% if (mod(m,50)==0)
%     f=figure;
%     plot(1:N_samples,wvfm./max(wvfm),'-b');
%     hold on;
%     plot(start_sample:stop_sample,wvfm(start_sample:stop_sample)./max(wvfm),'*r');
%     xlabel('Samples'); ylabel('Norm. Amplitude');
%     title(strcat('Waveform portion selection: #',num2str(m)));
%     legend
%     print('-dpng ',[path_Results,'plots/portion_selec_waveforms/',L1B_filename(17:47),'_wvfm_',num2str(m),'.png']);
%     close(f)
% end


end %end function


% %% --------------- Threshold retracker isolation --------------------------
% %Take the maximum peak in waveform: isolate the samples left and right to
% %the peak above a very low percentage of the peak
% [peak_pow,idx_max_peak]=max(wvfm);
% samples_above_threshold_left=find(wvfm>peak_pow*cnf_p.wvfm_portion_selec_l_thres/100);
% samples_above_threshold_right=find(wvfm<=peak_pow*cnf_p.wvfm_portion_selec_r_thres/100);
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

