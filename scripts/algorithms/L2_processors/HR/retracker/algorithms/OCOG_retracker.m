function [fit_res]=OCOG_retracker(data,cnf_p,chd_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs an OCOG retracker based on Technical note CryoVal-LI
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
% - 
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0:
% V1.1: 
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
    
    
    if cnf_p.wvfm_portion_selec
        %--------------------------------------------------------------------------
        %---------------------- waveform portion selection ------------------------
        %--------------------------------------------------------------------------
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
        n1=start_sample;
        n2=stop_sample;        
    else
        n1=ceil(cnf_p.OCOG_retracker.n1); %first sample to be used
        n2=ceil(cnf_p.OCOG_retracker.n2); %last sample to be used
        flag_nadir_return_in_win=NaN;
    end
    

    
    %----------------- Discarding samples beginning and end -----------------
    if cnf_p.wvfm_discard_samples
        n1=max(n1,1+cnf_p.wvfm_discard_samples_begin);
        n2=min(n2,data.N_samples-cnf_p.wvfm_discard_samples_end);
    end

        % added EM: 20.09.2016
    %checking waveform valid due to L1B processing related to noise
    %filtering outer beams
    if ~any(~isnan(data.HRM.power_wav_filtered(m,n1:n2))) || ~any(data.HRM.power_wav_filtered(m,n1:n2))
        fit_res.Epoch(m)    =   NaN;
        fit_res.Pu(m)       =   NaN;        
        fit_res.COG(m)      =   NaN;
        fit_res.A(m)        =   NaN;
        fit_res.W(m)        =   NaN;
        fit_res.flag_validity_L1B(m) = 0; %due to L1B processing
        fit_res.flag_nadir_return_in_win(m)=0;
        continue;
    end
    
    offset=cnf_p.OCOG_retracker.offset; %offset
    vec_range=n1:n2;
   
    switch cnf_p.OCOG_retracker.param_comp_method
        case 0
            % using the squares of the power waveform a la Frappart
            %--------------------- COG ---------------------------------------------
            fit_res.COG(m)=(sum(vec_range.*(data.HRM.power_wav_filtered(m,n1:n2)).^2)/sum((data.HRM.power_wav_filtered(m,n1:n2)).^2));
            
            %-------------- Width --------------------------------------------------
            fit_res.W(m)= (sum((data.HRM.power_wav_filtered(m,n1:n2)).^2)).^2/sum((data.HRM.power_wav_filtered(m,n1:n2)).^4);
        case 1
             % using the power waveform a la Wingham definition
            %--------------------- COG ---------------------------------------------
            fit_res.COG(m)=(sum(vec_range.*(data.HRM.power_wav_filtered(m,n1:n2)))/sum((data.HRM.power_wav_filtered(m,n1:n2))));
            
            %-------------- Width --------------------------------------------------
            fit_res.W(m)= (sum(data.HRM.power_wav_filtered(m,n1:n2))).^2/sum((data.HRM.power_wav_filtered(m,n1:n2)).^2);
    end
   
   
   %----------- Computing the epoch ---------------------------------------
   switch cnf_p.OCOG_retracker.implementation_method
       case 0
           %----------- thresholding the amplitude --------------------------------
           % amplitude computed as Frappart: from n1 to n2 as for COG
           switch cnf_p.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((data.HRM.power_wav_filtered(m,n1:n2)).^4)/sum((data.HRM.power_wav_filtered(m,n1:n2)).^2));
               case 1
                   fit_res.A(m) = (sum((data.HRM.power_wav_filtered(m,n1:n2)).^2)/sum((data.HRM.power_wav_filtered(m,n1:n2))));
           end
           threshold_pow=cnf_p.OCOG_retracker.percentage_pow_OCOG/100.0*fit_res.A(m);
           dumm = (n1-1)+find(data.HRM.power_wav_filtered(m,n1:n2) > threshold_pow,1);
           if ~isempty(dumm)
               %position by simple linear interpolation
               if dumm>1
                   idx_leading_edge_ISD=offset+dumm-1+((threshold_pow-data.HRM.power_wav_filtered(m,dumm-1))/(data.HRM.power_wav_filtered(m,dumm)-data.HRM.power_wav_filtered(m,dumm-1)));
                   if idx_leading_edge_ISD<0 || idx_leading_edge_ISD>data.N_samples
                       idx_leading_edge_ISD=offset+dumm; %it could happen that power at dumm and dumm-1 is almost the same (linear interpolation does not work properly)
                       %leading to a really high value of epoch.
                   end
               else
                   idx_leading_edge_ISD=offset+dumm;
               end
           else
               %if there is no leading edge or the waveform has displaced that
               %much from the window to the left select the peak as leading
               %edge
               [~,idx_max_peak]=max(data.HRM.power_wav_filtered(m,:),[],2);
               idx_leading_edge_ISD=idx_max_peak;
           end          
       case 1
           %----------- thresholding the amplitude --------------------------------
           % amplitude computed within the window centered ot OCOG
           n1_bis=round(fit_res.COG(m)-fit_res.W(m)/2);
           n2_bis=round(fit_res.COG(m)+fit_res.W(m)/2);
           switch cnf_p.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((data.HRM.power_wav_filtered(m,n1_bis:n2_bis)).^4)/sum((data.HRM.power_wav_filtered(m,n1_bis:n2_bis)).^2));
               case 1
                   fit_res.A(m) = (sum((data.HRM.power_wav_filtered(m,n1_bis:n2_bis)).^2)/sum((data.HRM.power_wav_filtered(m,n1_bis:n2_bis))));
           end
           threshold_pow=cnf_p.OCOG_retracker.percentage_pow_OCOG/100.0*fit_res.A(m);
           dumm = n1_bis-1+find(data.HRM.power_wav_filtered(m,n1_bis:n2_bis) > threshold_pow,1);
           if ~isempty(dumm)
               %position by simple linear interpolation
               if dumm>1
                   idx_leading_edge_ISD=offset+dumm-1+((threshold_pow-data.HRM.power_wav_filtered(m,dumm-1))/(data.HRM.power_wav_filtered(m,dumm)-data.HRM.power_wav_filtered(m,dumm-1)));
                   if idx_leading_edge_ISD<0 || idx_leading_edge_ISD>data.N_samples
                       idx_leading_edge_ISD=offset+dumm; %it could happen that power at dumm and dumm-1 is almost the same (linear interpolation does not work properly)
                       %leading to a really high value of epoch.
                   end
               else
                   idx_leading_edge_ISD=offset+dumm;
               end
           else
               %if there is no leading edge or the waveform has displaced that
               %much from the window to the left select the peak as leading
               %edge
               [~,idx_max_peak]=max(data.HRM.power_wav_filtered(m,:),[],2);
               idx_leading_edge_ISD=idx_max_peak;
           end

       case 2
           %----------- Using OCOG only --------------------------------
           %epoch=offset+(COG-W/2)
           switch cnf_p.OCOG_retracker.param_comp_method
               case 0
                   fit_res.A(m) = sqrt(sum((data.HRM.power_wav_filtered(m,n1:n2)).^4)/sum((data.HRM.power_wav_filtered(m,n1:n2)).^2));
               case 1
                   fit_res.A(m) = (sum((data.HRM.power_wav_filtered(m,n1:n2)).^2)/sum((data.HRM.power_wav_filtered(m,n1:n2))));
           end
           idx_leading_edge_ISD=offset+(fit_res.COG(m)-fit_res.W(m)/2);
   end
   fit_res.Epoch(m)    =   idx_leading_edge_ISD + cnf_p.IFmask_N*cnf_p.ZP - 1;
   fit_res.Pu(m)       =   10*log10(fit_res.A(m)); %dB using the amplitude of the OCOG to be discussed
   
   if cnf_p.plot_fits_flag && ((mod(m,cnf_p.plot_fits_downsampling)==0) || (m==1)) && (data.GEO.LAT(m)>=cnf_p.plot_fits_lat_range(1) && data.GEO.LAT(m)<=cnf_p.plot_fits_lat_range(2))
       f1=figure;
       p1=plot(1:data.N_samples,data.HRM.power_wav_filtered(m,:)./max(data.HRM.power_wav_filtered(m,:)),'-b'); %whole waveform
       hold on;
       p2=plot(n1:n2,data.HRM.power_wav_filtered(m,n1:n2)./max(data.HRM.power_wav_filtered(m,:)),'+r'); %portion selection
       hold on;
       %COG
       plot(fit_res.COG(m),fit_res.A(m)/2./max(data.HRM.power_wav_filtered(m,:)),'*k','MarkerSize',9);
       %window COG
       plot([fit_res.COG(m)-fit_res.W(m)/2,fit_res.COG(m)-fit_res.W(m)/2,fit_res.COG(m)+fit_res.W(m)/2,fit_res.COG(m)+fit_res.W(m)/2],...
           [0,fit_res.A(m)./max(data.HRM.power_wav_filtered(m,:)),fit_res.A(m)./max(data.HRM.power_wav_filtered(m,:)),0],'-g');
       %plot the epoch point
       plot([fit_res.Epoch(m)+1,fit_res.Epoch(m)+1],[0,1],'--m');
       h_leg=legend([p1,p2],'L1B-Waveform', 'L1B-Waveform portion','Location','northeastoutside');
       pos_leg=get(h_leg,'Position');
       annotation('textbox', [pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],...
           'String',{['Epoch = ', num2str(fit_res.Epoch(m),4), ' [r.b]'],...
           ['Pu = ', num2str(fit_res.Pu(m),4), ' [ dB ]']},...
           'FitBoxToText','on');
       grid on
       xlabel('range bin'); ylabel('Norm. amplitude');
        title(['wav # = ', num2str(m), ' (LAT: ',num2str(data.GEO.LAT(m)),' [deg])']);
       print('-dpng ',[path_Results,'plots',filesep,'fitted_waveforms',filesep,L1B_filename(17:47),'_OCOG_wvfm_',num2str(m),'.png']);
       close(f1)
   end
   
   
   
   fit_res.A(m)=10*log10(fit_res.A(m));
   fit_res.flag_validity_L1B(m) = 1;
   fit_res.flag_nadir_return_in_win(m)=flag_nadir_return_in_win;
   
   

   
end

end

