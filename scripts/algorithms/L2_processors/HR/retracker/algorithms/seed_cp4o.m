function data = seed_cp4o(data,cnf_p,chd_p,cst_p)
%seed_cp4o Production of seed for S3 retracking
%   The reference gate is found for the retracking solution
%   1. track sections are defined
%   2. window delay (WD) along each section is smoothed
%   3. seed is found accounting for the WD residuals
%
% This version is for SCOOP project --> S3B in closed loop over Cuba
%

tic
%% aux files data extraction
c=cst_p.c_cst; % light speed in m/s
BW=chd_p.bw_rx_ku_chd; % BW in Hz

for i=1:1 % fake loop for hiding old code section
    % %% create raw ocean sections (only 100m off-shore limitation)
    % % surf type = 0 --> ocean
    % % surf type = 3 --> land
    % % track_section: rows = ocean track section, columns = start&end of section
    %
    % flag_jumps=diff(data.surf_type_flag);
    %
    % % setting sections limits: track starts in OCEAN
    % if data.surf_type_flag(1)==0 % track starts in OCEAN
    %     section_id=1; track_section(section_id,1)=1; % track_section(section_id,recordlim1:recordlim2)
    % else
    %     section_id=0;
    % end
    %
    % % setting sections limits: intermediate sections
    % for i_wfm=1:length(flag_jumps) % waveforms
    %     if flag_jumps(i_wfm)==3    % from sea to land
    %         track_section(section_id,2)=i_wfm;
    %     elseif flag_jumps(i_wfm)==-3 % from land to sea
    %         section_id=section_id+1;
    %         track_section(section_id,1)=i_wfm;
    %     end
    % end
    %
    % % setting sections limits: track ends in OCEAN
    % if data.surf_type_flag(end)==0 % track ends in  OCEAN
    %     track_section(section_id,2)=length(data.surf_type_flag);
    % end
    %
    % % delete sea sections of less than 5 records
    % if size(track_section,1) > 1 % if more than 1 section
    %     i_sec=1;
    %     while i_sec < size(track_section,1) % loop through sections
    %         if track_section(i_sec,2)-track_section(i_sec,1) < 5 % if track section lenght < 5
    %             track_section(i_sec,:)=[]; % delete current section
    %         else
    %             i_sec = i_sec + 1;
    %         end
    %     end
    % end
    %
    % % delete land sections of less than 5 records
    % if size(track_section,1) > 1 % if more than 1 section
    %     i_sec=1;
    %     while i_sec < size(track_section,1)-1 % loop through sections
    %         if track_section(i_sec+1,1)-track_section(i_sec,2) < 5 % if (next section start - current section end) lenght < 5
    %             track_section(i_sec,2) = track_section(i_sec+1,2); % merge this section with the next
    %             track_section(i_sec+1,:)=[]; % delete next section
    %         else
    %             i_sec = i_sec + 1;
    %         end
    %     end
    % end
end

% ***old conditions***: --- ONLY CLOSED LOOP PRODUCTS STARTING&ENDING IN OCEAN ARE VALID
% mode_id: [0,1,2] --> 'closed_loop, open_loop, open_loop_fixed_gain'
% surf_type_flag: [0,1,2,3] --> 'open_ocean, sea_ice, lead, unclassified'
% if ~(data.surf_type_flag(1)==0 && data.surf_type_flag(end)==0)
%     disp('------------ IT IS NOT POSIBLE TO RUN THIS COASTAL PROCESSING IN PRODUCTS NOT STARTING FROM OCEAN ------------')
%     return
% end
if ~max(data.GEO.mode_id)==0
%     disp('------------ IT IS NOT POSIBLE TO RUN THIS COASTAL PROCESSING IN PRODUCTS IN OPEN LOOP TRACKING MODE ------------')
    disp('------------ RUNNING THIS COASTAL PROCESSING IN PRODUCT IN OPEN LOOP TRACKING MODE ------------')
%     return
end


%% find window elevation residuals (wrt poynomial fitting) by sections
window_elev=data.GEO.H - data.MEA.win_delay*c/2; % altitude - Windel ... in m
window_elevsamp=(data.GEO.H - data.MEA.win_delay*c/2)/(c/2/BW); % altitude - Windel ... in samples units
[~,m2]=max(data.HRM.power_wav,[],1);% detect max sample in the waveform
window_elevsamp_lepmx=window_elevsamp-m2-2; % LEP altitude - Windel in samples units
window_elev_lepmx=window_elevsamp_lepmx*c/2/BW; % LEP altitude - Windel in m

open_ocean_lim = 1e4; % distance to coast for Open Ocean

bias_MSS_windel=mean(data.COR.mss2_20(data.GEO.dist_coast_20>open_ocean_lim)-window_elev_lepmx(data.GEO.dist_coast_20>open_ocean_lim)); % bias MSS2 to windel, based in Open Ocean differences
MSS2windel=data.COR.mss2_20-bias_MSS_windel; % Elevation of L2 MSS CNES/CLS15 wrt ellipsoid, unbiased wrt windel elevation ... in m

windeljumps_m=window_elev - MSS2windel; % jumps of window delay in meters wrt unbiased MSS2
windeljumps_samp=windeljumps_m/(c/2/BW); % jumps of window delay in samples wrt unbiased MSS2

seed=round(windeljumps_samp)-5; % empirical bias of a number of samples
% if seed is out of bounds, get the maximum position
seed(seed < cnf_p.wvfm_portion_selec_l_samples | seed > 128 - cnf_p.wvfm_portion_selec_r_samples)=m2(seed < cnf_p.wvfm_portion_selec_l_samples | seed > 128 - cnf_p.wvfm_portion_selec_r_samples);

% wfms=data.HRM.power_wav;
% for i=565:572
%     figure; hold all; plot(wfms(:,i),'.-'); plot(round(seed(i)),wfms(round(seed(i)),i),'o','markersize', 12);
% end

for i=1:1 % fake loop for hiding old code section

% % find sections for extrapolation
% 
% % first zone will be from record 1 forward, searching openloop&ocean --> closeloop&ocean
% % once closeloop&ocean zone is over, extrapolation section is closed
% % openloop=1, closeloop=0 // ocean=0, land=3
% extrapzone_forward(1)=0;
% if (data.GEO.mode_id(1)==1 & data.surf_type_flag(1)==0)
%     disp('product starts in OCEAN & OPEN LOOP');
%     extrapzone_forward(1)=1;
%     for i=2:length(data.GEO.mode_id)
%         if (data.GEO.mode_id(i)==0 & data.surf_type_flag(i)==0)
%             extrapzone_forward(2)=i-1;
%             for j=i+1:length(data.GEO.mode_id)
%                 if ~(data.GEO.mode_id(j)==0 & data.surf_type_flag(j)==0)
%                     extrapzone_forward(3)=j-1;
%                     break
%                 end
%             end
%             break
%         end
%     end
% end
% disp(['Extrapolation forward from ' num2str(extrapzone_forward(1)) '-' num2str(extrapzone_forward(2)) ' to ' num2str(extrapzone_forward(2)+1) '-' num2str(extrapzone_forward(3)) '.'])
% 
% % second zone will be from last record backward, same search
% extrapzone_backward(1)=0;
% if (data.GEO.mode_id(end)==1 & data.surf_type_flag(end)==0)
%     disp('product ends in OCEAN & OPEN LOOP');
%     extrapzone_backward(1)=length(data.GEO.mode_id);
%     for i=length(data.GEO.mode_id):-1:1
%         if (data.GEO.mode_id(i)==0 & data.surf_type_flag(i)==0)
%             extrapzone_backward(2)=i+1;
%             for j=i-1:-1:1
%                 if ~(data.GEO.mode_id(j)==0 & data.surf_type_flag(j)==0)
%                     extrapzone_backward(3)=j+1;
%                     break
%                 end
%             end
%             break
%         end
%     end
% end
% disp(['Extrapolation backward from ' num2str(extrapzone_backward(1)) '-' num2str(extrapzone_backward(2)) ' to ' num2str(extrapzone_backward(2)+1) '-' num2str(extrapzone_backward(3)) '.'])
% 
% track_section(1,1)=1; track_section(1,2)=extrapzone_forward(2);
% track_section(2,1)=extrapzone_backward(2); track_section(2,2)=extrapzone_backward(1);
% track_section_to(1,1)=extrapzone_forward(2)+1; track_section_to(1,2)=extrapzone_forward(3);
% track_section_to(2,1)=extrapzone_backward(3); track_section_to(2,2)=extrapzone_backward(2)-1;
% 
% % models to sample units
% elev_L2geoid_samp=data.COR.geoid_20'/(c/2/BW); % Elevation of L2 EGM2008 Geoid (wrt ellipsoid), in samples units
% elev_L2mss01_samp=data.COR.mss1_20'/(c/2/BW);  % Elevation of L2 MSS CLS15 (wrt ellipsoid), in samples units
% elev_L2mss02_samp=data.COR.mss2_20'/(c/2/BW);  % Elevation of L2 MSS CNES/CLS15 (wrt ellipsoid), in samples units
% elev_geoid96_samp=geoidheight(data.GEO.LAT,data.GEO.LON,'EGM96')'/(c/2/BW); % Elevation of EGM96 Geoid (wrt ellipsoid), in samples units
% 
% % polynomial interpolation of each section of order polydeg
% for i_sec=1:2
%     if track_section(i_sec,2)-track_section(i_sec,1) < 250 % section with less than 250 records
%         polydeg(i_sec)=2;
%     elseif track_section(i_sec,2)-track_section(i_sec,1) < 500 % section with less than 500 records
%         polydeg(i_sec)=3;
%     elseif track_section(i_sec,2)-track_section(i_sec,1) < 900 % section with less than 900 records
%         polydeg(i_sec)=4;
%     else % section with more than 1100 records
%         polydeg(i_sec)=5;
%     end
%     
%     % Y-axis
%     elev_WDmaxpw_samp=window_elevsamp_lepmx(track_section(i_sec,1):track_section(i_sec,2))'; % Y axis is windel elevation corrected by maxpower position, in samples units
%     % center & scale the X-axis
%     x_section=track_section(i_sec,1):track_section(i_sec,2); % X axis
%     center=mean(x_section); scale=max(x_section)-min(x_section);
%     x1=(track_section(i_sec,1)-center)/scale; x2=(track_section(i_sec,2)-center)/scale;
%     xstep=1/(track_section(i_sec,2)-track_section(i_sec,1));
%     polyxaxis=x1:xstep:x2;
%     polyxaxis1=(track_section_to(i_sec,1)-center)/scale:xstep:(track_section_to(i_sec,2)-center)/scale; % x axis in zone extrapto
%     % polynomial interpolation of WD
%     not_corrupted_recs=find(elev_WDmaxpw_samp~=0);
%     valpol_section_WD=polyfit(polyxaxis(not_corrupted_recs)',elev_WDmaxpw_samp(not_corrupted_recs),polydeg(i_sec));% extrapolating the polynomy
%     fitsection_WD(track_section_to(i_sec,1):track_section_to(i_sec,2)) = polyval(valpol_section_WD,polyxaxis1); % evaluating the extrapolation
%     
%     if i_sec==1 % first zone
%         % unbiasing models wrt WD elevation on the limit record of zone 1
%         record_ref=track_section_to(i_sec,1);
%     else % second zone
%         % unbiasing models wrt WD elevation on the limit record of zone 2
%         record_ref=track_section_to(i_sec,2);
%     end
%     bias_Geoid96toWD=elev_geoid96_samp(record_ref)-fitsection_WD(record_ref); % for unbiasing the Matlab geoid wrt WD elevation
%     bias_GeoidL2toWD=elev_L2geoid_samp(record_ref)-fitsection_WD(record_ref); % for unbiasing the L2 geoid wrt WD elevation
%     bias_MSS01toWD  =elev_L2mss01_samp(record_ref)-fitsection_WD(record_ref); % for unbiasing the L2 MSS CNESCLS wrt WD elevation
%     bias_MSS02toWD  =elev_L2mss02_samp(record_ref)-fitsection_WD(record_ref); % for unbiasing the L2 MSS DTU wrt WD elevation
%     % matching models with open ocean WD elevation
%     WD_zone=window_elevsamp(track_section_to(i_sec,1):track_section_to(i_sec,2)); % windel elevation in samples units
%     WDmxpw_zone=fitsection_WD(track_section_to(i_sec,1):track_section_to(i_sec,2)); % extrapolation result of windel elevation corrected by maxpower position
%     Geoid96_zone=elev_geoid96_samp(track_section_to(i_sec,1):track_section_to(i_sec,2))-bias_Geoid96toWD; % match of Matlab geoid to WD profile
%     GeoidL2_zone=elev_L2geoid_samp(track_section_to(i_sec,1):track_section_to(i_sec,2))-bias_GeoidL2toWD; % match of L2 geoid to WD profile
%     MSS01_zone=elev_L2mss01_samp(track_section_to(i_sec,1):track_section_to(i_sec,2))-bias_MSS01toWD; % match of L2 MSS CNESCLS to WD profile
%     MSS02_zone=elev_L2mss02_samp(track_section_to(i_sec,1):track_section_to(i_sec,2))-bias_MSS02toWD; % match of L2 MSS DTU to WD profile
%     % computing the jump of WD --> seed  **********************************************************************************************************************
%     LEP_zone_WD(track_section_to(i_sec,1):track_section_to(i_sec,2))        = WD_zone-WDmxpw_zone;
%     LEP_zone_Geoid96(track_section_to(i_sec,1):track_section_to(i_sec,2))   = WD_zone-Geoid96_zone';
%     LEP_zone_GeoidL2(track_section_to(i_sec,1):track_section_to(i_sec,2))   = WD_zone-GeoidL2_zone';
%     LEP_zone_MSS01(track_section_to(i_sec,1):track_section_to(i_sec,2))     = WD_zone-MSS01_zone';
%     LEP_zone_MSS02(track_section_to(i_sec,1):track_section_to(i_sec,2))     = WD_zone-MSS02_zone';
%     % plotting the unbiased models, WD and WD+maxpw
%     figure; hold all; title(['Extrapolation result of zone ' num2str(i_sec) ' of product ' data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0(17:31)]); xlabel('records'); ylabel('windel elevation [range bins]');
%     plot(track_section(i_sec,1):track_section(i_sec,2),window_elevsamp_lepmx(track_section(i_sec,1):track_section(i_sec,2)),'DisplayName','WD&maxpw elev')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),WDmxpw_zone,'DisplayName','extrapolation of WD&maxpw elev')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),Geoid96_zone,'DisplayName','unbiased Geoid96')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),GeoidL2_zone,'DisplayName','unbiased L2 Geoid')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),MSS01_zone,'DisplayName','unbiased L2 MSS01')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),MSS02_zone,'DisplayName','unbiased L2 MSS02')
%     plot(track_section_to(i_sec,1):track_section_to(i_sec,2),window_elevsamp_lepmx(track_section_to(i_sec,1):track_section_to(i_sec,2)),'DisplayName','WD&maxpw elev original')
%     legend('show'); set(gca,'FontSize', 18);
%     saveas(gcf,['Extrapolation_zone' num2str(i_sec) '_prod_' data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0(17:31) '.jpg']);
%     saveas(gcf,['Extrapolation_zone' num2str(i_sec) '_prod_' data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0(17:31) '.fig']);
%     close;
%     clear ('elev_WDmaxpw_samp','valpol_section', 'center', 'scale','x1','x2','xstep','polyxaxis','polyxaxis1','fitsection','x_section');
% end
% 
% % mapping the wfms in one track, correcting by Window Delay Elevation (and LEP)
% WDHsamp=window_elevsamp-max(window_elevsamp); % WindelElev -- in samples units and referenced to maximum WDelev
% m3=m2;
% m3(track_section_to(1,1):track_section_to(1,2))=LEP_zone_MSS02(track_section_to(1,1):track_section_to(1,2));
% m3(track_section_to(2,1):track_section_to(2,2))=LEP_zone_MSS02(track_section_to(2,1):track_section_to(2,2));
% for i_rec=1:length(WDHsamp)
%     wfm_corrWDH(i_rec,(1:128)-round(WDHsamp(i_rec)))=wfms(1:128,i_rec); % radargram corrected by window delay elevation
%     wfm_corrWDHcp4o(i_rec,(1:128)+100-(round(WDHsamp(i_rec)+m3(i_rec)-43)))=wfms(1:128,i_rec); % radargram corrected by window delay elevation + residual of poynomial fitting
%     wfm_corrWDHmaxpw(i_rec,(1:128)+100-(round(WDHsamp(i_rec)+m2(i_rec)-43)))=wfms(1:128,i_rec); % radargram corrected by window delay elevation + location of maximum power in the wfm
% end
% figure; imagesc(wfm_corrWDH'); caxis([0 2e3]); % wfms corrected by WD, elev

end

data.CP4O.seed_Elev=data.GEO.H - (data.MEA.win_delay*BW + seed)*c/2/BW; % elevation in meters of seed position
data.CP4O.seedpw=seed; % seed in window samples
data.CP4O.window_elev_samp=window_elevsamp; % window elevation in samples

save(['C:\Users\Pablo\Desktop\S3MPC\Routine\code\waterlevel\edu\product_' data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0(17:31) '_run_' datestr(now,'ddmmmyyyy_HHMMSS')], 'data');
kmlwritepoint(['C:\Users\Pablo\Desktop\S3MPC\Routine\code\waterlevel\edu\product_' data.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level0(17:31) '_run_' datestr(now,'ddmmmyyyy_HHMMSS')],data.GEO.LAT,data.GEO.LON);
toc

end