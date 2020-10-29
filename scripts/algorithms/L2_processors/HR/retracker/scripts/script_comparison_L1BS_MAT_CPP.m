%%% -------------- Comparison Stacks L1BS s-6: MAT and CPP ---------------------
% Heavy computational load (this script can block your computer)
common_path = '\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Integration_test\L1BS\';
filename_L1BS_CPP = strcat(common_path,'inputs\S6A_OS20_P4__HR__1S_20170305T065123_20170305T065233_0001.nc');
filename_L1BS_MAT = strcat(common_path,'inputs\S6A_OS20_P4__RAW_1BS_00000000T000000_99999999T999999_0001.mat');
resultPath        = strcat(common_path,'results\');
mkdir(resultPath);

%FILTER THE SURFACES TO BE CONSIDERED IN THE AVERAGED STACKS COMPARISON
N_beams_threshold = 469;

% DON'T PLOT THE STACK COMPARISON FOR SURFACES ABOVE 150 (OTHERWISE TOO BIG)
surf_index_limit = 150;

%FOR QUANTIZATION EMULATION FROM MATLAB
quantization = 8; % [8 or 16]

figure_format='jpg';
res_fig='-r150';

switch lower(figure_format)
    case 'pdf'
        file_ext='.pdf';
        print_file='-dpdf';
    case 'eps'
        file_ext='.eps';
        print_file='-depsc';
    case 'png'
        file_ext='.png';
        print_file='-dpng';
    case 'jpg'
        file_ext='.jpg';
        print_file='-djpeg';        
end

%% load MAT L1Bs
load(filename_L1BS_MAT,'beams_rng_cmpr','beams_rng_cmpr_I',...
                        'beams_rng_cmpr_Q','N_beams_stack','N_beams_multilooking');
beams_rng_cmpr_MAT=beams_rng_cmpr;
beams_rng_cmpr_I_MAT=beams_rng_cmpr_I;
beams_rng_cmpr_Q_MAT=beams_rng_cmpr_Q;
s=size(permute(beams_rng_cmpr_I_MAT,[3 2 1]));
clear beams_rng_cmpr beams_rng_cmpr_I beams_rng_cmpr_Q;


%% Load L1BS CPP:
%L1B-S
look_index_20_ku=double(ncread(filename_L1BS_CPP,'data_20/ku/look_index')).';
iq_scale_factor = double(ncread(filename_L1BS_CPP,'data_20/ku/iq_scale_factor'));
look_i_samples_20_ku=double(ncread(filename_L1BS_CPP,'data_20/ku/look_i_samples'));
look_q_samples_20_ku=double(ncread(filename_L1BS_CPP,'data_20/ku/look_q_samples'));

N_samples_range = s(1);
N_beams_max = max(N_beams_stack);
N_surf_total = s(3);
wfm_geo_corr_CPP = zeros(N_surf_total,max(N_beams_stack),N_samples_range);
for i_surf=1:N_surf_total
    for i_beam=1:N_beams_max
        wfm_geo_corr_CPP(i_surf,i_beam,:)=(look_i_samples_20_ku(:,i_beam,i_surf)+...
            1i*look_q_samples_20_ku(:,i_beam,i_surf)).*iq_scale_factor(i_beam,i_surf);
    end
end
%clear look_i_samples_20_ku;
%clear look_q_samples_20_ku;
%clear iq_scale_factor;

beams_rng_cmpr_CPP=abs(wfm_geo_corr_CPP).^2;
clear wfm_geo_corr_CPP;

%% Quantization of the stack for Matlab
s=size(beams_rng_cmpr_I_MAT);
N_beams_max = s(2);
N_samples_range=s(3);


if quantization==8
    
    lookISamples=int8(NaN(N_surf_total,N_beams_max,N_samples_range));
    lookQSamples=int8(NaN(N_surf_total,N_beams_max,N_samples_range));
    iqScaleFactor = single(zeros(N_surf_total,N_beams_max));
    for i_surf=1:N_surf_total
        for i_beam=1:N_beams_max
            iqScaleFactor(i_surf,i_beam) = single(max([max(abs(beams_rng_cmpr_I_MAT(i_surf,i_beam,:))),...
                max(abs(beams_rng_cmpr_Q_MAT(i_surf,i_beam,:)))])/(2^7-1));
            lookISamples(i_surf,i_beam,:)=int8(round(squeeze(beams_rng_cmpr_I_MAT(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
            lookQSamples(i_surf,i_beam,:)=int8(round(squeeze(beams_rng_cmpr_Q_MAT(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
        end
    end
    
elseif quantization==16
    
    lookISamples=int16(NaN(N_surf_total,N_beams_max,N_samples_range));
    lookQSamples=int16(NaN(N_surf_total,N_beams_max,N_samples_range));
    iqScaleFactor = single(zeros(N_surf_total,N_beams_max));
    for i_surf=1:N_surf_total
        for i_beam=1:N_beams_max
            iqScaleFactor(i_surf,i_beam) = single(max([max(abs(beams_rng_cmpr_I_MAT(i_surf,i_beam,:))),...
                max(abs(beams_rng_cmpr_Q_MAT(i_surf,i_beam,:)))])/(2^15-1));
            lookISamples(i_surf,i_beam,:)=int16(round(squeeze(beams_rng_cmpr_I_MAT(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
            lookQSamples(i_surf,i_beam,:)=int16(round(squeeze(beams_rng_cmpr_Q_MAT(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
        end
    end
end
clear beams_rng_cmpr_I_MAT beams_rng_cmpr_Q_MAT;

% -------------------- Decoding -------------------------------------------
beams_rng_cmpr_MAT_QUAN=zeros(N_surf_total,N_beams_max,N_samples_range);

for i_surf=1:N_surf_total
    for i_beam=1:N_beams_max
        beams_rng_cmpr_MAT_QUAN(i_surf,i_beam,:) = (double(lookISamples(i_surf,i_beam,:))+...
            1i*double(lookQSamples(i_surf,i_beam,:))).*double(iqScaleFactor(i_surf,i_beam));
    end
end

beams_rng_cmpr_MAT_QUAN = (abs(beams_rng_cmpr_MAT_QUAN)).^2;


%% Comparison 
max_image=10*log10(max([max(beams_rng_cmpr_MAT(:)),max(beams_rng_cmpr_CPP(:)),...
                   max(beams_rng_cmpr_MAT_QUAN(:))]));
min_image=max_image-40.0;

max_image_error_abs=max([max(abs(beams_rng_cmpr_MAT(:)-beams_rng_cmpr_CPP(:))),...
                         max(abs(beams_rng_cmpr_MAT_QUAN(:)-beams_rng_cmpr_CPP(:)))]);
min_image_error_abs=0;

max_image_error_rel=100;
min_image_error_rel=0;

VISIBLE=0;
set_default_plot;
set(0,'defaultFigureVisible','on');
colormapATDD ='jet';

stack_MAT=NaN(N_surf_total,N_beams_max,N_samples_range);
stack_MAT_QUAN=NaN(N_surf_total,N_beams_max,N_samples_range);
stack_CPP=NaN(N_surf_total,N_beams_max,N_samples_range);

for i_surf=1:N_surf_total
    
    
    stack_MAT(i_surf,1:N_beams_stack(i_surf),:)=squeeze(beams_rng_cmpr_MAT(i_surf,1:N_beams_stack(i_surf),:));
    stack_MAT_QUAN(i_surf,1:N_beams_stack(i_surf),:)=squeeze(beams_rng_cmpr_MAT_QUAN(i_surf,1:N_beams_stack(i_surf),:));
    stack_CPP(i_surf,1:N_beams_stack(i_surf),:)=squeeze(beams_rng_cmpr_CPP(i_surf,1:N_beams_stack(i_surf),:));
    
    if i_surf<=surf_index_limit
        f1=figure;
        h(1)=subplot(3,3,1);
        surfc(10*log10(squeeze(stack_MAT(i_surf,:,:))),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('MAT (no QUAN.)'));
        
        h(2)=subplot(3,3,2);
        surfc(10*log10(squeeze(stack_MAT_QUAN(i_surf,:,:))),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('MAT (QUAN.',num2str(quantization),')'));
        
        h(3)=subplot(3,3,3);
        surfc(10*log10(squeeze(stack_CPP(i_surf,:,:))),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('GPP '));
        
        
        h(4)=subplot(3,3,4);
        surfc(squeeze(abs(stack_MAT(i_surf,:,:)-stack_CPP(i_surf,:,:))),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image_error_abs,max_image_error_abs]); ylabel(c,'[W]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('Absolute error MAT (no. QUAN.) & GPP'));
        
        h(5)=subplot(3,3,5);
        surfc(squeeze(abs((stack_MAT(i_surf,:,:)-stack_CPP(i_surf,:,:))./stack_MAT(i_surf,:,:)).*100),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image_error_rel,max_image_error_rel]); ylabel(c,'[%]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('Relative error MAT (no. QUAN.) & GPP'));
                
        h(6)=subplot(3,3,7);
        surfc(squeeze(abs(stack_MAT_QUAN(i_surf,:,:)-stack_CPP(i_surf,:,:))),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image_error_abs,max_image_error_abs]); ylabel(c,'[W]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('Absolute error MAT (QUAN. ',num2str(quantization),') & GPP'));
        
        h(7)=subplot(3,3,8);
        surfc(squeeze(abs((stack_MAT_QUAN(i_surf,:,:)-stack_CPP(i_surf,:,:))./stack_MAT_QUAN(i_surf,:,:)).*100),'FaceColor','texturemap','EdgeColor','none');
        view(2); grid off;
        colormap(colormapATDD); c=colorbar; caxis([min_image_error_rel,max_image_error_rel]); ylabel(c,'[%]');
        xlabel('Samples'); ylabel('Beams');
        title(strcat('Relative error MAT (QUAN. ',num2str(quantization),') & GPP'));
        
        pos = get(h,'Position');
        x_pos_middle_1 = mean(cellfun(@(v)v(1),pos(1:2)));
        set(h(4),'Position',[x_pos_middle_1,pos{4}(2:end)])
        x_pos_middle_2 = mean(cellfun(@(v)v(1),pos(2:3)));
        set(h(5),'Position',[x_pos_middle_2,pos{5}(2:end)])
        set(h(6),'Position',[x_pos_middle_1,pos{6}(2:end)])
        set(h(7),'Position',[x_pos_middle_2,pos{7}(2:end)])
        
        [axT,hT]=suplabel(strcat('Comparison MAT & GPP #',num2str(i_surf)),'t');
        
        print(print_file,res_fig,strcat(resultPath,'Stack_comparison_MAT_CPP_',num2str(i_surf,'%05d'),file_ext));
        close(f1);
    end
end

clear beams_rng_cmpr_MAT beams_rng_cmpr_MAT_QUAN beams_rng_cmpr_CPP

%% Comparison of the average stacks 
% Filter out the non-complete stacks (below num beams)
idx_stacks_int=N_beams_multilooking>=N_beams_threshold;

stack_MAT=stack_MAT(idx_stacks_int,:,:);
stack_MAT_QUAN=stack_MAT_QUAN(idx_stacks_int,:,:);
stack_CPP=stack_CPP(idx_stacks_int,:,:);

stack_MAT_AVG=nanmean(stack_MAT,1);
stack_MAT_QUAN_AVG=nanmean(stack_MAT_QUAN,1);
stack_CPP_AVG=nanmean(stack_CPP,1);

clear stack_MAT stack_MAT_QUAN stack_CPP;

% max_image=10*log10(max([max(stack_MAT_AVG(:)),max(stack_MAT_QUAN_AVG(:)),...
%                    max(stack_CPP_AVG(:))]));
% min_image=max_image-40.0;

max_image_error_abs=max([max(abs(stack_MAT_AVG(:)-stack_CPP_AVG(:))),...
                         max(abs(stack_MAT_QUAN_AVG(:)-stack_CPP_AVG(:)))]);
min_image_error_abs=0;

max_image_error_rel=100;
min_image_error_rel=0;

clear h pos;

f1=figure;
h(1)=subplot(3,3,1);
surfc(10*log10(squeeze(stack_MAT_AVG)),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
xlabel('Samples'); ylabel('Beams');
title(strcat('MAT (no QUAN.)'));

h(2)=subplot(3,3,2);
surfc(10*log10(squeeze(stack_MAT_QUAN_AVG)),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
xlabel('Samples'); ylabel('Beams');
title(strcat('MAT (QUAN.',num2str(quantization),')'));

h(3)=subplot(3,3,3);
surfc(10*log10(squeeze(stack_CPP_AVG)),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
xlabel('Samples'); ylabel('Beams');
title(strcat('GPP'));

h(4)=subplot(3,3,4);
surfc(squeeze(abs(stack_MAT_AVG-stack_CPP_AVG)),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image_error_abs,max_image_error_abs]); ylabel(c,'[W]');
xlabel('Samples'); ylabel('Beams');
title(strcat('Absolute error MAT (no. QUAN.) & GPP'));

h(5)=subplot(3,3,5);
surfc(squeeze(abs((stack_MAT_AVG-stack_CPP_AVG)./stack_MAT_AVG).*100),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image_error_rel,max_image_error_rel]); ylabel(c,'[%]');
xlabel('Samples'); ylabel('Beams');
title(strcat('Relative error MAT (no. QUAN.) & GPP'));


h(6)=subplot(3,3,7);
surfc(squeeze(abs(stack_MAT_QUAN_AVG-stack_CPP_AVG)),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image_error_abs,max_image_error_abs]); ylabel(c,'[W]');
xlabel('Samples'); ylabel('Beams');
title(strcat('Absolute error MAT (QUAN. ',num2str(quantization),') & GPP'));


h(7)=subplot(3,3,8);
surfc(squeeze(abs((stack_MAT_QUAN_AVG-stack_CPP_AVG)./stack_MAT_QUAN_AVG).*100),'FaceColor','texturemap','EdgeColor','none');
view(2); grid off;
colormap(colormapATDD); c=colorbar; caxis([min_image_error_rel,max_image_error_rel]); ylabel(c,'[%]');
xlabel('Samples'); ylabel('Beams');
title(strcat('Relative error MAT (QUAN. ',num2str(quantization),') & GPP'));

pos = get(h,'Position');
x_pos_middle_1 = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(4),'Position',[x_pos_middle_1,pos{4}(2:end)])
x_pos_middle_2 = mean(cellfun(@(v)v(1),pos(2:3)));
set(h(5),'Position',[x_pos_middle_2,pos{5}(2:end)])
set(h(6),'Position',[x_pos_middle_1,pos{6}(2:end)])
set(h(7),'Position',[x_pos_middle_2,pos{7}(2:end)])

[axT,hT]=suplabel('Comparison AVERAGED STACKS MAT & GPP','t');

print(print_file,res_fig,strcat(resultPath,'Stack_AVERAGED_comparison_MAT_CPP',file_ext));
close(f1);


% %% Comparison stacks MAT with quantization
% s=size(beams_rng_cmpr_I);
% N_beams_max = s(2);
% N_samples_range=s(3);
% quantization = 16;
% 
% if quantization==8
%     
%     lookISamples=int8(NaN(N_total_surf_loc,N_beams_max,N_samples_range));
%     lookQSamples=int8(NaN(N_total_surf_loc,N_beams_max,N_samples_range));
%     iqScaleFactor = single(zeros(N_total_surf_loc,N_beams_max));
%     for i_surf=1:N_total_surf_loc
%         for i_beam=1:N_beams_max
%             iqScaleFactor(i_surf,i_beam) = single(max([max(abs(beams_rng_cmpr_I(i_surf,i_beam,:))),...
%                 max(abs(beams_rng_cmpr_Q(i_surf,i_beam,:)))])/(2^7-1));
%             lookISamples(i_surf,i_beam,:)=int8(round(squeeze(beams_rng_cmpr_I(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
%             lookQSamples(i_surf,i_beam,:)=int8(round(squeeze(beams_rng_cmpr_Q(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
%         end
%     end
%     
% elseif quantization==16
%     
%     lookISamples=int16(NaN(N_total_surf_loc,N_beams_max,N_samples_range));
%     lookQSamples=int16(NaN(N_total_surf_loc,N_beams_max,N_samples_range));
%     iqScaleFactor = single(zeros(N_total_surf_loc,N_beams_max));
%     for i_surf=1:N_total_surf_loc
%         for i_beam=1:N_beams_max
%             iqScaleFactor(i_surf,i_beam) = single(max([max(abs(beams_rng_cmpr_I(i_surf,i_beam,:))),...
%                 max(abs(beams_rng_cmpr_Q(i_surf,i_beam,:)))])/(2^15-1));
%             lookISamples(i_surf,i_beam,:)=int16(round(squeeze(beams_rng_cmpr_I(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
%             lookQSamples(i_surf,i_beam,:)=int16(round(squeeze(beams_rng_cmpr_Q(i_surf,i_beam,:)./iqScaleFactor(i_surf,i_beam))));
%         end
%     end
% end
% 
% % -------------------- Decoding -------------------------------------------
% beams_rng_cmpr_decoded=zeros(N_total_surf_loc,N_beams_max,N_samples_range);
% 
% for i_surf=1:N_total_surf_loc
%     for i_beam=1:N_beams_max
%         beams_rng_cmpr_decoded(i_surf,i_beam,:) = (double(lookISamples(i_surf,i_beam,:))+...
%             1i*double(lookQSamples(i_surf,i_beam,:))).*double(iqScaleFactor(i_surf,i_beam));
%     end
% end
% 
% % ------------------- Ploting ---------------------------------------------
% beams_rng_cmpr_decoded = (abs(beams_rng_cmpr_decoded)).^2;
% 
% max_image=10*log10(max([max(beams_rng_cmpr(:)),max(beams_rng_cmpr_decoded(:))]));
% min_image=max_image-40;
% 
% max_image_error_abs=max(abs(beams_rng_cmpr(:)-beams_rng_cmpr_decoded(:)));
% min_image_error_abs=0;
% 
% max_image_error_rel=20;
% min_image_error_rel=0;
% 
% 
% for i_surf=1:N_total_surf_loc
%     
%     stack_MAT=NaN(N_beams_max,N_samples_range);
%     stack_MAT_decoded=NaN(N_beams_max,N_samples_range);
%     
%     stack_MAT(1:N_beams_stack(i_surf),:)=squeeze(beams_rng_cmpr(i_surf,1:N_beams_stack(i_surf),:));
%     stack_MAT_decoded(1:N_beams_stack(i_surf),:)=squeeze(beams_rng_cmpr_decoded(i_surf,1:N_beams_stack(i_surf),:));
%     
%     f1=figure;
%     subplot(2,2,1);
%     surfc(10*log10(stack_MAT),'FaceColor','texturemap','EdgeColor','none');
%     view(2); grid off;
%     colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
%         xlabel('Samples'); ylabel('Beams');
%     title(strcat('MAT #',num2str(i_surf)));
%     subplot(2,2,2);
%     surfc(10*log10(stack_MAT_decoded),'FaceColor','texturemap','EdgeColor','none');
%     view(2); grid off;
%     colormap(colormapATDD); c=colorbar; caxis([min_image,max_image]); ylabel(c,'[dBW]');
%     xlabel('Samples'); ylabel('Beams');
%     title(strcat('MAT quant=',num2str(quantization),' #',num2str(i_surf)));
%     subplot(2,2,3);
%     surfc(abs(stack_MAT-stack_MAT_decoded),'FaceColor','texturemap','EdgeColor','none');
%     view(2); grid off;
%     colormap(colormapATDD); c=colorbar; caxis([min_image_error_abs,max_image_error_abs]); ylabel(c,'[W]');
%     xlabel('Samples'); ylabel('Beams');
%     title(strcat('Absolute error #',num2str(i_surf)));
%     subplot(2,2,4);
%     surfc(abs((stack_MAT-stack_MAT_decoded)./stack_MAT).*100,'FaceColor','texturemap','EdgeColor','none');
%     view(2); grid off;
%     colormap(colormapATDD); c=colorbar; caxis([min_image_error_rel,max_image_error_rel]); ylabel(c,'[%]');
%     xlabel('Samples'); ylabel('Beams');
%     title(strcat('Relative error #',num2str(i_surf)));
%     [axT,hT]=suplabel('Comparison MAT & MAT quantization','t');
%     print(print_file,res_fig,strcat(resultPath,'Stack_comparison_MAT_quantization_',num2str(quantization),'_',num2str(i_surf),file_ext)); 
%     close(f1)
% end


