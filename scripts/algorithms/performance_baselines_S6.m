function performance_baselines_S6(SSH, SWH, sigma0, COR, filename_L2_ISR, index_inside_mask,...
                                            name_bs, path_comparison_results,cnf_tool,varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for the comparison of the geophsycial retrievals of
% different input baselines of isardSAT:
% exctracting the performance in terms of noise metric (as std or/and RMSE)
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 30/08/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       - SSH, SWH, sigma0, COR = N_baselines x num. records array with retrieved parameters from L2 product
%       - filename_L2_ISR = cell array containg L2 file name for each baseline
%       - name_bs = cell array containing baselines names
%       - path_comparison_results: path to output folder which will contain
%                                   plots
%       - index_inside_mask: cell array containing index of recordings inside mask for each baseline
% 
% OUTPUT:
%       
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% Code based on /retracker/validation/ by Eduard Markhoul
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:

                                            
%% --------------- HANDLING OPTIONAL INPUTS ---------------------------
if(nargin<13)
    error('Wrong number of input parameters');
end

N_baselines = size(SWH,1);

p = inputParser;
p.addParamValue('win_size_detrending',20);
p.addParamValue('flag_outliers_removal',1);
p.addParamValue('type_outliers_removal','tukey_fence',@(x)ischar(x));
p.addParamValue('IQR_times', 1.5); %number of IQR 
%--------using hampel filter
p.addParamValue('smooth_param',1);
p.addParamValue('annotation_box_active',1);
p.addParamValue('generate_plot_SSH',1);
p.addParamValue('generate_plot_SWH',1);
p.addParamValue('generate_plot_sigma0',1);
p.addParamValue('generate_plot_COR',1);
% p.addParamValue('generate_plot_Misfit',1);
default_processor_ID=cell(1,N_baselines);
default_processor_ID(:)={'MAT'};
p.addParamValue('ref_SSH',12); %reference SSH to compute error for S6 sim
p.addParamValue('ref_SWH',2); %reference SWH to compute error for S6 sim
p.addParamValue('ref_sigma0',12); %reference sigma0 to compute error for S6
default_L1B_processor=cell(1,N_baselines);
default_L1B_processor(:)={'ISD'};
p.addParamValue('L1B_processor',default_L1B_processor); %indicate L1B is being ingested from GPOD (for each baseline)
% p.addParamValue('LineStyle', 'none', @(x)ischar(x));


p.parse(varargin{:});
win_size_detrending=p.Results.win_size_detrending;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=p.Results.type_outliers_removal;
IQR_times=p.Results.IQR_times;
smooth_param=p.Results.smooth_param;
annotation_box_active=p.Results.annotation_box_active;
generate_plot_SSH=p.Results.generate_plot_SSH;
generate_plot_SWH=p.Results.generate_plot_SWH;
generate_plot_sigma0=p.Results.generate_plot_sigma0;
generate_plot_COR=p.Results.generate_plot_COR;
% generate_plot_Misfit=p.Results.generate_plot_Misfit;
ref_SSH=p.Results.ref_SSH;
ref_SWH=p.Results.ref_SWH;
ref_sigma0=p.Results.ref_sigma0;
L1B_processor = p.Results.L1B_processor;
% LineStyle = p.Results.LineStyle;
clear p;


%----------------- Define linstyles for bulk comparison -------------------
switch lower(cnf_tool.figure_format)
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

color_bs = cnf_tool.color_bs;
marker_bs = cnf_tool.marker_bs;

%% 

N_baselines = length(SWH);
                            
%% --- define some output variables as empty if not available the information -----------

if ~smooth_param
    SSH_smooth=[];
    SWH_smooth=[];
    sigma0_smooth=[];
    COR_smooth=[];  
%     Misfit_smooth = [];
    SSH_RMSE=[];
    SWH_RMSE=[];
    sigma0_RMSE=[];
    COR_RMSE=[];
%     Misfit_RMSE=[];
end

indx_LR = 0;


%% ----------------------------- Run the comparison -----------------------                                      
for i_baseline=1:N_baselines   
    
    if ~isempty(find(ismember(strsplit(filename_L2_ISR{i_baseline}, '_'), 'LR')))
        indx_LR = i_baseline;
    end
    
    % Filtering estimated geophysical parameters outside the
    % cnf_tool.performance_latitude_range specified (default: [-inf, inf])
    switch cnf_tool.L2proc{i_baseline}
        case 'ISD'
            lat_surf{i_baseline} = double(ncread(filename_L2_ISR{i_baseline},'lat_20_ku')).';
            lon_surf{i_baseline} = double(ncread(filename_L2_ISR{i_baseline},'lon_20_ku')).';
        case 'GPP'
            lat_surf{i_baseline} = double(ncread(filename_L2_ISR{i_baseline},'data_20/ku/latitude')).';
    end
    
    lat_mask{i_baseline} = find(lat_surf{i_baseline}>cnf_tool.performance_latitude_range(1) & lat_surf{i_baseline}<cnf_tool.performance_latitude_range(2));    

    lat_surf{i_baseline} = lat_surf{i_baseline}(lat_mask{i_baseline});
    lon_surf{i_baseline} = lon_surf{i_baseline}(lat_mask{i_baseline});
    SSH{i_baseline}=SSH{i_baseline}(lat_mask{i_baseline});
    SWH{i_baseline}=SWH{i_baseline}(lat_mask{i_baseline});

%     SWH{i_baseline}(SWH{i_baseline}<4)=NaN;   
    
    sigma0{i_baseline}=sigma0{i_baseline}(lat_mask{i_baseline});
    if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
        COR{i_baseline}=COR{i_baseline}(lat_mask{i_baseline});
    end
    ISD_num_surfaces = length(SWH{i_baseline});

%     fprintf('Performance analysis of %s: range latitude from %.4g to %.4g\n', name_bs{i_baseline}, lat_surf{i_baseline}(1), lat_surf{i_baseline}(end));

    %% ------------------- Outliers filtering -----------------------------
    %----------------------------------------------------------------------
    if flag_outliers_removal                                
        % compute the errors w.r.t fitting on the data using a smooth
        % function
        %------------------- SSH ----------------------------------
        [SSH{i_baseline},idx_outliers_SSH{i_baseline}]=remove_outliers(SSH{i_baseline},'type_outliers_removal', type_outliers_removal, 'IQR_times', IQR_times);
        idx_nooutliers_SSH=find(~idx_outliers_SSH{i_baseline});
        %----------------- sigma0 ---------------------------------
        [sigma0{i_baseline},idx_outliers_sigma0{i_baseline}]=remove_outliers(sigma0{i_baseline},'type_outliers_removal', type_outliers_removal, 'IQR_times',IQR_times);
        idx_nooutliers_sigma0=find(~idx_outliers_sigma0{i_baseline});
        %----------------- SWH ------------------------------------
        [SWH{i_baseline},idx_outliers_SWH{i_baseline}]=remove_outliers(SWH{i_baseline},'type_outliers_removal', type_outliers_removal, 'IQR_times', IQR_times);
        idx_nooutliers_SWH{i_baseline}=find(~idx_outliers_SWH{i_baseline});
        if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
            %----------------- COR ------------------------------------
            [COR{i_baseline},idx_outliers_COR{i_baseline}]=remove_outliers(COR{i_baseline},'type_outliers_removal', type_outliers_removal, 'IQR_times', IQR_times);
            idx_nooutliers_COR=find(~idx_outliers_COR{i_baseline});
        end

    else
        %------------------ SSH ---------------------------
        idx_outliers_SSH{i_baseline}=zeros(1,ISD_num_surfaces);
        idx_nooutliers_SSH=find(~idx_outliers_SSH{i_baseline});
        %----------------- sigma0 ---------------------------------
        idx_outliers_sigma0{i_baseline}= zeros(1,ISD_num_surfaces);
        idx_nooutliers_sigma0=find(~idx_outliers_sigma0{i_baseline});
        %----------------- SWH ------------------------------------
        idx_outliers_SWH{i_baseline}=zeros(1,ISD_num_surfaces);
        idx_nooutliers_SWH{i_baseline}=find(~idx_outliers_SWH{i_baseline});
        if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
            %----------------- COR ------------------------------------
            idx_outliers_COR{i_baseline}=zeros(1,ISD_num_surfaces);
            idx_nooutliers_COR=find(~idx_outliers_COR{i_baseline});
        end        
        
    end
     
    
    %% --------------- Smoothing geophysical retrievals ---------------
    if smooth_param
        %------------------ SSH ---------------------------
        SSH_smooth{i_baseline} = smooth(SSH{i_baseline},win_size_detrending).';
        %----------------- sigma0 ---------------------------------
        sigma0_smooth{i_baseline} = smooth(sigma0{i_baseline},win_size_detrending).';
        %----------------- SWH ------------------------------------
        SWH_smooth{i_baseline} = smooth(SWH{i_baseline},win_size_detrending).';
        if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
            %----------------- COR ------------------------------------
            COR_smooth{i_baseline} = smooth(COR{i_baseline},win_size_detrending).';
        end
    end
    
    
    %% ------------------- Compute the std BLOCK-WISE ---------------------
    %----------------------------------------------------------------------
    num_boxes=floor(ISD_num_surfaces/win_size_detrending);
    for i_box=1:(num_boxes+1)
        
        init_sample=max([(i_box-1)*win_size_detrending+1,1]);
        last_sample=min([(i_box-1)*win_size_detrending+win_size_detrending,ISD_num_surfaces]);
        
        SSH_std{i_baseline}(i_box)=nanstd(detrend(SSH{i_baseline}(idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
        SWH_std{i_baseline}(i_box)=nanstd(detrend(SWH{i_baseline}(idx_nooutliers_SWH{i_baseline}(idx_nooutliers_SWH{i_baseline}>=init_sample & idx_nooutliers_SWH{i_baseline}<=last_sample))));
        sigma0_std{i_baseline}(i_box)=nanstd(detrend(sigma0{i_baseline}(idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
        
        SSH_mean{i_baseline}(i_box)=nanmean((SSH{i_baseline}(idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
        SWH_mean{i_baseline}(i_box)=nanmean((SWH{i_baseline}(idx_nooutliers_SWH{i_baseline}(idx_nooutliers_SWH{i_baseline}>=init_sample & idx_nooutliers_SWH{i_baseline}<=last_sample))));
        sigma0_mean{i_baseline}(i_box)=nanmean((sigma0{i_baseline}(idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));

        if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
            COR_std{i_baseline}(i_box)=nanstd(detrend(COR{i_baseline}(idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
            COR_mean{i_baseline}(i_box)=nanmean((COR{i_baseline}(idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
        end        
    end

    
    %% --------- Compute the RMSE FITTING ---------------------------------
    if smooth_param     
        SSH_RMSE(i_baseline)=sqrt(nanmean((SSH{i_baseline}-SSH_smooth{i_baseline}).^2));
        SWH_RMSE(i_baseline)=sqrt(nanmean((SWH{i_baseline}-SWH_smooth{i_baseline}).^2));
        sigma0_RMSE(i_baseline)=sqrt(nanmean((sigma0{i_baseline}-sigma0_smooth{i_baseline}).^2));
        
        if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
            COR_RMSE(i_baseline)=sqrt(nanmean((COR{i_baseline}-COR_smooth{i_baseline}).^2));
        end
    end
    
    
    %% --------- Compute the mean std over track equivalent to RMSE--------

    SSH_std_mean(i_baseline)=nanmean(SSH_std{i_baseline});
    SWH_std_mean(i_baseline)=nanmean(SWH_std{i_baseline});
    sigma0_std_mean(i_baseline)=nanmean(sigma0_std{i_baseline});
    
    if strcmp(cnf_tool.L2proc{i_baseline}, 'ISD')
        COR_std_mean(i_baseline)=nanmean(COR_std{i_baseline});
    end
    
    clear idx_nooutliers_SSH idx_nooutliers_sigma0 idx_nooutliers_COR; %idx_nooutliers_Misfit;
    clear idx_nooutliers_SSH_ESA idx_nooutliers_sigma0_ESA;
    clear idx_nooutliers_SSH_ISR_ESA idx_nooutliers_sigma0_ISR_ESA;
    clear idx_nooutliers_SSH_STL idx_nooutliers_sigma0_STL idx_nooutliers_SWH_STL;
    clear idx_nooutliers_SSH_GPOD idx_nooutliers_sigma0_GPOD idx_nooutliers_SWH_GPOD;
    clear idx_nooutliers_amp_fit idx_nooutliers_max_wvfm idx_nooutliers_Pu;
end %loop over the channels

clear ISD_i2q2_meas;
ext_baselines_comp=strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_');


%% ---------------- PLOTING RESULTS -----------------------------------
%----------------------------------------------------------------------
%--------------------- SSH --------------------------------------------

if N_baselines > 1
    index_baseline_compare_aux = 1:N_baselines;
    for i_baseline=1:N_baselines
        if i_baseline == indx_LR
            index_baseline_compare_aux(i_baseline) = [];
        end
    end
    index_baseline_compare = nchoosek(index_baseline_compare_aux,2);
end

[~,file_id,~]=fileparts(filename_L2_ISR{1});
aux=strsplit(file_id, '_');
text_title = {aux{1}, 'L1', aux{2}};

% finfo = ncinfo(filename_L2_ISR{1});
% if any(ismember({finfo.Attributes.Name}, {'cycle_number'})) 
%     cycle_number = ncreadatt(filename_L2_ISR{1}, '/','cycle_number');
%     pass_number = ncreadatt(filename_L2_ISR{1}, '/','pass_number');
%     text_title = [text_title, sprintf(' - cycle %d pass %d', cycle_number, pass_number)]; 
% end


text_interpreter=get(0, 'defaultAxesTickLabelInterpreter'); %cnf_p.text_interpreter;
plot_baseline = 1;

% write statistics in .txt file
fid_stat_param=fopen([path_comparison_results,file_id,'_stat_param.txt'],'w');
fprintf(fid_stat_param,'------------------------------------------------------------------------------------\n');
fprintf(fid_stat_param,'----------- Computed statistics of retrieved geophysical parameters ----------------\n');
fprintf(fid_stat_param,'------------------------------------------------------------------------------------\n');
fprintf(fid_stat_param,'%35s%s%7s%s%7s%s\n','', 'RMSE', '', '<std>', '', '<mean>');

if generate_plot_SSH 
    
    %--------------------- SSH --------------------------------------------
    f1=figure;

    legend_text={''};
    text_in_textbox={''};

    plot_baseline = 1;
    for b=1:N_baselines
        bias_compensation=0.0;%nanmean(SSH(b,:))-12.0;
          
        coef_width=1*0.6^(plot_baseline-1);
        plot_baseline = plot_baseline + 1;
        plt=plot(lat_surf{b},SSH{b}-bias_compensation,'Marker',char(marker_bs(b)),'Color',color_bs(b,:),...
            'LineStyle',cnf_tool.LineStyle{b}, 'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
        plt.Color(4) = 1; % transparency
        hold on;
        legend_text=[legend_text,name_bs(b)];
%         text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [m]\nstd = %.4g [m]\nBias = %.4g [m]', ...
%             SSH_RMSE(b), SSH_std_mean(b), nanmean(SSH_mean{b}-bias_compensation)-ref_SSH)];
        text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [m]\nstd = %.4g [m]\nmean = %.4g [m]', ...
            SSH_RMSE(b), SSH_std_mean(b), nanmean(SSH_mean{b}))];
        if b ~= N_baselines
           text_in_textbox = [text_in_textbox, sprintf('\n')]; 
        end
        
        fprintf(fid_stat_param,'%s [%s] %15s %.5f %3s %.5f %3s %.5f\n','SSH', char(name_bs(b)), '', SSH_RMSE(b), '', SSH_std_mean(b), '', nanmean(SSH_mean{b}));
    end
    plot_baseline = 1;
    title(strjoin(text_title), 'Interpreter',text_interpreter);
    leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside');
    pos_leg=get(leg,'Position');
    xlabel('Latitude [deg.]','Interpreter',text_interpreter); ylabel(strcat('SSH',' [m]'),'Interpreter',text_interpreter);
    text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
    xlim([min(lat_surf{1}), max(lat_surf{1})]);
    if annotation_box_active
        h=annotation('textbox',[0.999*pos_leg(1),pos_leg(2)-1.4*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
            'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
        h.LineWidth = 0.5;
    end
    ax=gca;
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual'; 

    addlogoisardSAT('plot');
%     addlogoESA('plot');
    
    print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_SSH', file_ext]); % figure_format
    if cnf_tool.save_figure_format_fig
       savefig([path_comparison_results,file_id,'_SSH', '.fig']) 
    end
        
    %--------------------- SSH difference ----------------------------------------

    if N_baselines > 1
        f2=figure;

        legend_text={''};
        text_in_textbox={''};

        for b=1:size(index_baseline_compare,1) 
            b1=index_baseline_compare(b,1);
            b2=index_baseline_compare(b,2);

            bias_compensation=0.0;%nanmean(SSH(b,:))-12.0;

            coef_width=1*0.6^(plot_baseline-1);
            plot_baseline = plot_baseline + 1;
            
            plt=plot(lat_surf{b1},SSH{b1}-SSH{b2},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),...
                'LineStyle',cnf_tool.LineStyle{b1}, 'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
            hold on;
            legend_text=[legend_text,strcat(sprintf('%s vs %s', char(name_bs(b1)), char(name_bs(b2))))];
            text_in_textbox=[text_in_textbox, sprintf('%s vs %s:', char(name_bs(b1)), char(name_bs(b2))), sprintf('mean = %.4g [m]\nstd = %.4g [m]', ...
                nanmean(SSH_mean{b1}-SSH_mean{b2}), nanstd(SSH_mean{b1}-SSH_mean{b2}))];
            if b ~= size(index_baseline_compare,1) 
               text_in_textbox = [text_in_textbox, sprintf('\n')]; 
            end
            
            fprintf(fid_stat_param,'%s [%s-%s] %5s %s %2s %.5f %2s %.5f\n','SSH', char(name_bs(b1)), char(name_bs(b2)), '', 'NaN', '', nanstd(SSH_mean{b1}-SSH_mean{b2}), '', nanmean(SSH_mean{b1}-SSH_mean{b2}));

        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]','Interpreter',text_interpreter); ylabel(strcat('SSH difference',' [m]'),'Interpreter',text_interpreter);
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        xlim([min(lat_surf{1}), max(lat_surf{1})]);
        if annotation_box_active
            h=annotation('textbox',[pos_leg(1),pos_leg(2)-3*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
                'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
            h.LineWidth = 0.5;
        end
        ax=gca;
        ax.XTickMode = 'manual';
        ax.YTickMode = 'manual';
        ax.ZTickMode = 'manual';
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual'; 

        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differenceSSH', file_ext]); % figure_format
        if cnf_tool.save_figure_format_fig
           savefig([path_comparison_results,file_id,'_differenceSSH', '.fig']) 
        end
        close(f2)
    end
    close(f1)
    
%     % % map parameters value on a geographical map - see Pablo's plots
%     f1=figure;
%     axesm eckert4;
%     framem; gridm;
%     axis off
%     scatterm(lat_dat, lon_dat, 30, res_P4_yaw, 'filled');
%     geoshow('landareas.shp', 'FaceColor', 'none', 'EdgeColor', 'black');
%     hcb = colorbar('southoutside');
%     title(['My title bla bla bla']);
%     set(get(hcb,'Xlabel'),'String','Yaw residual [deg]')
%     addlogoisardSAT('map_attitude') (edited) 
    
end

if generate_plot_SWH
%--------------------- SWH --------------------------------------------
    f1=figure;
    legend_text={''};
    text_in_textbox={''};

    for b=1:N_baselines
        
        coef_width=1*0.6^(plot_baseline-1);
        plot_baseline = plot_baseline + 1;
        plt=plot(lat_surf{b},SWH{b},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),'LineStyle',cnf_tool.LineStyle{b},...
            'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
        plt.Color(4) = 1; % transparency
        hold on;
        legend_text=[legend_text,name_bs(b)];
%         text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [m]\nstd = %.4g [m]\nBias = %.4g [m]',...
%             SWH_RMSE(b), SWH_std_mean(b), nanmean(SWH_mean{b})-ref_SWH)];
        text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [m]\nstd = %.4g [m]\nmean = %.4g [m]',...
            SWH_RMSE(b), SWH_std_mean(b), nanmean(SWH_mean{b}))];
        if b ~= N_baselines
           text_in_textbox = [text_in_textbox, sprintf('\n')]; 
        end
        
        fprintf(fid_stat_param,'%s [%s] %15s %.5f %3s %.5f %3s %.5f\n','SWH', char(name_bs(b)), '', SWH_RMSE(b), '', SWH_std_mean(b), '', nanmean(SWH_mean{b}));

    end
    plot_baseline = 1;
    title(strjoin(text_title), 'Interpreter',text_interpreter);
    leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
    pos_leg=get(leg,'Position');

    xlabel('Latitude [deg.]','Interpreter',text_interpreter); ylabel('SWH [m]','Interpreter',text_interpreter);
    xlim([min(lat_surf{1}), max(lat_surf{1})]);
    text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
    if annotation_box_active
        h=annotation('textbox',[0.999*pos_leg(1),pos_leg(2)-1.4*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
            'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
        h.LineWidth = 0.5;
    end
    ax=gca;
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    
    addlogoisardSAT('plot');

    print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_SWH', file_ext]); % figure_format
    if cnf_tool.save_figure_format_fig
       savefig([path_comparison_results,file_id,'_SWH', '.fig']) 
    end
    
    %--------------- SWH difference --------------------------------------------
    if N_baselines > 1
        f2=figure;
        legend_text={''};
        text_in_textbox={''};

        for b=1:size(index_baseline_compare,1) 
            b1=index_baseline_compare(b,1);
            b2=index_baseline_compare(b,2);

            coef_width=1*0.6^(plot_baseline-1);
            plot_baseline = plot_baseline + 1;
            plt=plot(lat_surf{b1},SWH{b1}-SWH{b2},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),'LineStyle',cnf_tool.LineStyle{b1}, ...
                'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
                    plt.Color(4) = 1; % transparency
            hold on;
            legend_text=[legend_text,sprintf('%s vs %s', char(name_bs(b1)), char(name_bs(b2)))];
            text_in_textbox=[text_in_textbox, sprintf('%s vs %s:', char(name_bs(b1)), char(name_bs(b2))), sprintf('mean = %.4g [m]\nstd = %.4g [m]', ...
                nanmean(SWH_mean{b1}-SWH_mean{b2}), nanstd(SWH_mean{b1}-SWH_mean{b2}))];
            if b ~= size(index_baseline_compare,1) 
               text_in_textbox = [text_in_textbox, sprintf('\n')]; 
            end
                        
            fprintf(fid_stat_param,'%s [%s-%s] %5s %s %2s %.5f %2s %.5f\n','SWH', char(name_bs(b1)), char(name_bs(b2)), '', 'NaN', '', ...
                nanstd(SWH_mean{b1}-SWH_mean{b2}), '', nanmean(SWH_mean{b1}-SWH_mean{b2}));

        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
        xlim([min(lat_surf{1}), max(lat_surf{1})]);
        xlabel('Latitude [deg.]','Interpreter',text_interpreter); ylabel('SWH difference [m]','Interpreter',text_interpreter);
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        if annotation_box_active
            h=annotation('textbox',[pos_leg(1),pos_leg(2)-3*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
                'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
            h.LineWidth = 0.5;
        end
        ax=gca;
        ax.XTickMode = 'manual';
        ax.YTickMode = 'manual';
        ax.ZTickMode = 'manual';
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual';

        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differenceSWH', file_ext]); % figure_format
        if cnf_tool.save_figure_format_fig
           savefig([path_comparison_results,file_id,'_differenceSWH', '.fig']) 
        end
        close(f2);
        
        % -------------- SWH difference vs SWH --------------------
        f2=figure;
        legend_text={''};
        for b=1:size(index_baseline_compare,1) 
            b1=index_baseline_compare(b,1);
            b2=index_baseline_compare(b,2);

            coef_width=1*0.6^(plot_baseline-1);
            plot_baseline = plot_baseline + 1;
            plt=plot(SWH{b1},abs(SWH{b1}-SWH{b2}),'Marker',char(marker_bs(b)),'Color',color_bs(b,:),'LineStyle','none', ...
                'MarkerSize', cnf_tool.default_markersize);
            plt.Color(4) = 1; % transparency
            hold on;
            legend_text=[legend_text,sprintf('%s vs %s', char(name_bs(b1)), char(name_bs(b2)))];
        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
        xlim([0 11]); ylim([0,nanmean(abs(SWH{b1}-SWH{b2}))+2*nanstd(abs(SWH{b1}-SWH{b2}))]);
        xlabel('SWH [m]','Interpreter',text_interpreter); ylabel('|SWH difference| [m]','Interpreter',text_interpreter);
        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differenceSWH_vs_SWH', file_ext]); % figure_format
        close(f2);
        
        % -------------- SWH difference histogram for each baseline  --------------------
        SWH_original = SWH;  
        f2=figure;
        legend_text={''};
        text_in_textbox={''};   
        gleq = '<';
        for ii=1:2  % set SWH>4 and SWH<4 to NaN and recompute mean and std
            clear SWH_mean SWH
            for b=1:N_baselines        
                SWH{b}= SWH_original{b};
                
%                 [SWH{b},idx_outliers_SWH{b}]=remove_outliers(SWH{b},'type_outliers_removal', type_outliers_removal, 'IQR_times', IQR_times);
%                 idx_nooutliers_SWH{b}=find(~idx_outliers_SWH{b});
                
                if ii==1 
                    SWH{b}(SWH{b}>=4)=NaN; % Hs < 4
                else 
                    SWH{b}(SWH{b}<4)=NaN; % Hs>=4
                    gleq = '>=';
                end
                
                num_boxes=floor(ISD_num_surfaces/win_size_detrending);
                for i_box=1:(num_boxes+1)
                    init_sample=max([(i_box-1)*win_size_detrending+1,1]);
                    last_sample=min([(i_box-1)*win_size_detrending+win_size_detrending,ISD_num_surfaces]);
                    SWH_mean{b}(i_box)=nanmean((SWH{b}(idx_nooutliers_SWH{b}(idx_nooutliers_SWH{b}>=init_sample & idx_nooutliers_SWH{b}<=last_sample))));
                end      
            end
            
            % plot
            for b=1:size(index_baseline_compare,1)
                b1=index_baseline_compare(b,1);
                b2=index_baseline_compare(b,2);

                plt=histogram(SWH{b1}-SWH{b2},'BinWidth', 0.01, 'Normalization','probability', 'FaceColor', color_bs(ii,:), 'EdgeColor', 'none');
%                 plt.Color(4) = 1; % transparency
                hold on;
                legend_text=[legend_text,sprintf('%s vs %s, H_s %s 4', char(name_bs(b1)), char(name_bs(b2)), gleq)];
                text_in_textbox=[text_in_textbox, sprintf('%s vs %s:', char(name_bs(b1)), char(name_bs(b2))), sprintf('mean = %.4g [m]\nstd = %.4g [m]', ...
                    nanmean(SWH_mean{b1}-SWH_mean{b2}), nanstd(SWH_mean{b1}-SWH_mean{b2}))]; 
                if b ~= size(index_baseline_compare,1) 
                   text_in_textbox = [text_in_textbox, sprintf('\n\n')]; 
                end

            end
        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeast','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        if annotation_box_active
            h=annotation('textbox',[pos_leg(1),pos_leg(2)-3*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
                'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
            h.LineWidth = 0.5;
        end
        xlim([-1 1]);
%         xlim([-(nanmean(abs(SWH{b1}-SWH{b2}))+2*nanstd(abs(SWH{b1}-SWH{b2}))),nanmean(abs(SWH{b1}-SWH{b2}))+2*nanstd(abs(SWH{b1}-SWH{b2}))]);
        xlabel('SWH difference [m]','Interpreter',text_interpreter); ylabel('frequency','Interpreter',text_interpreter);
        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differenceSWH_histogram', file_ext]); % figure_format
        close(f2);
        
    end
    close(f1);
    
end

if generate_plot_sigma0
    %--------------------- sigma0 --------------------------------------------
    f1=figure;
    legend_text={''};
    text_in_textbox={''};

    for b=1:N_baselines
        bias_compensation=0.0;%mean(sigma0(b,:))-12.0;
        coef_width=1*0.6^(plot_baseline-1);
        plot_baseline = plot_baseline + 1;
        plt=plot(lat_surf{b},sigma0{b}-bias_compensation,'Marker',char(marker_bs(b)),'Color',color_bs(b,:),...
            'LineStyle',cnf_tool.LineStyle{b}, 'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
        plt.Color(4) = 1; % transparency
        hold on;
        legend_text=[legend_text,name_bs(b)];
        text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [dB]\nstd = %.4g [dB]\nMean = %.4g [dB]',...
            sigma0_RMSE(b), sigma0_std_mean(b), nanmean(sigma0_mean{b}))];
%         text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), sprintf('RMSE = %.4g [dB]\nstd = %.4g [dB]\nmean = %.4g [dB]',...
%             sigma0_RMSE(b), sigma0_std_mean(b), nanmean(sigma0_mean{b}))];
        if b ~= N_baselines
           text_in_textbox = [text_in_textbox, sprintf('\n')]; 
        end

        fprintf(fid_stat_param,'%s [%s] %15s %.5f %3s %.5f %3s %.5f\n','sigma0', char(name_bs(b)), '', sigma0_RMSE(b), '', sigma0_std_mean(b), '', ...
                                                                                                    nanmean(sigma0_mean{b}));

    end
    plot_baseline = 1;
    title(strjoin(text_title), 'Interpreter',text_interpreter);
    leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
    pos_leg=get(leg,'Position');
    xlim([min(lat_surf{1}), max(lat_surf{1})]);
    xlabel('Latitude [deg.]','Interpreter',text_interpreter); 
    if strcmp(text_interpreter, 'latex')
        ylabel('$\sigma^0$ [dB]','Interpreter',text_interpreter);    
    else
        ylabel('\sigma^0 [dB]','Interpreter',text_interpreter);
    end
    text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));

    if annotation_box_active
       h=annotation('textbox',[0.999*pos_leg(1),pos_leg(2)-1.4*pos_leg(4),pos_leg(3),pos_leg(4)],'String',...
           text_in_textbox,'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
        h.LineWidth = 0.5;
    end
    ax=gca;
    ax.XTickMode = 'manual';
    ax.YTickMode = 'manual';
    ax.ZTickMode = 'manual';
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    
    addlogoisardSAT('plot');
    
    print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_sigma0', file_ext]); % figure_format
    if cnf_tool.save_figure_format_fig
        savefig([path_comparison_results,file_id,'_sigma0', '.fig']) 
    end
    
    %--------------------- sigma0 difference ------------------------------------
    if N_baselines > 1
        f2=figure;
        legend_text={''};
        text_in_textbox={''};

        for b=1:size(index_baseline_compare,1) 
            b1=index_baseline_compare(b,1);
            b2=index_baseline_compare(b,2);

            bias_compensation=0.0;%nanmean(SSH(b,:))-12.0;
            coef_width=1*0.6^(plot_baseline-1);
            plot_baseline = plot_baseline + 1;
            plt=plot(lat_surf{b1},sigma0{b1}-sigma0{b2},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),...
                'LineStyle',cnf_tool.LineStyle{b1}, 'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
                    plt.Color(4) = 1; % transparency
            hold on;
            legend_text=[legend_text,sprintf('%s vs %s', char(name_bs(b1)), char(name_bs(b2)))];
            text_in_textbox=[text_in_textbox, sprintf('%s vs %s:', char(name_bs(b1)), char(name_bs(b2))), sprintf('Mean = %.4g [dB]\nstd = %.4g [dB]', ...
                nanmean(sigma0_mean{b1}-sigma0_mean{b2}), nanstd(sigma0_mean{b1}-sigma0_mean{b2}))];
            if b ~= size(index_baseline_compare,1) 
               text_in_textbox = [text_in_textbox, sprintf('\n')]; 
            end
            
            fprintf(fid_stat_param,'%s [%s-%s] %5s %s %2s %.5f %2s %.5f\n','sigma0', char(name_bs(b1)), char(name_bs(b2)), '', 'NaN', '', ...
                nanstd(sigma0_mean{b1}-sigma0_mean{b2}), '', nanmean(sigma0_mean{b1}-sigma0_mean{b2}));
        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
        xlim([min(lat_surf{1}), max(lat_surf{1})]);
        xlabel('Latitude [deg.]','Interpreter',text_interpreter); 
        if strcmp(text_interpreter, 'latex')
            ylabel('$\sigma^0$ difference [dB]','Interpreter',text_interpreter);
        else
             ylabel('\sigma^0 difference [dB]','Interpreter',text_interpreter);
        end
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));

        if annotation_box_active
           h=annotation('textbox',[pos_leg(1),pos_leg(2)-3*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
               'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
            h.LineWidth = 0.5;
        end
        ax=gca;
        ax.XTickMode = 'manual';
        ax.YTickMode = 'manual';
        ax.ZTickMode = 'manual';
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual';

        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differencesigma0', file_ext]); % figure_format
        if cnf_tool.save_figure_format_fig
            savefig([path_comparison_results,file_id,'_differencesigma0', '.fig']) 
        end    
        close(f2)
    end
    close(f1)
end

if generate_plot_COR && ~any(ismember(cnf_tool.L2proc,'GPP'))
%--------------------- COR --------------------------------------------
    f1=figure;
    legend_text={''};
    text_in_textbox={''};

    if strcmp(text_interpreter, 'latex')
        percent = ' [\%]';
    else
        percent = ' [%]';
    end
    for b=1:N_baselines
        
        coef_width=1*0.6^(plot_baseline-1);
        plot_baseline = plot_baseline + 1;
        plt=plot(lat_surf{b},COR{b},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),'LineStyle',cnf_tool.LineStyle{b}, ...
            'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
                plt.Color(4) = 1; % transparency
        hold on;
        legend_text=[legend_text,name_bs(b)];
        text_in_textbox=[text_in_textbox, strcat(char(name_bs(b)), ':'), strcat(sprintf('RMSE = %.4g', COR_RMSE(b)), percent), ...
            strcat(sprintf('std = %.4g', COR_std_mean(b)), percent), strcat(sprintf('mean = %.4g',nanmean(COR_mean{b})), percent)];
        if b ~= N_baselines
           text_in_textbox = [text_in_textbox, sprintf('\n')]; 
        end
        
        fprintf(fid_stat_param,'%s [%s] %15s %.5f %3s %.5f %3s %.5f\n','COR', char(name_bs(b)), '', COR_RMSE(b), '', COR_std_mean(b), '', ...
                                                                                                    nanmean(COR_mean{b}));
    end
    plot_baseline = 1;
    title(strjoin(text_title), 'Interpreter',text_interpreter);
    leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
    pos_leg=get(leg,'Position');
    ylim([0.995*min([COR{:}]),min(100, 1.005*max([COR{:}]))]);
    xlim([min(lat_surf{1}), max(lat_surf{1})]);
    xlabel('Latitude [deg.]','Interpreter',text_interpreter); 
    if strcmp(text_interpreter, 'latex')
        ylabel('$\rho$ [\%]','Interpreter',text_interpreter);
    else
        ylabel('\rho [%]','Interpreter',text_interpreter);
    end
    text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
    if annotation_box_active
        h=annotation('textbox',[0.999*pos_leg(1),pos_leg(2)-1.4*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
            'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
        h.LineWidth = 0.5;
    end
    
    addlogoisardSAT('plot');
    
    print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_COR', file_ext]); % figure_format
    if cnf_tool.save_figure_format_fig
        savefig([path_comparison_results,file_id,'_COR', '.fig']) 
    end
    
    %--------------------- COR difference --------------------------------
    if N_baselines > 1
        f2=figure;
        legend_text={''};
        text_in_textbox={''};

        for b=1:size(index_baseline_compare,1) 
            b1=index_baseline_compare(b,1);
            b2=index_baseline_compare(b,2);

            bias_compensation=0.0;%nanmean(SSH(b,:))-12.0;
            coef_width=1*0.6^(plot_baseline-1);
            plot_baseline = plot_baseline + 1;
            plt=plot(lat_surf{b1},COR{b1}-COR{b2},'Marker',char(marker_bs(b)),'Color',color_bs(b,:),'LineStyle',cnf_tool.LineStyle{b1}, ...
                'MarkerSize', cnf_tool.default_markersize, 'LineWidth', coef_width*cnf_tool.default_linewidth);
                    plt.Color(4) = 1; % transparency
            hold on;
            legend_text=[legend_text,sprintf('%s vs %s', char(name_bs(b1)), char(name_bs(b2)))];
            text_in_textbox=[text_in_textbox, sprintf('%s vs %s:', char(name_bs(b1)), char(name_bs(b2))), strcat(sprintf('Mean = %.4g', ...
                nanmean(COR_mean{b1}-COR_mean{b2})), percent, sprintf('\nstd = %.4g', ...
                nanstd(COR_mean{b1}-COR_mean{b2})), percent)];
            if b ~= size(index_baseline_compare,1) 
               text_in_textbox = [text_in_textbox, sprintf('\n')]; 
            end
            
            fprintf(fid_stat_param,'%s [%s-%s] %5s %s %2s %.5f %2s %.5f\n','COR', char(name_bs(b1)), char(name_bs(b2)), '', 'NaN', '', ...
                nanstd(COR_mean{b1}-COR_mean{b2}), '', nanmean(COR_mean{b1}-COR_mean{b2}));
        end
        plot_baseline = 1;
        title(strjoin(text_title), 'Interpreter',text_interpreter);
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',cnf_tool.legend_fontsize);
        pos_leg=get(leg,'Position');
    %     ylim([99,100]);
        xlim([min(lat_surf{1}), max(lat_surf{1})]);
        xlabel('Latitude [deg.]','Interpreter',text_interpreter); 
        if strcmp(text_interpreter, 'latex')
            ylabel('$\rho$ difference [\%]','Interpreter',text_interpreter);
        else
            ylabel('\rho difference [%]','Interpreter',text_interpreter);
        end
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        if annotation_box_active
            h=annotation('textbox',[pos_leg(1),pos_leg(2)-3*pos_leg(4),pos_leg(3),pos_leg(4)],'String',text_in_textbox,...
                'FitBoxToText','on','FontSize',cnf_tool.textbox_fontsize, 'Interpreter',text_interpreter);
            h.LineWidth = 0.5;
        end

        addlogoisardSAT('plot');

        print(print_file,cnf_tool.res_fig,[path_comparison_results,file_id,'_differenceCOR', file_ext]); % figure_format
        if cnf_tool.save_figure_format_fig
            savefig([path_comparison_results,file_id,'_COR', '.fig']) 
        end

        close(f2);
    end
    close(f1);
end

fclose(fid_stat_param);

% save(strcat(path_comparison_results,file_id,'_Perf_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_')));

end