function [lat_surf, lon_surf, SSH, SWH, sigma0, COR, SSH_RMSE, SWH_RMSE, sigma0_RMSE, COR_RMSE, SSH_std_mean, SWH_std_mean, sigma0_std_mean, COR_std_mean, ...
                SSH_mean, SWH_mean, sigma0_mean, COR_mean]=performance_baselines_S6_dynamic(SSH, SWH, sigma0, COR, filename_L1_ISR, index_inside_mask,...
                                            name_bs,cnf_tool,varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for the comparison of the geophsycial retrievals of
% different input baselines of isardSAT:
% exctracting the performance in terms of noise metric (as std or/and RMSE)
% model_fit_to_power_...._dynamic calls this function to plot estimates
% below fitted curve and produce a .gif
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 10/12/2020
% This software is built within the Sentinel-6 P4 L1 GPP project -
% commissioning analysis
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       - SSH, SWH, sigma0, COR = N_baselines x num. records array with retrieved parameters from L2 product
%       - filename_L1_ISR = cell array containg L1 file name for each baseline
%       - name_bs = cell array containing baselines names
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
if(nargin<8)
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
    
    if ~isempty(find(ismember(strsplit(filename_L1_ISR{i_baseline}, '_'), 'LR')))
        indx_LR = i_baseline;
    end
    
    % Filtering estimated geophysical parameters outside the
    % cnf_tool.performance_latitude_range specified (default: [-inf, inf])
    lat_surf{i_baseline} = double(ncread(filename_L1_ISR{i_baseline},'data_20/ku/latitude')).';
    lon_surf{i_baseline} = double(ncread(filename_L1_ISR{i_baseline},'data_20/ku/longitude')).';
    lat_mask{i_baseline} = find(lat_surf{i_baseline}>cnf_tool.performance_latitude_range(1) & lat_surf{i_baseline}<cnf_tool.performance_latitude_range(2));    

    lat_surf{i_baseline} = lat_surf{i_baseline}(lat_mask{i_baseline});
    lon_surf{i_baseline} = lon_surf{i_baseline}(lat_mask{i_baseline});
    SSH{i_baseline}=SSH{i_baseline}(lat_mask{i_baseline});
    SWH{i_baseline}=SWH{i_baseline}(lat_mask{i_baseline});
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
        idx_nooutliers_SWH=find(~idx_outliers_SWH{i_baseline});
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
        idx_nooutliers_SWH=find(~idx_outliers_SWH{i_baseline});
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
        SWH_std{i_baseline}(i_box)=nanstd(detrend(SWH{i_baseline}(idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
        sigma0_std{i_baseline}(i_box)=nanstd(detrend(sigma0{i_baseline}(idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
        
        SSH_mean{i_baseline}(i_box)=nanmean((SSH{i_baseline}(idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
        SWH_mean{i_baseline}(i_box)=nanmean((SWH{i_baseline}(idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
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
    
    clear idx_nooutliers_SSH idx_nooutliers_SWH idx_nooutliers_sigma0 idx_nooutliers_COR; %idx_nooutliers_Misfit;
    clear idx_nooutliers_SSH_ESA idx_nooutliers_sigma0_ESA;
    clear idx_nooutliers_SSH_ISR_ESA idx_nooutliers_sigma0_ISR_ESA;
    clear idx_nooutliers_SSH_STL idx_nooutliers_sigma0_STL idx_nooutliers_SWH_STL;
    clear idx_nooutliers_SSH_GPOD idx_nooutliers_sigma0_GPOD idx_nooutliers_SWH_GPOD;
    clear idx_nooutliers_amp_fit idx_nooutliers_max_wvfm idx_nooutliers_Pu;
end %loop over the channels

clear ISD_i2q2_meas;
ext_baselines_comp=strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_');

end