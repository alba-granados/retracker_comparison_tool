%% Comparison of GR for different scenarios
%path pointed to where the L2 data (specific structure exist subfolders HR_RAW, HR_RMC, LR_RMC) and ground
%truth of slopped scenairo for OS16 scenario
common_path ='\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Performance_S6\geophysical_retrievals\inputs\'; 
%path where to save the output figures
resultPath  ='\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Performance_S6\geophysical_retrievals\results\'; 

%common_path_bis='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1_MAT/';


mkdir(resultPath);
% A template of the excel output file shall exist in the resultPath 
%(color code of excel: red not wihtin limits performacne and green acceptable are configured in the format of 
%the template based on the numerical information being updated by the current matlab script)      
file_xls=strcat(resultPath,'performance_AR_GPP.xlsx');
row_std_SSH  = '3';
row_std_sig0 = '5';
row_std_SWH  = '7';

row_bias_SSH  = '10';
row_bias_sig0 = '12';
row_bias_SWH  = '14';
        



%% S6A_OS10
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS10_P4__HR__L2_20170305T065123_20170305T065233_0002_isd.nc',...
                     'HR_RMC/data/S6A_OS10_P4__HR__L2_20170305T065123_20170305T065233_0002_isd.nc',...
                     'LR_RMC/data/S6A_OS10_P4__LRRML2_20170305T065123_20170305T065233_0002_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS10_0001_1400RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS10_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                 
ID_GR_file_string = 'S6A_OS10';
ref_SSH    = 12; % m
ref_SWH    = 1; % m
ref_sigma0 = 11; % dB
S6A_OS10=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);   

column_xls='B';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS10.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS10.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS10.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS10.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS10.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS10.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end

%% S6A_OS16
ID_baselines={'RAW','RMC','LR-RMC'};%,'LR-RMC ad. Dopp.'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS16_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS16_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS16_P4__LRRML2_20170305T065123_20170305T065238_0001_isd.nc'});
%input_files(4)={'C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\AR\results\L2_L1_MAT\S6A_OS16_0001_1500RC_2017\LR_RMC_slope_corr\data\S6A_OS16_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'};                    
ID_GR_file_string = 'S6A_OS16';
% ref_SSH_whole    = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','ground_truth_SSH_20_ku').'; % m
% ref_SSH_lat      = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','latitude_20_ku').'; % deg
ref_SSH_whole    = ncread(strcat(common_path,'S6A_OS11_0001_slope_ground_truth_file.nc'),'ground_truth_SSH_20_ku').'; % m
ref_SSH_lat      = ncread(strcat(common_path,'S6A_OS11_0001_slope_ground_truth_file.nc'),'latitude_20_ku').'; % deg
ref_SWH    = 1; % m
ref_sigma0 = 11; % dB
S6A_OS16=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0,...
            'ref_SSH_whole',ref_SSH_whole,'lat_SSH_whole',ref_SSH_lat);   
        
column_xls='L';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS16.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS16.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS16.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS16.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS16.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS16.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end




%% S6A_OS20
%common_path ='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1_GPP/'; 
%resultPath  ='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1_GPP/GR_performance_OSTS2018/'; 
ID_baselines={'RAW','RMC','LR-RMC'};

%ID_baselines={'RAW','RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS20_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS20_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS20_P4__LRRML2_20170305T065123_20170305T065233_0001_isd.nc'});
%input_files(4)={'C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1_MAT/S6A_OS20_0001_1400RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS20_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'};                  
title_str = 'SSH = 12 m, SWH = 2 m, \sigma^0 = 11 dB';
ID_GR_file_string = 'S6A_OS20';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS20=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0,'title_str',title_str); 
        
        
column_xls='D';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS20.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS20.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS20.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS20.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS20.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS20.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end
        
%% S6A_OS21
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS21_P4__HR__L2_20170305T070423_20170305T070628_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS21_P4__HR__L2_20170305T070423_20170305T070628_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS21_P4__LRRML2_20170305T070423_20170305T070628_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS21_0001_2500RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS21_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                  
ID_GR_file_string = 'S6A_OS21';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS21=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);     

column_xls='E';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS21.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS21.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS21.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS21.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS21.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS21.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end        
        
%% S6A_OS22
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS22_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS22_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS22_P4__LRRML2_20170305T065123_20170305T065233_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS22_0001_1400RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS22_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                   
ID_GR_file_string = 'S6A_OS22';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS22=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);    
        

column_xls='F';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS22.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS22.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS22.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS22.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS22.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS22.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end          
        
        
%% S6A_OS23
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS23_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS23_P4__HR__L2_20170305T065123_20170305T065233_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS23_P4__LRRML2_20170305T065123_20170305T065233_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS23_0001_1400RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS23_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                 
ID_GR_file_string = 'S6A_OS23';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS23=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);       

column_xls='G';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS23.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS23.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS23.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS23.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS23.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS23.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end          
            
        
%% S6A_OS24
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS24_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS24_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS24_P4__LRRML2_20170305T065123_20170305T065238_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS24_0001_1400RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS24_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                   
ID_GR_file_string = 'S6A_OS24';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS24=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0); 
        
column_xls='H';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS24.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS24.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS24.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS24.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS24.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS24.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end                 
        
%% S6A_OS30
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS30_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS30_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS30_P4__LRRML2_20170305T065123_20170305T065238_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS30_0001_1500RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS30_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                  
ID_GR_file_string = 'S6A_OS30';
ref_SSH    = 12; % m
ref_SWH    = 5; % m
ref_sigma0 = 11; % dB
S6A_OS30=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);
        
column_xls='I';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS30.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS30.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS30.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS30.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS30.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS30.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end       
        
%% S6A_OS40
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS40_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS40_P4__HR__L2_20170305T065123_20170305T065238_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS40_P4__LRRML2_20170305T065123_20170305T065238_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS40_0001_1500RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS40_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                                   
ID_GR_file_string = 'S6A_OS40';
ref_SSH    = 12; % m
ref_SWH    = 8; % m
ref_sigma0 = 11; % dB
S6A_OS40=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);    
        
column_xls='J';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS40.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS40.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS40.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS40.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS40.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS40.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end

%% S6A_OS41
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,...
                    {'HR_RAW/data/S6A_OS41_P4__HR__L2_20170305T070423_20170305T070538_0001_isd.nc',...
                     'HR_RMC/data/S6A_OS41_P4__HR__L2_20170305T070423_20170305T070538_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS41_P4__LRRML2_20170305T070423_20170305T070538_0001_isd.nc'});
%input_files(3)=strcat(common_path_bis,{'S6A_OS41_0001_1500RC_2017/LR_RMC_FROM_RAW_AD/data/S6A_OS41_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});                     
ID_GR_file_string = 'S6A_OS41';
ref_SSH    = 12; % m
ref_SWH    = 8; % m
ref_sigma0 = 11; % dB
S6A_OS41=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);    
        
column_xls='K';
for i_baseline=1:length(ID_baselines)
    % std
    xlswrite(file_xls,S6A_OS41.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
    xlswrite(file_xls,S6A_OS41.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
    xlswrite(file_xls,S6A_OS41.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
    
    % bias
    xlswrite(file_xls,abs(S6A_OS41.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
    xlswrite(file_xls,abs(S6A_OS41.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
    xlswrite(file_xls,abs(S6A_OS41.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
end
        
              