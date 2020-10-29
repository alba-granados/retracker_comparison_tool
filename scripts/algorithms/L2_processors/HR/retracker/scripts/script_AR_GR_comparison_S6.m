%% Comparison of GR for different scenarios
common_path ='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1B_GPP/'; 
resultPath  ='C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/results/L2_L1B_GPP/GR_performance/'; 
%ID_baselines={'RAW','RMC','LR-RMC'};



mkdir(resultPath);
       
file_xls=strcat(resultPath,'performance_AR.xlsx');
row_std_SSH  = '3';
row_std_sig0 = '5';
row_std_SWH  = '7';

row_bias_SSH  = '10';
row_bias_sig0 = '12';
row_bias_SWH  = '14';
        

%% S6A_OS16 (Inter-alignment)
ID_baselines={'RAW','RMC','LR-RMC'}; %
input_files = strcat(common_path,'S6A_OS16_0001_1500RC_2017/',...
                    {'RAW/data/S6A_OS10_P4__HR__L2_20170305T065123_20170305T065233_RAW1_isd.nc',...
                     'RMC/data/S6A_OS10_P4__HR__L2_20170305T065123_20170305T065233_RMC1_isd.nc',...
                     'LR_RMC/data/S6A_OS10_P4__LRRML2_20170305T065123_20170305T065233_0001_isd.nc'});
ID_GR_file_string = 'S6A_OS16';
ref_SSH_whole    = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','ground_truth_SSH_20_ku').'; % m
ref_SSH_lat      = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','latitude_20_ku').'; % deg
ref_SWH    = 1; % m
ref_sigma0 = 11; % dB
S6A_OS16=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0,...
            'ref_SSH_whole',ref_SSH_whole,'lat_SSH_whole',ref_SSH_lat,...
            'plot_figure',0); 

column_xls='L';

for i_baseline=1:length(ID_baselines)
    if strfind(char(ID_baselines(i_baseline)),'LR-RMC RAW')
        % std
        xlswrite(file_xls,S6A_OS16.SSH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SSH));
        xlswrite(file_xls,S6A_OS16.sigma0_mean_std(i_baseline),char(ID_baselines(i_baseline)),strcat(column_xls,row_std_sig0));
        xlswrite(file_xls,S6A_OS16.SWH_mean_std(i_baseline)*100.0,char(ID_baselines(i_baseline)),strcat(column_xls,row_std_SWH));
        
        % bias
        xlswrite(file_xls,abs(S6A_OS16.SSH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SSH));
        xlswrite(file_xls,abs(S6A_OS16.sigma0_mean_bias(i_baseline)),char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_sig0));
        xlswrite(file_xls,abs(S6A_OS16.SWH_mean_bias(i_baseline)).*100,char(ID_baselines(i_baseline)),strcat(column_xls,row_bias_SWH));
    end
end

%% S6A_OS10
input_files = strcat(common_path,'S6A_OS10_0001_1400RC_2017/',...
                    {'RAW/L2/data/S6A_OS10_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS10_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS10_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS10_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS16_0001_1500RC_2017/',...
                    {'RAW/L2/data/S6A_OS16_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS16_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS16_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW_ADAPTIVE_CENTROID_SLOPE/data/S6A_OS16_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
                     %'LR_RMC_FROM_RAW/data/S6A_OS16_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
ID_GR_file_string = 'S6A_OS16';
ref_SSH_whole    = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','ground_truth_SSH_20_ku').'; % m
ref_SSH_lat      = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','latitude_20_ku').'; % deg
ref_SWH    = 1; % m
ref_sigma0 = 11; % dB
S6A_OS16=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0,...
            'ref_SSH_whole',ref_SSH_whole,'lat_SSH_whole',ref_SSH_lat);   
        
column_xls='C';
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
input_files = strcat(common_path,'S6A_OS20_0001_1400RC_2017/',...
                    {'RAW/L2/data/S6A_OS20_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS20_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS20_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS20_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
ID_GR_file_string = 'S6A_OS20';
ref_SSH    = 12; % m
ref_SWH    = 2; % m
ref_sigma0 = 11; % dB
S6A_OS20=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
            'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0); 
        
        
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
input_files = strcat(common_path,'S6A_OS21_0001_2500RC_2017/',...
                    {'RAW/L2/data/S6A_OS21_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS21_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS21_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS21_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS22_0001_1400RC_2017/',...
                    {'RAW/L2/data/S6A_OS22_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS22_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS22_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS22_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS23_0001_1400RC_2017/',...
                    {'RAW/L2/data/S6A_OS23_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS23_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS23_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS23_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS24_0001_1400RC_2017/',...
                    {'RAW/L2/data/S6A_OS24_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS24_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS24_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS24_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS30_0001_1500RC_2017/',...
                    {'RAW/L2/data/S6A_OS30_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS30_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS30_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS30_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS40_0001_1500RC_2017/',...
                    {'RAW/L2/data/S6A_OS40_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS40_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS40_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS40_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
input_files = strcat(common_path,'S6A_OS41_0001_1500RC_2017/',...
                    {'RAW/L2/data/S6A_OS41_P4__RAW_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'RMC/L2/data/S6A_OS341_P4__RMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC/data/S6A_OS341_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc',...
                     'LR_RMC_FROM_RAW/data/S6A_OS41_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
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
% %% S6A_OS10
% input_files = strcat(common_path,'S6A_OS10_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS10_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T212452_isd.nc',...
%                      'RMC/L2/data/S6A_OS10_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T212536_isd.nc',...
%                      'LR_RMC/data/S6A_OS10_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234024_isd.nc'});
% ID_GR_file_string = 'S6A_OS10';
% ref_SSH    = 12; % m
% ref_SWH    = 1; % m
% ref_sigma0 = 11; % dB
% S6A_OS10=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);       
%         
% %% S6A_OS16
% input_files = strcat(common_path,'S6A_OS16_0001_1500RC_2017/',...
%                     {'RAW/L2/data/S6A_OS16_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T212736_isd.nc',...
%                      'RMC/L2/data/S6A_OS16_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T212916_isd.nc',...
%                      'LR_RMC/data/S6A_OS16_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234104_isd.nc'});
% ID_GR_file_string = 'S6A_OS16';
% ref_SSH_whole    = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','ground_truth_SSH_20_ku').'; % m
% ref_SSH_lat      = ncread('C:/Users/eduard.makhoul/isardSAT/projects/Sentinel-6/data/AR/inputs/S6A_OS16_0001_1500RC_2017/S6A_OS11_0001_slope_ground_truth_file.nc','latitude_20_ku').'; % deg
% ref_SWH    = 1; % m
% ref_sigma0 = 11; % dB
% S6A_OS16=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0,...
%             'ref_SSH_whole',ref_SSH_whole,'lat_SSH_whole',ref_SSH_lat);   
%         
%         
% %% S6A_OS20
% input_files = strcat(common_path,'S6A_OS20_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS20_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T211945_isd.nc',...
%                      'RMC/L2/data/S6A_OS20_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211802_isd.nc',...
%                      'LR_RMC/data/S6A_OS20_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234006_isd.nc'});
% ID_GR_file_string = 'S6A_OS20';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS20=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0); 
%         
%         
% %% S6A_OS20_LR_RMC_RAW_fixed
% input_files = strcat(common_path,'S6A_OS20_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS20_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T211945_isd.nc',...
%                      'RMC/L2/data/S6A_OS20_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211802_isd.nc',...
%                      'LR_RMC_FROM_RAW_fixed/data/S6A_OS20_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
% ID_GR_file_string = 'S6A_OS20_LR_RMC_RAW_fixed';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS20=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);         
%         
%         
% %% S6A_OS20_LR_RMC_RAW_adaptive
% input_files = strcat(common_path,'S6A_OS20_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS20_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T211945_isd.nc',...
%                      'RMC/L2/data/S6A_OS20_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211802_isd.nc',...
%                      'LR_RMC_FROM_RAW_adaptive/data/S6A_OS20_P4__LRRMC_L2_00000000T000000_99999999T999999_0001_isd.nc'});
% ID_GR_file_string = 'S6A_OS20_LR_RMC_RAW_adaptive';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS20=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);        
%         
% %% S6A_OS21
% input_files = strcat(common_path,'S6A_OS21_0001_2500RC_2017/',...
%                     {'RAW/L2/data/S6A_OS21_P4__RAW_L2_00000000T000000_99999999T999999_0001isd.nc',...
%                      'RMC/L2/data/S6A_OS21_P4__RMC_L2_00000000T000000_99999999T999999_0001isd.nc',...
%                      'LR_RMC/data/S6A_OS21_P4__LRRMC_L2_00000000T000000_99999999T999999_0001isd.nc'});
% ID_GR_file_string = 'S6A_OS21';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS21=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);     
% 
% %% S6A_OS22
% input_files = strcat(common_path,'S6A_OS22_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS22_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T212316_isd.nc',...
%                      'RMC/L2/data/S6A_OS22_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T212241_isd.nc',...
%                      'LR_RMC/data/S6A_OS22_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234006_isd.nc'});
% ID_GR_file_string = 'S6A_OS22';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS22=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);    
%         
%         
% %% S6A_OS23
% input_files = strcat(common_path,'S6A_OS23_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS23_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T212001_isd.nc',...
%                      'RMC/L2/data/S6A_OS23_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211958_isd.nc',...
%                      'LR_RMC/data/S6A_OS23_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T233957_isd.nc'});
% ID_GR_file_string = 'S6A_OS23';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS23=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);       
%         
% %% S6A_OS24
% input_files = strcat(common_path,'S6A_OS24_0001_1400RC_2017/',...
%                     {'RAW/L2/data/S6A_OS24_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T212548_isd.nc',...
%                      'RMC/L2/data/S6A_OS24_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T212606_isd.nc',...
%                      'LR_RMC/data/S6A_OS24_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234047_isd.nc'});
% ID_GR_file_string = 'S6A_OS24';
% ref_SSH    = 12; % m
% ref_SWH    = 2; % m
% ref_sigma0 = 11; % dB
% S6A_OS24=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0); 
%         
% %% S6A_OS30
% input_files = strcat(common_path,'S6A_OS30_0001_1500RC_2017/',...
%                     {'RAW/L2/data/S6A_OS30_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T211429_isd.nc',...
%                      'RMC/L2/data/S6A_OS30_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211418_isd.nc',...
%                      'LR_RMC/data/S6A_OS30_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234037_isd.nc'});
% ID_GR_file_string = 'S6A_OS30';
% ref_SSH    = 12; % m
% ref_SWH    = 5; % m
% ref_sigma0 = 11; % dB
% S6A_OS30=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);
%         
% %% S6A_OS40
% input_files = strcat(common_path,'S6A_OS40_0001_1500RC_2017/',...
%                     {'RAW/L2/data/S6A_OS40_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T211309_isd.nc',...
%                      'RMC/L2/data/S6A_OS40_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211303_isd.nc',...
%                      'LR_RMC/data/S6A_OS40_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234053_isd.nc'});
% ID_GR_file_string = 'S6A_OS40';
% ref_SSH    = 12; % m
% ref_SWH    = 8; % m
% ref_sigma0 = 11; % dB
% S6A_OS40=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);     
%         
% %% S6A_OS41
% input_files = strcat(common_path,'S6A_OS41_0001_1500RC_2017/',...
%                     {'RAW/L2/data/S6A_OS41_P4__L2__1B_00000000T000000_99999999T999999_0001_20180307T210931_isd.nc',...
%                      'RMC/L2/data/S6A_OS341_P4__RMC_1B_00000000T000000_99999999T999999_0001_20180307T211305_isd.nc',...
%                      'LR_RMC/data/S6A_OS341_P4__LRRMC_1B_00000000T000000_99999999T999999_0001_20180307T234059_isd.nc'});
% ID_GR_file_string = 'S6A_OS41';
% ref_SSH    = 12; % m
% ref_SWH    = 8; % m
% ref_sigma0 = 11; % dB
% S6A_OS41=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,...
%             'ref_SSH',ref_SSH,'ref_SWH',ref_SWH,'ref_sigma0',ref_sigma0);                   