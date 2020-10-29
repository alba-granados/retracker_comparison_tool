%% Comparison of TRP performance for different scenarios
%path that includes the different .mat metric files for the TRP over the different
%scenarios considered in the AR
performance_results_path ='\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Performance_S6\Transponder\inputs\example_AR_output_GPP_TRP_performance\'; 
%Excel file used to include the geophysical retrievals performacne and to be updated with the information of the TRP
file_xls = '\\N8800PRO\share\ISARDS\emak\Sentinel-6\test_data\Performance_S6\geophysical_retrievals\results\performance_AR_GPP.xlsx';

%ID_baselines={'RAW','RMC','LR-RMC'};
ID_baselines={'RAW','RMC'}; %
ID_method_processing={'exact'};



      
row_std_error_stack  = '3';
row_mean_error_stack = '5';
row_range_variation_stack = '7';
row_error_L1B  = '10';
row_datation_error = '12';


ID_GR_file_string = {'S6A_PT10','S6A_PT11','S6A_PT12'};
column_xls = {'P','Q'; 'R','S'; 'T', 'U'};
for i_scenario=1:length(ID_GR_file_string)
    for i_baseline=1:length(ID_baselines)
        for i_method=1:length(ID_method_processing)
            file=dir([performance_results_path '*' ...
                char(ID_GR_file_string(i_scenario)) '*' ...
                char(ID_baselines(i_baseline)) '*' ...
                char(ID_method_processing(i_method)) '.mat']);
            if ~isempty(file)
                load([performance_results_path char(file.name)]);
            else
                continue;
            end
            xlswrite(file_xls,abs(Req.std_error_stack),char(ID_baselines(i_baseline)),strcat(char(column_xls(i_scenario,i_method)),row_std_error_stack));
            xlswrite(file_xls,abs(Req.mean_error_stack),char(ID_baselines(i_baseline)),strcat(char(column_xls(i_scenario,i_method)),row_mean_error_stack));
            xlswrite(file_xls,abs(Req.range_variation_stack),char(ID_baselines(i_baseline)),strcat(char(column_xls(i_scenario,i_method)),row_range_variation_stack));
            xlswrite(file_xls,abs(Req.error_L1B),char(ID_baselines(i_baseline)),strcat(char(column_xls(i_scenario,i_method)),row_error_L1B));
            xlswrite(file_xls,abs(Req.datation_error),char(ID_baselines(i_baseline)),strcat(char(column_xls(i_scenario,i_method)),row_datation_error));
            
            
        end
    end
end

              