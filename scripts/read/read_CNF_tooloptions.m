function [cnf_p]=read_CNF_tooloptions(cnf_file)
    %use global variable structure instead of return data
    %global cnf_p
    %global reduced_set_cnf
    struct=loadjson(cnf_file);
    
    %% -------------------- PLOTING OPTIONS -----------------------------------
    cnf_p.mission    =   struct.mission.value;              % 0 - not plot; 1 - plot fit results
    cnf_p.name_bs = struct.name_bs.value;%[-91,91]; %latitude range of waveforms fit to be plotted (by default between -91 and 91 to force all of them to be plotted)
    cnf_p.run_L2 =struct.run_L2.value;
%     cnf_p.visible_figures=struct.visible_figures.value;
%     if ~exist('struct.text_interpreter')
%         struct.text_interpreter = 'latex';
%     end
%     cnf_p.text_interpreter=struct.text_interpreter;
    cnf_p.L2proc=struct.L2proc.value;

    % ----------- OUTPUT PLOTS OPTIONS ----------------------------------
    % fitted waveforms
    cnf_p.figure_format=struct.figure_format.value;
    cnf_p.res_fig=struct.res_fig.value;
    cnf_p.color_bs=cell2mat(struct.color_bs.value);
    cnf_p.plot_downsampling=struct.plot_downsampling.value;
    cnf_p.text_interpreter =struct.text_interpreter.value;
    cnf_p.textbox_fontsize =struct.textbox_fontsize.value;
    cnf_p.legend_fontsize =struct.legend_fontsize.value;
    cnf_p.default_fontsize =struct.default_fontsize.value;
    cnf_p.overlay_baselines =struct.overlay_baselines.value; 
    cnf_p.default_fontname=struct.default_fontname.value;
    cnf_p.default_figuresize=struct.default_figuresize.value;
    cnf_p.default_linewidth=struct.default_linewidth.value;
    cnf_p.default_markersize=struct.default_markersize.value;
    cnf_p.save_figure_format_fig=struct.save_figure_format_fig.value;
        
    % retrieved parameters performance analysis
    cnf_p.LineStyle=struct.LineStyle.value;
    cnf_p.marker_bs=struct.marker_bs.value;
    cnf_p.win_size_detrending=struct.win_size_detrending.value;
    cnf_p.flag_outliers_removal=struct.flag_outliers_removal.value;
    cnf_p.type_outliers_removal=struct.type_outliers_removal.value;
    cnf_p.smooth_param=struct.smooth_param.value;
    cnf_p.generate_plots=struct.generate_plots.value;
    cnf_p.generate_plot_SSH=struct.generate_plot_SSH.value;
    cnf_p.generate_plot_SWH=struct.generate_plot_SWH.value;
    cnf_p.generate_plot_sigma0=struct.generate_plot_sigma0.value;
    cnf_p.generate_plot_COR=struct.generate_plot_COR.value;
    % cnf_p.generate_plot_Misfit=struct.1.value;
    cnf_p.filter_ISR_baselines_mask=struct.filter_ISR_baselines_mask.value;
    cnf_p.performance_latitude_range=struct.performance_latitude_range.value;
    cnf_p.plot_fitted_waveforms_bs=struct.plot_fitted_waveforms_bs.value;
    
    % paths
    cnf_p.input_path_L1_ISR_bs=struct.input_path_L1_ISR_bs.value;
    cnf_p.input_path_L2_ISR_bs=struct.input_path_L2_ISR_bs.value;
    cnf_p.output_data_path=struct.output_data_path.value;
    cnf_p.path_to_L2_processor=struct.path_to_L2_processor.value;
    cnf_p.cnf_chd_cst_path=struct.cnf_chd_cst_path.value;

    clear struct;

end
