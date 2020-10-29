function []=call_performance_global(input_path_L2_ISR_bs,name_bs,path_comparison_results,i_proc,varargin)

    %% --------------- HANDLING OPTIONAL INPUTS ---------------------------
    if(nargin<3 || nargin>(3+44*2))
        error('Wrong number of input parameters');
    end
    p = inputParser;
    p.addParamValue('figure_format','',@(x)ischar(x));
    p.addParamValue('res_fig','',@(x)ischar(x));
    p.addParamValue('win_size_detrending',20);
    p.addParamValue('step_SWH',0.2);
    p.addParamValue('step_SSH',0.5);
    p.addParamValue('step_hr',0.5);
    p.addParamValue('sh_name_nc','ssh');
    p.addParamValue('flag_outliers_removal',0);
    p.addParamValue('type_outliers_removal','',@(x)ischar(x));
    p.addParamValue('smooth_param',0);
    p.addParamValue('define_min_max_SWH',0);
    p.addParamValue('min_SWH',1.5);
    p.addParamValue('max_SWH',4.5);
    p.addParamValue('input_path_L1_ISR_bs',{''},@(x)iscellstr(x));
    p.addParamValue('input_path_L2_ESA',{''},@(x)iscellstr(x));
    p.addParamValue('input_path_L2_GPOD',{''},@(x)iscellstr(x));
    p.addParamValue('name_bs_L2_ESA',{'L2 ESA'},@(x)iscellstr(x));
    p.addParamValue('name_bs_L2_GPOD',{'L2 GPOD'},@(x)iscellstr(x));
    p.addParamValue('input_path_L2_ISR_L1B_ESA',{''},@(x)iscellstr(x));
    p.addParamValue('name_bs_L2_ISR_L1B_ESA',{'L1B-ESA'},@(x)iscellstr(x));
    p.addParamValue('input_path_L2_STL',{''},@(x)iscellstr(x));
    p.addParamValue('name_bs_L2_STL',{'SCOOP STL'},@(x)iscellstr(x));
    p.addParamValue('input_path_L1_ESA',{''},@(x)iscellstr(x));
    p.addParamValue('annotation_box_active',1);
    p.addParamValue('filename_mask_KML',{''},@(x)iscellstr(x));
    p.addParamValue('num_pools',1);
    p.addParamValue('plot_downsampling',50);
    p.addParamValue('generate_plots',0);
    p.addParamValue('generate_plot_SSH',1);
    p.addParamValue('generate_plot_SWH',1);
    p.addParamValue('generate_plot_sigma0',1);
    p.addParamValue('generate_plot_nb',1);
    p.addParamValue('generate_plot_COR',1);
    p.addParamValue('generate_plot_misspointing',1);
    p.addParamValue('generate_kml',0);
    p.addParamValue('filter_land_surf_type',0);
    p.addParamValue('compute_comparison_performance_tracks',0);
    p.addParamValue('filter_ISR_baselines_mask',0);
    %inclusion of DEM
    p.addParamValue('DEM_inclusion',0); %incorporate the height of surface w.r.t elliposid in the different plots
    p.addParamValue('DEM_ref','SRTM',@(x)ischar(x)); %indicate type of DEM used: currently only the SRTM is available
    p.addParamValue('DEM_dir','.',@(x)ischar(x)); %directory where to save the downloaded DEM for SRTM
    p.addParamValue('mission','CS2',@(x)ischar(x)); %directory where to save the downloaded DEM for SRTM
    p.addParamValue('ref_SSH',12); %reference SSH to compute error for S6 sim
    p.addParamValue('ref_SWH',2); %reference SWH to compute error for S6 sim
    p.addParamValue('ref_sigma0',12); %reference sigma0 to compute error for S6
    
    p.parse(varargin{:});
    figure_format=p.Results.figure_format;
    res_fig=p.Results.res_fig;
    win_size_detrending=p.Results.win_size_detrending;
    step_SWH=p.Results.step_SWH;
    step_SSH=p.Results.step_SSH;
    step_hr=p.Results.step_hr;
    sh_name_nc=p.Results.sh_name_nc;
    flag_outliers_removal=p.Results.flag_outliers_removal;
    type_outliers_removal=p.Results.type_outliers_removal;
    smooth_param=p.Results.smooth_param;
    define_min_max_SWH=p.Results.define_min_max_SWH;
    min_SWH=p.Results.min_SWH;
    max_SWH=p.Results.max_SWH;
    input_path_L1_ISR_bs=p.Results.input_path_L1_ISR_bs;
    input_path_L2_ESA=p.Results.input_path_L2_ESA;
    input_path_L2_GPOD=p.Results.input_path_L2_GPOD;
    input_path_L2_ISR_L1B_ESA=p.Results.input_path_L2_ISR_L1B_ESA;
    name_bs_L2_ISR_L1B_ESA=p.Results.name_bs_L2_ISR_L1B_ESA;
    input_path_L1_ESA=p.Results.input_path_L1_ESA;
    input_path_L2_STL=p.Results.input_path_L2_STL;
    annotation_box_active=p.Results.annotation_box_active;
    filename_mask_KML=p.Results.filename_mask_KML;
    num_pools=p.Results.num_pools;
    plot_downsampling=p.Results.plot_downsampling;
    generate_plots=p.Results.generate_plots;
    generate_plot_SSH=p.Results.generate_plot_SSH;
    generate_plot_SWH=p.Results.generate_plot_SWH;
    generate_plot_sigma0=p.Results.generate_plot_sigma0;
    generate_plot_nb=p.Results.generate_plot_nb;
    generate_plot_COR=p.Results.generate_plot_COR;
    generate_plot_misspointing=p.Results.generate_plot_misspointing;
    generate_kml=p.Results.generate_kml;
    filter_land_surf_type=p.Results.filter_land_surf_type;
    name_bs_L2_STL=p.Results.name_bs_L2_STL;
    name_bs_L2_ESA=p.Results.name_bs_L2_ESA;
    name_bs_L2_GPOD=p.Results.name_bs_L2_GPOD;
    compute_comparison_performance_tracks=p.Results.compute_comparison_performance_tracks;
    DEM_inclusion=p.Results.DEM_inclusion;
    DEM_ref=p.Results.DEM_ref;
    DEM_dir=p.Results.DEM_dir;
    filter_ISR_baselines_mask=p.Results.filter_ISR_baselines_mask;
    mission=p.Results.mission;
    ref_SSH=p.Results.ref_SSH;
    ref_SWH=p.Results.ref_SWH;
    ref_sigma0=p.Results.ref_sigma0;
    clear p

    %% -------------------------- CALL PERFORMANCE VALIDATION -------------
    
    performance_global_ISR_baselines_paralelization_new_GPOD...
        (input_path_L2_ISR_bs(i_proc,:),...
        name_bs(i_proc,:),char(path_comparison_results(i_proc)),...
        'figure_format',figure_format,...
        'res_fig',res_fig,...
        'win_size_detrending',win_size_detrending,...
        'step_SWH',step_SWH,...
        'sh_name_nc',sh_name_nc,...
        'flag_outliers_removal',flag_outliers_removal,...
        'type_outliers_removal',type_outliers_removal,...
        'smooth_param',smooth_param,...
        'input_path_L2_ESA',char(input_path_L2_ESA(i_proc)),...
        'input_path_L1_ISR_bs',input_path_L1_ISR_bs(i_proc,:),...
        'input_path_L1_ESA',char(input_path_L1_ESA(i_proc)),...
        'input_path_L2_STL',char(input_path_L2_STL(i_proc)),...
        'input_path_L2_GPOD',char(input_path_L2_GPOD(i_proc)),...
        'filename_mask_KML',char(filename_mask_KML(i_proc)),...
        'num_pools',num_pools,...
        'generate_plots',generate_plots,...
        'plot_downsampling',plot_downsampling,...
        'define_min_max_SWH',define_min_max_SWH,...
        'min_SWH',min_SWH,...
        'max_SWH',max_SWH,...
        'generate_kml',generate_kml,...
        'filter_land_surf_type',filter_land_surf_type,...
        'compute_comparison_performance_tracks',compute_comparison_performance_tracks,...
        'generate_plot_SSH',generate_plot_SSH,'generate_plot_SWH',generate_plot_SWH,...
        'generate_plot_sigma0',generate_plot_sigma0,'generate_plot_nb',generate_plot_nb,...
        'generate_plot_COR',generate_plot_COR,...
        'generate_plot_misspointing',generate_plot_misspointing,...
        'filter_ISR_baselines_mask',filter_ISR_baselines_mask,...
        'mission',mission,...
        'ref_SSH',ref_SSH(i_proc),...
        'ref_SWH',ref_SWH(i_proc),...
        'ref_sigma0',ref_sigma0(i_proc));

end