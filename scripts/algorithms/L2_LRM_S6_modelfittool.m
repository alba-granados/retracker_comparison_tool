function L2_LRM_S6_modelfittool(input_path_L1B,output_path_L2, cnf_chd_cst_path, cnf_tool, varargin)

% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This is the main file to run isardSAT's LRM retracker for model fit tool WP1700 within the S6 GPP
% P4 L1 project.
%
% -------------------------------------------------------------------------
% 
% Author:           Alba Granados / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Alba Granados / isardSAT V1 07/09/2020
% This software is built within the Sentinel-6 P4 L1 GPP project - CCN 3 - WP 1700
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -input_path_L1B    =   path containing the L1B file/s to be
%       processed
%       -output_path_L2 = output folder to save the different L2 products
%       (will be written whenever the corresponding flag cnf_p.write_output (in the cnf_file.m) is activated )
%       -cnf_chd_cst_path: path containing the different configuration
%       cnf_file.m
%       (configuration processing parameters for L2), characterization
%       chd_file.m
%       (different missions parameters) and the constant definition
%       cst_file.m; the Look Up Tables (LUTs) for the f0 and f1 functions
%       to be used by the single look waveform model are also included in this folder
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - output_data_generation_L2_LR: write output L2 data in netcdf
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% Code based on /retracker/algorithms/LRM/Main_S6_LRM_RTRK_LR.m by Cristina Martin-Puig / isardSAT 
% 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:

                                    
%% --------------- HANDLING OPTIONAL INPUTS ---------------------------
% if(nargin<13)
%     error('Wrong number of input parameters');
% end

p = inputParser;
p.addParamValue('MODE','LRM',@(x)ischar(x));  %LR, LROS_RAW, LROS_RMC
%p.addParamValue('SCENARIO',20);
p.addParamValue('BAND','ku',@(x)ischar(x));
p.addParamValue('SAVE',1);
p.addParamValue('filter_input_file','.nc',@(x)ischar(x));
p.addParamValue('proc_bsln_id',{''},@(x)ischar(x));

p.parse(varargin{:});
MODE = p.Results.MODE;
%SCENARIO = p.Results.SCENARIO;
BAND = p.Results.BAND;
SAVE = p.Results.SAVE;
filter_input_file = p.Results.filter_input_file;
proc_bsln_id = p.Results.proc_bsln_id;
clear p;

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

colors = cnf_tool.color_bs;

% global FILTERING
% 
% FILTERING = 0;


%% --------- Define Paths -------------------------------------------------
filesBulk.inputPath       =   input_path_L1B;
filesBulk.resultPath      =   output_path_L2;
if ~exist(filesBulk.resultPath, 'dir')
    mkdir(filesBulk.resultPath);
end
if ~exist([filesBulk.resultPath 'data' filesep], 'dir')
    mkdir([filesBulk.resultPath 'data' filesep]);
    mkdir([filesBulk.resultPath 'plots' filesep]);
    mkdir([filesBulk.resultPath 'plots' filesep 'fitted_waveforms' filesep]);    
end

path_Results = filesBulk.resultPath;


%% --------------- Configuration--------------------
inputFiles      =   dir(cnf_chd_cst_path);
aux=struct2cell(inputFiles); aux=aux(1,:); 
if ~isempty(proc_bsln_id)
    filesBulk.CNF_file=[cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file','_',proc_bsln_id])))).name];
else
    filesBulk.CNF_file=[cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file'])))).name];
end
%new read configuration files based on json files
reduced_set_cnf.SCOOP=0;
reduced_set_cnf.SHAPE=0;
cnf_p=read_CNF_json(filesBulk.CNF_file,reduced_set_cnf);

filesBulk.CST_file=[cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cst_file','.json'])))).name];
cst_p=read_CST_json(filesBulk.CST_file); %CST --> provide the global constant variables

%% Editing

cnf_p.retracker_name = {'LR'};
cnf_p.mode = MODE;

discard_records_begin=cnf_tool.discard_records_begin;

% Thermal noise or noise floor estimation window configuration
switch MODE
    case {'LROS-RMC', 'LROS-RAW'}
        cnf_p.TNini             =   25; % window positioning initial range bin for thermal noise estimation
        cnf_p.TNlast            =   35; % window positioning last range bin for thermal noise estimation
    case {'LRM'}
        cnf_p.TNini             =   10; % window positioning initial range bin for thermal noise estimation
        cnf_p.TNlast            =   20; % window positioning last range bin for thermal noise estimation
end

cnf_p.ini_Epoch             =   50; % Alba: it is better 50 as an approx. 110;    % Epoch initial values --> this can be estimated from the max peak detection or leave it like this
cnf_p.ini_Hs                =   2; %1;  % SWH initial value
cnf_p.ini_Pu                =   1;    % Initial amplitude
% % note we assume xi = misspointing = to zero , otherwise we would have to
% % input it as additional model parameter
cnf_p.plot_fits_flag              =   1; % 1 - yes, 0 -no

%% Main path

cnf_p.MainPath          = filesBulk.inputPath; %['/media/alba/DATA/isardSAT/coding/data/S6_GPP/retracker/S6A_OS20/inputs/L1B_inputs_data_Euribia/LR-RMC/'];

%% Characterisation

nf_p.semi_major_axis    =   6378137;            % Equatorial radius [m] ['float32']
nf_p.semi_minor_axis    =   6356752.3142;
nf_p.c                  =   299792458;
if strcmp(BAND,'c')
    nf_p.beamwidth      =   (3.4*pi/180);
else
    nf_p.beamwidth      =   (1.35*pi/180);
end
nf_p.BWsampl            =   395e6;
nf_p.BWClock            =   320e6;

nf_p.pi                 =   3.141592653589790;  % pi_cst from the CST file


%% ------------------ Input L1B files filtering ---------------------------
% based on /retracker/L2_bulk_processsing.m
inputFiles      =   dir(filesBulk.inputPath);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
          
filter=filter_input_file;

filterDATAFILES=(~cellfun(@isempty,strfind(lower(aux),filter)));
indexFilesL1B=find(filterDATAFILES);
    
filesBulk.nFilesL1B=length(indexFilesL1B);
filesBulk.L1BFiles=inputFiles(indexFilesL1B);

fprintf('\nTotal number of LR-L1B files to be processed: %d\n', filesBulk.nFilesL1B);


%% ------------- Loop per file to be processed ------------------------
for i_fileL1B_input=1:filesBulk.nFilesL1B

    
    %% DEFINITION OF STRUCTURE
    % -----------------------------------------------------------------
    fit_params_ini          =   [cnf_p.ini_Epoch cnf_p.ini_Hs cnf_p.ini_Pu];
    

    %% DATA STRUCTURE
    
    cnf_p.Filename = filesBulk.L1BFiles(i_fileL1B_input).name;
    
    data = readL1B_S6_ISD([cnf_p.MainPath,cnf_p.Filename],cnf_p,cst_p); 
    
    [~,file_id,fileext_L1B]=fileparts(cnf_p.Filename);
    disp(strcat('Processing L1B: ',file_id,fileext_L1B));
    
    if strcmp(MODE,'LRM')
        cnf_p.ZP            =   1;
    else
        cnf_p.ZP            =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'/global/ku/range_oversampling_factor'));  % specify if Zero Padding
    end

    ncid=netcdf.open([cnf_p.MainPath,cnf_p.Filename],'NC_NOWRITE'); % alba, from /retracker/
%         dimid=netcdf.inqDimID(ncid,'ns');
    dimid=netcdf.inqDimID(ncid,'samples');
%     dimid = 1;
    [dimname, N_samples] = netcdf.inqDim(ncid,dimid);
    if strcmp(dimname, 'samples') % new
        altitude_name = '/altitude';
    elseif strcmp(dimname, 'ns') % old
        altitude_name = '/com_altitude';
    end

    
    netcdf.close(ncid);

    nf_p.rt                 =   1/(nf_p.BWClock*cnf_p.ZP(1));
%     nf_p.Ns                 =   N_samples*cnf_p.ZP;
    nf_p.Ns                 =   N_samples; % Alba: N_samples alread oversampled
    % nf_p.Ns                 =   256*cnf_p.ZP;
    cnf_p.ZP_sample         =   nf_p.Ns/2 + 1; % Mid window location to compensate Epoch and range
    
    % -----------------------------------------------------------------
    data.GEO.TAI.total             =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/time_tai']));
    data.GEO.LAT        =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/latitude']));
    data.GEO.LON        =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/longitude']));
    data.GEO.H      	=   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND, altitude_name]));
    data.MEA.win_delay         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/tracker_range_calibrated']));

    data.Re         =   sqrt(nf_p.semi_major_axis^2*cos(data.GEO.LAT).^2+nf_p.semi_minor_axis^2*sin(data.GEO.LAT).^2);
    data.GEO.pitch      =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/off_nadir_pitch_angle_pf']));
    data.xi         =   data.GEO.pitch*pi/180; %mispointing
    data.alpha_c    =   1 + (data.GEO.H)./data.Re;
    data.h          =   data.GEO.H.*data.alpha_c;

    switch MODE
        case {'LRM'}
            wfm_init    = ones(1,length(data.GEO.LAT));
            wfm_last    = ones(1,length(data.GEO.LAT)) * nf_p.Ns;
            %         wfm_lrm         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/power_waveform']));
            i2q2_meas         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/power_waveform'])).';
        case {'LROS-RAW', 'LROS-RMC'}
            wfm_init    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_os_first_valid_sample']))+1;
            wfm_last    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_os_last_valid_sample']))+1;
            i2q2_meas         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/power_waveform'])).';
        case {'RMC'}
            wfm_init    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_os_first_valid_sample']))+1;
            wfm_last    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_os_last_valid_sample']))+1;
            wfm_lrm         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_rmc_power_waveform']));
%     i2q2_meas         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/lr_rmc_power_waveform']));
    end
    
    wfm_sc_ft       =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/waveform_scale_factor']));
    wfm_lrm=(i2q2_meas.*repmat(wfm_sc_ft,1,nf_p.Ns)).';     %apply scaling factor to waveforms
    data.HRM.s0_sf    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],['/data_20/',BAND,'/sig0_scaling_factor']));
    
    N_records = length(data.GEO.LAT);
    data.N_records = N_records;
    
    %t           =   1:cnf_p.Nfit_lim;
    options     =   optimset('Algorithm', 'levenberg-marquardt','Display','off');

    out.Epoch       =   zeros(1,length(data.GEO.LAT));
    out.Hs          =   zeros(1,length(data.GEO.LAT));
    out.Pu          =   zeros(1,length(data.GEO.LAT));
    out.xi          =   zeros(1,length(data.GEO.LAT));
    out.flag        =   zeros(1,length(data.GEO.LAT));
    out.cor         =   zeros(1,length(data.GEO.LAT));
    out.SSH         =   zeros(1,length(data.GEO.LAT));
    out.sigma0      =   zeros(1,length(data.GEO.LAT));
    wfm_lrm_norm    =   zeros(nf_p.Ns,length(data.GEO.LAT));

    for m = 1:N_records

        wfm_lrm_norm(:,m) = wfm_lrm(:,m)/max(wfm_lrm(:,m));

    %     disp(['Fitting wfm # = ',num2str(m)]);
        nf_p.xi         =   data.xi(m);
        nf_p.h          =   data.h(m);
    %     index           =   find(wfm_lrm(:,m) >= 1.3*mean(wfm_lrm(1:10,m)));
        switch MODE
            case {'RAW'} % Alba: ?
                noise_th(m)     =   2*mean(wfm_lrm_norm(cnf_p.TNini*cnf_p.ZP:cnf_p.TNlast*cnf_p.ZP,m));
                le_index        =   find(wfm_lrm_norm(1:nf_p.Ns/2,m) >= noise_th(m));
            case {'RMC'} % Alba: ?
                noise_th(m)     =   2*mean(wfm_lrm_norm(cnf_p.TNini*cnf_p.ZP:cnf_p.TNlast*cnf_p.ZP,m));
                le_index        =   find(wfm_lrm_norm(1:nf_p.Ns/2,m) >= noise_th(m));
            otherwise
                noise_th(m)     =   1.4*mean(wfm_lrm_norm(cnf_p.TNini*cnf_p.ZP:cnf_p.TNlast*cnf_p.ZP,m));
                le_index        =   find(wfm_lrm_norm(:,m) >= noise_th(m));
        end
        if (~isempty(find(diff(le_index)>1,1)))
            le_index = le_index(find(diff(le_index)>1)+1:length(le_index)); %to make sure there are not some higher samples at the beginning
                                                                %of the waveform and we have an index vector like [1,2,104,105,...]
        end

        le_index           =   le_index(1);
        

    %%  Waveform fitting

    %     index = le_index - 20;

        index = le_index;
        cnf_p.Nfit_width    =   wfm_last(m) - wfm_init(m) + 1;

        fit_wfm             =   wfm_lrm_norm(wfm_init(m):wfm_last(m),m).';


        if strcmp(MODE,'RMC') || strcmp(MODE,'RAW')
            nf_p.TN         =   noise_th(m)/2;
        else
            nf_p.TN         =   noise_th(m)/1.4;
        end
        %     nf_p.TN         =   mean(fit_wfm(cnf_p.TNini*cnf_p.ZP:cnf_p.TNlast*cnf_p.ZP));


        mpfun               =   @(fit_params,t)sl_lrm_wave_gen_LA(t, fit_params, nf_p);

        [fit_params,~,res,flagfit]     =   lsqcurvefit(mpfun,fit_params_ini,(1:length(fit_wfm)),fit_wfm,[],[],options);          

        if m < discard_records_begin
           fit_params = [NaN, NaN, NaN];  
        else        
            fit_params_ini      = fit_params;
        end
        out.Epoch(m)        = fit_params(1) + wfm_init(m) - 1;
    %     out.Epoch(m)    =   fit_params(1)*395/320;
        out.Hs(m)           =   fit_params(2) * nf_p.BWClock/nf_p.BWsampl;
    %     out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(1:index+cnf_p.Nfit_width*cnf_p.ZP,m))*wfm_sc_ft(m)); %when starting at 1
%         out.Pu(m)           =   10*log10(fit_params(3)*max(wfm_lrm(wfm_init:wfm_last,m))*wfm_sc_ft(m)); %when starting at 1, ending at cnf_p.Nfit_widht
%         out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(wfm_init:wfm_last,m))*wfm_sc_ft(m)); %when starting at index
        out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(wfm_init:wfm_last,m))); %wfm_lrm, in this case, has the scaling factor already applied.
        out.sigma0(m)       =   out.Pu(m)+data.HRM.s0_sf(m);


        out.flag(m)         =   flagfit;
        sl_wfm              =   sl_lrm_wave_gen_LA((1:length(fit_wfm)),fit_params,nf_p);
        correlation_fit     =   corrcoef(fit_wfm,sl_wfm);
        out.COR(m)          =   correlation_fit(1,2);

        if (mod(m,cnf_tool.plot_downsampling)==0) || (m==1)
                        
%             mida = get(0,'ScreenSize');
%             mida(3:4)=cnf_tool.default_figuresize;
%             set(0,'defaultFigurePaperUnits','points'); set(0,'defaultFigurePosition',mida);
%             set(groot,'defaultAxesFontSize', cnf_tool.default_fontsize); set(groot,'defaultTextFontSize', cnf_tool.default_fontsize);
%             set(groot, 'defaultAxesTickLabelInterpreter',cnf_tool.text_interpreter); set(groot, 'defaultLegendInterpreter',cnf_tool.text_interpreter);
%             set(0,'defaultFigurePaperPosition', mida);
            
            f1 = figure;
            plt=plot(wfm_lrm_norm(:,m),'Color', colors(1,:), 'LineStyle', '-');
            plt.Color(4) = 0.3; % transparency
            hold on
            plot(wfm_init(m):wfm_last(m),sl_wfm, 'Color', colors(1,:), 'LineStyle', '-')
            h_leg=legend('L1b-Waveform', 'Analytical fit','Location','northeastoutside');
            pos_leg=get(h_leg,'Position'); 
            
            grid on
            xlabel('range bin','Interpreter',cnf_tool.text_interpreter);
            if strcmp(cnf_tool.text_interpreter, 'latex')
                title(sprintf('wav. \\# %d (LAT: %.4g deg)', m, data.GEO.LAT(m)), 'Interpreter',cnf_tool.text_interpreter); 
                textbox_string = {['Epoch = ', num2str(out.Epoch(m),4), ' [r.b]'],['SWH = ' ,num2str(abs(out.Hs(m)),4), ' [m]'],...
                                    ['Pu = ', num2str(out.Pu(m),4), ' [dB]'], ['$\rho$ = ', num2str(out.COR(m)*100,5), ' [\%]'], ...
                                    ['$sigma^0$ = ', num2str(out.Pu(m)+data.HRM.s0_sf(m)), ' [dB]'] };
            else
                title(sprintf('wav. # %d (LAT: %.4g deg)', m, data.GEO.LAT(m)), 'Interpreter',cnf_tool.text_interpreter); 
                textbox_string = {['Epoch = ', num2str(out.Epoch(m),4), ' [r.b]'],['SWH = ' ,num2str(abs(out.Hs(m)),4), ' [m]'],...
                                    ['Pu = ', num2str(out.Pu(m),4), ' [dB]'], ['\rho = ', num2str(out.COR(m)*100,5), ' [%]'], ...
                                    ['\sigma^0 = ', num2str(out.Pu(m)+data.HRM.s0_sf(m)), ' [dB]'] };
            end
            annotation('textbox', [pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],...
                                'String',textbox_string,...
                                'FontSize',cnf_tool.default_fontsize,'FitBoxToText','on',  'Interpreter',cnf_tool.text_interpreter);
            axis([1 size(wfm_lrm_norm,1) 0 1])
            
            addlogoisardSAT('plot');
            if SAVE
                print(print_file,cnf_tool.res_fig,[path_Results,'plots',filesep,'fitted_waveforms',filesep,file_id,'_W',num2str(m,'%04.0f'), file_ext])
                if cnf_tool.save_figure_format_fig
                   savefig([path_Results,'plots',filesep,'fitted_waveforms',filesep,file_id,'_W',num2str(m,'%04.0f'), '.fig']) 
                end
            end
            hold off

            close(f1);
        end
            
    end
    

    %% PLOTS SSH,Hs,Pu,Sigma0
    out.retracking_cor  = (out.Epoch'-cnf_p.ZP_sample) * nf_p.rt*nf_p.BWClock/nf_p.BWsampl * nf_p.c * 0.5;
    out.SSH             = data.GEO.H - data.MEA.win_delay - out.retracking_cor;
    % Filtering outliers
    TH          = 1.5;
    indexs      = ((abs(out.Hs) < (nanmean(abs(out.Hs)) + TH*nanmean(abs(out.Hs))))  &  (abs(out.Hs) > (nanmean(abs(out.Hs)) - TH*nanmean(abs(out.Hs)))));
%     % indexs      = out.COR > .9925;
%     % indexs      = logical(ones(1,length(data.GEO.LAT)));
%     % ssh_aux = abs(out.SSH(indexs));
%     % indexs      = ssh_aux > 11;
% %     out.SSH     = out.SSH(indexs);
%     out.SSH(indexs)     = NaN;
% %     out.SSH = (com_altitude_ku-altimeter_range_calibrated_ku - (out.Epoch-cnf_p.ZP_sample+index-1)/cnf_p.ZP*nf_p.c*0.5/nf_p.BWClock).';
% %     out.Hs      = abs(out.Hs(indexs));
    out.Hs      = abs(out.Hs);
%     out.Hs(indexs)      = NaN;
%     out.sigma0(indexs)  = NaN;
% %     data.GEO.LAT    = data.GEO.LAT(indexs);

    
%% -------------- ORGANIZE THE OUTPUT DATA --------------------------------
%--------------------------------------------------------------------------
%generating an output data structure and writing the corresponding product
%(either using .nc or .mat according to the configuration)
    
%     alba: from /retracker/algorithms/L2_processing.m :
retracker_results.LR = out;

file.filename_L1B_nopath=file_id;
file.fileext_L1B=fileext_L1B;
file.inputPath=filesBulk.inputPath;
file.resultPath=filesBulk.resultPath;
fprintf('Generating output L2 product...\n')
output_data_generation_L2_LRM(file,retracker_results,data,cnf_p,cst_p);
fprintf('Done\n\n')

    
end



