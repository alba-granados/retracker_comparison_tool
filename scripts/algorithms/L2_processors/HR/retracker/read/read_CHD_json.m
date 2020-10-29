%function to read the characterization file in JSON and adapted to MAT
%inputs
function [chd_p]=read_CHD_json(chd_file,cnf_p,cst_p,varargin)

    %load external values of the PTRs
    p = inputParser;
    p.addParamValue('PTR_along_external',[]); %reference SSH to compute error for S6 sim
    p.addParamValue('PTR_across_external',[]); %reference SWH to compute error for S6 sim
    
    p.parse(varargin{:});
    PTR_along_external       = p.Results.PTR_along_external;
    PTR_across_external       = p.Results.PTR_across_external;
    clear p;


    struct=loadjson(chd_file);

    % -----------------------------------------------------------------
    % OPERATION CONFIGURATION
    % -----------------------------------------------------------------
    % Total number of pulses in a burst this is equal for SAR and SARin
    % mode (uint8)
    chd_p.N_total_pulses_b_chd    =   struct.N_total_pulses_b_chd.value;
    % Number of samples per each SAR or SARin pulse (uint8)
    chd_p.N_samples_chd       =   struct.N_samples_chd.value;
    %Number of burst in a radar cycle
    % Number of burst in a SARM tracking cycle ('uint8')
    chd_p.N_bursts_cycle_chd  =   struct.N_bursts_cycle_chd.value;
    
    % bias in pitch and roll 
    %(BE CAREFUL WITH SUCH PARAMETERS: IF WE READ DIRECTLY L1B ESA OR L1B ISD GENERATED WITH FBR BASELINE-C SUCH BIAS ARE ALREADY CORRECTED)
    if cnf_p.SCOOP_flag==0 && cnf_p.SHAPE_flag==0
        chd_p.pitch_bias_chd                  =   struct.pitch_bias_chd.value;
    else
        chd_p.pitch_bias_chd                  =   0.0;
    end

    % Roll bias
    if cnf_p.SCOOP_flag==0 && cnf_p.SHAPE_flag==0
        chd_p.roll_bias_chd                   =   struct.roll_bias_chd.value;
    else
        chd_p.roll_bias_chd                   =   0.0;
    end

    %PTR width
    if cnf_p.SCOOP_flag==0 && cnf_p.SHAPE_flag==0
        chd_p.PTR_width_chd=struct.PTR_width_chd.value;
    else
        chd_p.PTR_width_chd=2.819e-09;
    end

    %TX power
    if cnf_p.SCOOP_flag==0 && cnf_p.SHAPE_flag==0
        chd_p.power_tx_ant_ku_chd = struct.power_tx_ant_ku_chd.value;
    else
        chd_p.power_tx_ant_ku_chd = 20;
    end
    
    % -----------------------------------------------------------------
    % INSTRUMENT PARAMETERS
    % -----------------------------------------------------------------
    % Central frequency in Ku band (float64)        
    chd_p.freq_ku_chd         =   struct.freq_ku_chd.value;
    chd_p.wv_length_ku        =   cst_p.c_cst/chd_p.freq_ku_chd;
    
    % PRF FOR SAR and SARin modes (float64)
    chd_p.prf_chd             =   struct.prf_chd.value;
    chd_p.brf_chd             =   struct.brf_chd.value;
    
    
    % Ku received Bandwidth SAR mode (float64)
    chd_p.bw_rx_ku_chd        =   struct.bw_rx_ku_chd.value;
    chd_p.fs_clock_ku_chd     =   struct.fs_clock_ku_chd.value; 
    
    % -----------------------------------------------------------------
    % ANTENNA Characterization information
    % -----------------------------------------------------------------
    chd_p.antenna_beamwidth_act_ku_chd    =   struct.antenna_beamwidth_act_ku_chd.value;
    chd_p.antenna_beamwidth_alt_ku_chd    =   struct.antenna_beamwidth_alt_ku_chd.value;
    if cnf_p.SCOOP_flag==0 && cnf_p.SHAPE_flag==0
        chd_p.antenna_gain_ku_chd             = 	struct.antenna_gain_ku_chd.value; %dB
    else
        chd_p.antenna_gain_ku_chd             = 	42.6; %dB
    end
    
    
    
    %----------------------------------------------------------------------
    % PTR settings
    %----------------------------------------------------------------------
    if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA')) 
        %azimuth
        if isempty(PTR_along_external)
            switch cnf_p.window_type_a
                case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
                    chd_p.A_s2Ga_chd=1.0196;
                    chd_p.alpha_ga_chd=1.0/(2*(struct.PTR_Boxcar.value).^2);
                case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
                    chd_p.A_s2Ga_chd=1.0101;
                    chd_p.alpha_ga_chd=1.0/(2*(struct.PTR_Hanning.value).^2);
                case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
                    chd_p.A_s2Ga_chd=1.0081;
                    chd_p.alpha_ga_chd=1.0/(2.0*(struct.PTR_Hamming.value).^2);
                case 'Forced' %taken from the cnf file
                    chd_p.A_s2Ga_chd=cnf_p.A_s2Ga_chd;
                    chd_p.alpha_ga_chd=cnf_p.alpha_ga_chd;
                otherwise
                    chd_p.A_s2Ga_chd=[];
                    chd_p.alpha_ga_chd=[];
            end
        else
            chd_p.A_s2Ga_chd=1;
            chd_p.alpha_ga_chd=1.0/(2.0*(PTR_along_external).^2);            
        end
        %disp(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_ga_chd)))))
        if isempty(PTR_across_external)
            %range
            switch cnf_p.window_type_r
                case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
                    chd_p.A_s2Gr_chd=1.0196;
                    chd_p.alpha_gr_chd=1.0/(2*(struct.PTR_Boxcar.value).^2);
                case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
                    chd_p.A_s2Gr_chd=1.0101;
                    chd_p.alpha_gr_chd=1.0/(2*(struct.PTR_Hanning.value).^2);
                case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
                    chd_p.A_s2Gr_chd=1.0081;
                    chd_p.alpha_gr_chd=1.0/(2.0*(struct.PTR_Hamming.value).^2);
                case 'Forced'%taken from the cnf file
                    chd_p.A_s2Gr_chd=cnf_p.A_s2Gr_chd;
                    chd_p.alpha_gr_chd=cnf_p.alpha_gr_chd;
                otherwise
                    chd_p.A_s2Gr_chd=[];
                    chd_p.alpha_gr_chd=[];
            end
            %disp(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_gr_chd)))))
        else
            chd_p.A_s2Gr_chd=1;
            chd_p.alpha_gr_chd=1.0/(2.0*(PTR_across_external).^2);            
        end
    end
    
    clear struct;
end
