function [data]   = sigma0_scaling_factor (data,cnf_p,chd_p,cst_p)


%   > Surface area (of every beam)
norm_vel_sat = data.GEO.V;
range_sat_surf = cst_p.c_cst/2.*data.MEA.win_delay;

azimuth_distance = (1+range_sat_surf./cst_p.earth_radius_cst).*chd_p.wv_length_ku .* range_sat_surf ./ (1./chd_p.prf_chd) ./ (2 * (norm_vel_sat) .* chd_p.N_total_pulses_b_chd);

range_distance = 2 * sqrt (cst_p.c_cst * range_sat_surf * (chd_p.PTR_width_chd).*...
    cst_p.earth_radius_cst./(cst_p.earth_radius_cst + range_sat_surf) );

switch lower(cnf_p.window_type_a)
    case 'hamming'
        wf_a=1.486*0.92;%widening factor as denoted by Dinardo
    case 'hanning'
        wf_a=1.0;% TBD
    otherwise
        wf_a=1.0;
end

switch lower(cnf_p.window_type_r)
    case 'hamming'
        wf_r=1.486*0.92;%widening factor as denoted by Dinardo
    case 'hanning'
        wf_r=1.0;% TBD
    otherwise
        wf_r=1.0;
end

switch cnf_p.mission
    case {'CS2','CR2'}
        zp_correction = -10*log10(cnf_p.ZP);
    otherwise
        zp_correction = 0.0;
end

surface_area = (wf_a.*azimuth_distance .* (wf_r.*range_distance)).*0.886;

%added by EM 30.11.2015: consider the surface and not the mean over
%the different beams & compensate for norm. in fft range
%compression & TBP or pulse compression gain
data.HRM.s0_sf = 10*log10(64) + 30*log10(cst_p.pi_cst)...
    + 40*log10(range_sat_surf) - 10*log10(chd_p.power_tx_ant_ku_chd) - 2*(chd_p.antenna_gain_ku_chd)...
    - 20*log10(chd_p.wv_length_ku) - 10*log10(surface_area) + zp_correction;




    
    
    
    
    
    
    
end