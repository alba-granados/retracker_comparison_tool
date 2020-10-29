function [res] = geo_corrections_computation(data,cnf_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and geophysical correctiosn to be applied to the
% retracked range based on the info available in the L1B product and
% depending on the surface being observed
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 17/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -data    =   structure of input data used in the L2 processing       
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       res        =   structure of Corrections to be included in the L2
%       product
% RESTRICTIONS: The type dependent surface corrections for each record
% independently is implemented only for CS2 data only for ocean (falg
% surface==0) and land ice (flag surface==2). Application of same
% geophysical corrections for all the records based on a single surface
% type is currently available only for "open ocean", "land ice" and "sea
% ice". No inclusion of the sea state bias corrections is available since
% there is no such available info at L1B.
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 



% -------- include the geophysical corrections to the output data -----
res=rmfield(data.COR,{'prop_GIM_ion','prop_Mod_ion','surf_dac','surf_invb','geop'});
res.surf_type_flag = data.surf_type_flag;
res.total_geo_corr        =   zeros(1,data.N_records);

%------- corrections flags --------------------------------------------
% flags indicating the inclusion 1 or not 0 of the related correction
res.flags.dry_trop                  =  zeros(1,data.N_records);
res.flags.wet_trop                  =  zeros(1,data.N_records);
res.flags.inv_bar                   =  zeros(1,data.N_records);
res.flags.dac                       =  zeros(1,data.N_records);
res.flags.gim_ion                   =  zeros(1,data.N_records);
res.flags.model_ion                 =  zeros(1,data.N_records);
res.flags.ocean_equilibrium_tide    =  zeros(1,data.N_records);
res.flags.ocean_longperiod_tide     =  zeros(1,data.N_records);
res.flags.ocean_loading_tide        =  zeros(1,data.N_records);
res.flags.solidearth_tide           =  zeros(1,data.N_records);
res.flags.geocentric_polar_tide     =  zeros(1,data.N_records);




if cnf_p.force_geocorr_surf_type
    %---------------------total as if all where same surface --------------
    % Different contributions depending on surface extracted from CS2
    % product handbook
    % force the same group of corrections to all records
    % independently of the surface flag included in the L1B product
    
    % ---------- Checking iono corrections to be applied ------------------
    if any(res.gim_ion)
        %from product handbook of Cryosat-2 GIM is nominal choice
        %and if not available use Bent model
        iono_corr=res.gim_ion;
        res.flags.gim_ion(:)=1;
    else
        iono_corr=res.model_ion;
        res.flags.model_ion(:)=1;
    end
    
    switch cnf_p.product_type_surface
        case 'open_ocean'
            res.total_geo_corr= res.ocean_equilibrium_tide+... %is referred in cryosat-2 as ocean tide or elastic ocean tide
                res.ocean_longperiod_tide+...
                res.ocean_loading_tide+...
                res.solidearth_tide+...
                res.geocentric_polar_tide+...
                res.dry_trop+...
                res.wet_trop+...
                iono_corr+...
                res.dac; %missing sea state bias
            
            res.flags.ocean_equilibrium_tide(:)=1;
            res.flags.ocean_longperiod_tide(:)=1;
            res.flags.ocean_loading_tide(:)=1;
            res.flags.solidearth_tide(:)=1;
            res.flags.geocentric_polar_tide(:)=1;
            res.flags.dry_trop(:)=1;
            res.flags.wet_trop(:)=1;
            res.flags.dac(:)=1;
        case 'sea_ice'
            res.total_geo_corr=res.ocean_equilibrium_tide+...
                res.ocean_longperiod_tide+...
                res.ocean_loading_tide+...
                res.solidearth_tide+...
                res.geocentric_polar_tide+...
                res.dry_trop+...
                res.wet_trop+...
                iono_corr+...
                res.inv_bar;
            res.flags.ocean_equilibrium_tide(:)=1;
            res.flags.ocean_longperiod_tide(:)=1;
            res.flags.ocean_loading_tide(:)=1;
            res.flags.solidearth_tide(:)=1;
            res.flags.geocentric_polar_tide(:)=1;
            res.flags.dry_trop(:)=1;
            res.flags.wet_trop(:)=1;
            res.flags.inv_bar(:)=1;
        case 'land_ice'
            res.total_geo_corr=res.ocean_loading_tide+...
                res.solidearth_tide+...
                res.geocentric_polar_tide+...
                res.dry_trop+...
                res.wet_trop+...
                iono_corr;
            res.flags.ocean_loading_tide(:)=1;
            res.flags.solidearth_tide(:)=1;
            res.flags.geocentric_polar_tide(:)=1;
            res.flags.dry_trop(:)=1;
            res.flags.wet_trop(:)=1;
        case 'land'
            res.total_geo_corr=res.solidearth_tide+...
                res.geocentric_polar_tide+...
                res.dry_trop+...
                res.wet_trop+...
                iono_corr;
            res.flags.solidearth_tide(:)=1;
            res.flags.geocentric_polar_tide(:)=1;
            res.flags.dry_trop(:)=1;
            res.flags.wet_trop(:)=1;            
            
    end
else
    %apply the corrections according to the types of surfaces
    %provided in the L1B product itself
    switch cnf_p.mission
        case {'CS2','S3'}
            % Based on CS2 and S3 formating
            %{0: open ocean or semi-enclosed seas, 1: enclosed seas or lakes, 2: continental ice, 3: land/transponder}
            %------ open ocean ------------------------------------------
            idx_int=res.surf_type_flag==0;
            %values
            %to use GIM or Bent model
            if any(res.gim_ion(idx_int))
                %from product handbook of Cryosat-2 GIM is nominal choice
                %and if not available use Bent model
                iono_corr=res.gim_ion(idx_int);
                res.flags.gim_ion(idx_int)=1;                
            else
                iono_corr=res.model_ion(idx_int);
                res.flags.model_ion(idx_int)=1;
            end
            res.total_geo_corr(idx_int)=res.ocean_equilibrium_tide(idx_int)+...
                res.ocean_longperiod_tide(idx_int)+...
                res.ocean_loading_tide(idx_int)+...
                res.solidearth_tide(idx_int)+...
                res.geocentric_polar_tide(idx_int)+...
                res.dry_trop(idx_int)+...
                res.wet_trop(idx_int)+...
                iono_corr+...
                res.dac(idx_int); %missing sea state bias
            %flags:
            res.flags.ocean_equilibrium_tide(idx_int)=1;
            res.flags.ocean_longperiod_tide(idx_int)=1;
            res.flags.ocean_loading_tide(idx_int)=1;
            res.flags.solidearth_tide(idx_int)=1;
            res.flags.geocentric_polar_tide(idx_int)=1;
            res.flags.dry_trop(idx_int)=1;
            res.flags.wet_trop(idx_int)=1;
            res.flags.dac(idx_int)=1;
            
            clear idx_int;
            %------ land ice ------------------------------------------
            idx_int=res.surf_type_flag==2;
            %values:
            %to use GIM or Bent model
            if any(res.gim_ion(idx_int))
                %from product handbook of Cryosat-2 GIM is nominal choice
                %and if not available use Bent model
                iono_corr=res.gim_ion(idx_int);
                res.flags.gim_ion(idx_int)=1; 
            else
                iono_corr=res.model_ion(idx_int);
                res.flags.model_ion(idx_int)=1;
            end
            res.total_geo_corr(idx_int)=res.ocean_loading_tide(idx_int)+...
                res.solidearth_tide(idx_int)+...
                res.geocentric_polar_tide(idx_int)+...
                res.dry_trop(idx_int)+...
                res.wet_trop(idx_int)+...
                iono_corr;
            %flags:
            res.flags.ocean_loading_tide(idx_int)=1;
            res.flags.solidearth_tide(idx_int)=1;
            res.flags.geocentric_polar_tide(idx_int)=1;
            res.flags.dry_trop(idx_int)=1;
            res.flags.wet_trop(idx_int)=1;
            
            %------ land ------------------------------------------
            idx_int=res.surf_type_flag==3;
            %values:
            %to use GIM or Bent model
            if any(res.gim_ion(idx_int))
                %from product handbook of Cryosat-2 GIM is nominal choice
                %and if not available use Bent model
                iono_corr=res.gim_ion(idx_int);
                res.flags.gim_ion(idx_int)=1; 
            else
                iono_corr=res.model_ion(idx_int);
                res.flags.model_ion(idx_int)=1;
            end
            res.total_geo_corr(idx_int)=...
                res.solidearth_tide(idx_int)+...
                res.geocentric_polar_tide(idx_int)+...
                res.dry_trop(idx_int)+...
                res.wet_trop(idx_int)+...
                iono_corr;
            %flags:
            res.flags.solidearth_tide(idx_int)=1;
            res.flags.geocentric_polar_tide(idx_int)=1;
            res.flags.dry_trop(idx_int)=1;
            res.flags.wet_trop(idx_int)=1;

        case 'S6'
            %{0: open ocean, 1: land, 2: continental_water, 3: acquatic vegetation, 4: continental_ice_snow, 5: floating_ice, 6: salted_basin}
    end
end


end

