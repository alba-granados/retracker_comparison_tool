function [data] = extract_ref_height (data, cnf_p,cst_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for extraction of the height & window_delay of surfaces beneath satellite from a DEM
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 13/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       data        =   structure of data as defined by our L2 processor
%       cnf_p       =   configuration parameters
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%       
%  
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Simply based on a interpolation of the height for the specified
% longitude and latitude from a given DEM

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(3+4*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('DEM_ref','SRTM',@(x)ischar(x));
p.addParamValue('dir_DEM','.',@(x)ischar(x));
p.addParamValue('path_Results',{''},@(x)ischar(x));
p.addParamValue('L1B_filename',{''},@(x)ischar(x));
p.parse(varargin{:});
DEM_ref=p.Results.DEM_ref;
dir_DEM=p.Results.dir_DEM;
path_Results=char(p.Results.path_Results);
L1B_filename=char(p.Results.L1B_filename);
clear p;

%% ----------------- Load DEM ---------------------------------------------
switch DEM_ref
    case 'SRTM'
        %format DEM.lat DEM.lon and DEM.z (height)
        min_LAT=min(data.GEO.LAT);
        max_LAT=max(data.GEO.LAT);
        min_LON=min(data.GEO.LON);
        max_LON=max(data.GEO.LON);
        if max_LAT-min_LAT<0.002 %tolerances
            max_LAT=min_LAT+0.002;
        end
        if max_LON-min_LON<0.002 %tolerances
            max_LON=min_LON+0.002;
        end
        DEM=readhgt([min_LAT,max_LAT,min_LON,max_LON],'outdir',dir_DEM);
        %DEM=readhgt([min(data.GEO.LAT),max(data.GEO.LAT),min(data.GEO.LON),max(data.GEO.LON)],'outdir',dir_DEM);
        DEM.z((DEM.z==-32768))=NaN;  
        data.GEO.N=geoidheight(data.GEO.LAT,data.GEO.LON,'EGM96');% Geoid height EGM96 w.r.t reference ellipsoid WGS84
    otherwise
        error('Not a valid type of DEM');
end

%% ---------------- Interpolate the heights -------------------------------
s=size(DEM.z);
lon_grid=repmat(DEM.lon,[s(1),1]);
lat_grid=repmat(DEM.lat,[1,s(2)]);
F = scatteredInterpolant(lon_grid(:),lat_grid(:),double(DEM.z(:))) ;
data.GEO.DEM_h = F(data.GEO.LON,data.GEO.LAT);
%extract the window delay referencing the height w.r.t ellipsoid 
%(include the geoid height w.r.t ellipsoid: N)
data.GEO.wd_ref_DEM=2.0*(data.GEO.H-(data.GEO.DEM_h+data.GEO.N))./cst_p.c_cst;

f1=figure;
plot(data.GEO.LAT,data.GEO.H-data.MEA.win_delay*cst_p.c_cst/2,'-b'); 
hold on; plot(data.GEO.LAT,data.GEO.H-data.GEO.wd_ref_DEM*cst_p.c_cst/2,'-r'); 
xlabel('Latitude [deg.]'); ylabel('Height [m]'); 
legend('SH: H_{orb}-R (Altimeter)','SH: DEM SRTM + geoid height')
print('-dpng ',[path_Results,'plots',filesep,L1B_filename(17:47),'_Height_surf_vs_height_DEM.png']);
close(f1)




end

