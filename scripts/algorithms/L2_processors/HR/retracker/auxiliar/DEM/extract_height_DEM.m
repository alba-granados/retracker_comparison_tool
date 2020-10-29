function [height_DEM_ell,height_DEM_geo] = extract_height_DEM (LAT,LON,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for extraction of the height of surfaces beneath
% satellite from a DEM (measured w.r.t reference ellipsoid)
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
%       LAT       =   latitudes of the points of interest
%       LON       =   longitude of the points of interest
% OUTPUT:
%       height_DEM_ell        =   surface height beneath the satellite from DEM
%       and referred to the reference elliposid.
%       height_DEM_geo        =   surface height beneath the satellite from DEM
%       and referred to the geoid.
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
% Assuming the use of the SRTM DEM
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Simply based on a interpolation of the height for the specified
% longitude and latitude from a given DEM

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(2+2*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('DEM_ref','SRTM',@(x)ischar(x));
p.addParamValue('DEM_dir','.',@(x)ischar(x));
p.parse(varargin{:});
DEM_ref=p.Results.DEM_ref;
DEM_dir=p.Results.DEM_dir;
clear p;

%% ----------------- Load DEM ---------------------------------------------
switch DEM_ref
    case 'SRTM'
        %format DEM.lat DEM.lon and DEM.z (height)
        min_LAT=min(LAT);
        max_LAT=max(LAT);
        min_LON=min(LON);
        max_LON=max(LON);
        if max_LAT-min_LAT<0.002 %tolerances
           max_LAT=min_LAT+0.002;
        end
        if max_LON-min_LON<0.002 %tolerances
            max_LON=min_LON+0.002;
        end
        DEM=readhgt([min_LAT,max_LAT,min_LON,max_LON],'outdir',DEM_dir);
        DEM.z((DEM.z==-32768))=NaN;  
        N=geoidheight(LAT,LON,'EGM96');% Geoid height EGM96 w.r.t reference ellipsoid WGS84
    otherwise
        error('Not a valid type of DEM');
end

%% ---------------- Interpolate the heights -------------------------------
s=size(DEM.z);
lon_grid=repmat(DEM.lon,[s(1),1]);
lat_grid=repmat(DEM.lat,[1,s(2)]);
F = scatteredInterpolant(lon_grid(:),lat_grid(:),double(DEM.z(:))) ;
height_DEM_geo = F(LON,LAT);
%extract the window delay referencing the height w.r.t ellipsoid 
%(include the geoid height w.r.t ellipsoid: N)
height_DEM_ell=height_DEM_geo+N;


end

