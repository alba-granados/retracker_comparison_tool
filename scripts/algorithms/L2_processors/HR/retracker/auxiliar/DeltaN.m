function [DeltaN] = DeltaN(n)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% This code generates the configuration file input to isardSATs
% re-tracker as described in isardSAT_Retracker_DPM_V1.a.docx internal
% document
% 
%
% ---------------------------------------------------------
% DeltaN: This code measures the sinc as it comes from the ATFFT
%
% Calling
%   DeltaN(n)
%
% Inputs
%   n - rank of values it determines number of samples
%
%
% Output
%   DeltaN fon to be fitted
%
%
% ----------------------------------------------------------
% 
% Author:   Cristina Martin-Puig / isardSAT
%
% Reviewer: Cristina Martin-Puig / isardSAT
%
% Last revision: Cristina Martin-Puig / isardSAT 25/02/2013
%
%
% This software is subject to the conditions 
% set forth in the ESA contract P4 Jason-CS GPP
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DeltaN = exp(1i*pi*(length(n)-1)*n/length(n)).*sinc(n);

end
