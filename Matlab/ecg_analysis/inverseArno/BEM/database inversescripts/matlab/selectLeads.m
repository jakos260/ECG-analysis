% function GEOM = selectLeads(GEOM,useLeads,useZeromean,wctleads)
%
% GEOM          = GEOM-structure including GEOM.BSM and GEOM.AMA
% useLeads      = leads included in analyses
% useZeromean   = [1] zeromean BSM and A-matrix, [2] apply WCT to A-matrix and BSM, 
%                 [-2] 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured.
% wctleads      = lead numbers in BSM to be used for WCT. [default = 1:3]
%
% Copy of selectLeads.m 15-jan-2015: version 1.0 AJ [added wctLeads]

function GEOM = selectLeads(GEOM,useLeads,useZeromean, wctleads)

if nargin < 4, wctleads    = 1:3; end

GEOM.BSMall = GEOM.BSM;
GEOM.AMAall = GEOM.AMA;

% check if highest number in useLeads is higher than number of BSM leads:
if max(useLeads) > size(GEOM.BSM,1),
    useBSMleads = 1:size(GEOM.BSM,1);
else
    useBSMleads = useLeads;
end

% Correct values with: [1] zeromean for BSM and A-matrix, [2] WCT for BSM and A-matrix, [-2] WCT for A-matrix only!
if useZeromean == 1,
    GEOM.BSM    = zeromean(GEOM.BSM(useBSMleads,:));
    GEOM.AMA    = zeromean(GEOM.AMA(useLeads,:));
elseif useZeromean == -2, % Assumed measured Wilson Central Terminal [rows 1:3], only correction A-matrix    
    GEOM.BSM    = GEOM.BSM(useBSMleads,:);
    GEOM.AMA    = GEOM.AMA(useLeads,:) - ones(length(useLeads),1)*mean(GEOM.AMA(wctleads,:));
    GEOM.AMAall = GEOM.AMAall - ones(size(GEOM.AMAall,1),1)*mean(GEOM.AMA(wctleads,:));
elseif useZeromean == 2, % Wilson Central Terminal   
    GEOM.BSM    = GEOM.BSM(useBSMleads,:) - ones(length(useBSMleads),1)*mean(GEOM.BSM(wctleads,:));
    GEOM.AMA    = GEOM.AMA(useLeads,:) - ones(length(useLeads),1)*mean(GEOM.AMA(wctleads,:));
    GEOM.AMAall = GEOM.AMAall - ones(size(GEOM.AMAall,1),1)*mean(GEOM.AMA(wctleads,:));
else
    GEOM.BSM    = GEOM.BSM(useBSMleads,:);
    GEOM.AMA    = GEOM.AMA(useLeads,:);
end

if strcmp(GEOM.type,'ventricles'), else end