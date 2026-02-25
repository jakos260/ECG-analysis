function BSM = standard12lead2BSM_cipsModel(modeldir, modelname, ECG12lead,varargin)


DATA = readGeomPeacsModel(modeldir,modelname);


elecs = DATA.VENTR.LEADPOS.A_standard12lead;

elecsI=zeros(length(elecs),1);
for i=1:length(elecsI)
    d= [DATA.GEOM.thorax.VER(:,1) - elecs(i,1)...
        DATA.GEOM.thorax.VER(:,2) - elecs(i,2)...
        DATA.GEOM.thorax.VER(:,3) - elecs(i,3)];
    d = norm3d(d);
    elecsI(i) = find(d==min(d));
end

disp('assumed 12 lead ECG used aVr, aVl , and aVf')
ECG12lead(4:6,:) = ECG12lead(4:6,:) / 1.5;

if ~isempty(varargin)
    wct_adapt = varargin{1};
    if size(wct_adapt,1)==1
        wct_adapt = wct_adapt';
    end
    elecsI(1:length(wct_adapt)) = wct_adapt;
end

   
T = intripol(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,elecsI);
BSM= T * ECG12lead(4:end,:); % the first leads are I, II , III


