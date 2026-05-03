function vec=gabornelson(normals,area,ECG)

ampl = (ECG .* (area * ones(1,size(ECG,2))))';

vec = [ampl * normals(:,1) ampl * normals(:,2) ampl * normals(:,3)]./length(normals);

