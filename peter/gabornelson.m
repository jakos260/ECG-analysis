function vec=gabornelson(normals,area,ECG)

for i=1:size(normals,1)
    normals(i,:) = normals(i,:) ./ norm3d( normals(i,:) );
end



ampl = (ECG .* (area * ones(1,size(ECG,2))))';

vec = [ampl * normals(:,1) ampl * normals(:,2) ampl * normals(:,3)]./length(normals);

