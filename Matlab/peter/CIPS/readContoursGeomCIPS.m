function contours = readContoursGeomCIPS(fn)
contours=[];
XML=XML2StructPeter(fn);
C=XML.Contours;
fl=fields(C);
for i=1:length(fl)
    flp = getfield(C,fl{i});
    if ~isempty( flp.Points )
        points = reshape(flp.Points,3,length(flp.Points) /3)';
        points = [-points(:,2) points(:,1) points(:,3)];
        contours. (fl{i}) = points;
    end
end