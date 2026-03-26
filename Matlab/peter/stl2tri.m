dirIn = dir('*.stl');
for i=1:length(dirIn)
    if isempty(strfind(dirIn(i).name,'HQ'))
        [VER,ITRI]=loadSTL([dirIn(i).name ]);
        if size(VER,1) < 6000    
            savetri([dirIn(i).name(1:end-3) 'tri'],VER,ITRI(:,[2 1 3]));
        end
    end
end