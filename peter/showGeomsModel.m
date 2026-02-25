function showGeomsModel(dirname,showThorax,crossval)

qtriplot('delete *');
pause (0.5);
dirname = fullfile(dirname);
a=strfind(dirname,filesep);
if a(end)== length(dirname)
    a(end)=[];
end
basename = [dirname filesep dirname(a(end)+1:end)];

[VVER,VITRI]= loadtri([basename '_ventricles.tri']);
[AVER,AITRI]= loadtri([basename '_atria.tri']);
[TVER,TITRI]= loadtri([basename '_thorax.tri']);
[LLVER,LLITRI]= loadtri([basename '_llung.tri']);
[RLVER,RLITRI]= loadtri([basename '_rlung.tri']);
[LVER,LITRI]= loadtri([basename '_lcav.tri']);
[RVER,RITRI]= loadtri([basename '_rcav.tri']);

[LADVER,LADITRI]= loadtri([basename '_lad.tri']);
[RCAVER,RCAITRI]= loadtri([basename '_rca.tri']);
[LCXVER,LCXITRI]= loadtri([basename '_lcx.tri']);

crosstxt = [ 'cross ' num2str(crossval) ' '];

liver=0;
if exist([basename '_liver.tri'],'file')
    [LIVER,LIITRI]= loadtri([basename '_liver.tri']);
    liver = 1;
end


qtriplot('bgdcolor white');
pause (0.5);

qtriplot(LLVER,LLITRI);
qtriplot('color darkgrey')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')

qtriplot(RLVER,RLITRI);
qtriplot('color darkgrey')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')

qtriplot(VVER,VITRI);
qtriplot('color heart')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')

qtriplot(AVER,AITRI);
qtriplot('color lightred')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')
pause (0.5);

qtriplot(RVER,RITRI);
qtriplot('color blue')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')

qtriplot(LVER,LITRI);
qtriplot('color red')
% qtriplot('cross fill')
% qtriplot(crosstxt)
% qtriplot('cross fill')
pause (0.5);

qtriplot(LADVER,LADITRI);
qtriplot('color red')
% qtriplot('cross fill')

qtriplot(RCAVER,RCAITRI);
qtriplot('color red')
% qtriplot('cross fill')

qtriplot(LCXVER,LCXITRI);
qtriplot('color red')

pause(0.5)
qtriplot('cross fill')
qtriplot(crosstxt)
pause(0.5)
if liver && diff(range(LIVER(:,1))) > 20
qtriplot(LIVER,LIITRI);
qtriplot(['color ' num2str([70/255 0 70/255])])
% qtriplot('cross fill')
qtriplot(crosstxt)
end

if showThorax
qtriplot(TVER,TITRI);
pause(0.5)
qtriplot('color flesh')
pause(0.5)
qtriplot('trans 0.6');
end

qtriplot(['png ' basename '.png']);