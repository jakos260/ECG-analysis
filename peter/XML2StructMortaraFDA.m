function varargout = XML2StructMortaraFDA(filename) 

A=xml2struct(filename);

C = A.AnnotatedECG.component.series.component.sequenceSet.component;


sampletime = str2num(C{1}.sequence.value.increment.Attributes.value);

ECG=[];
names=[];
for i=2:length(C)
    names{i-1} = C{i}.sequence.code.Attributes.code;
        
    
    value = str2num(C{i}.sequence.value.scale.Attributes.value);
    unit = C{i}.sequence.value.scale.Attributes.unit;
    factor = 1;
    if strcmp(unit,'uV')
        factor = 1000;
    end
    ecg = str2num(C{i}.sequence.value.digits.Text) * value/factor;

    ECG=[ECG; ecg];
end
varargout{1} = ECG;

if nargout >1
    varargout{2} = names;
end

