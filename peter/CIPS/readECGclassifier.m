function ECGClass = readECGclassifier(fn)

XML=xml2struct(fn);
ECGfiles = XML.ECGClassifier{2}.ECGClassifierDATA.ECG_cases.ECG_CASE.ECGfiles.ECGFiledata;
for i = 1:length(ECGfiles)
    
    ECGClass{i}.filename            = ECGfiles{i}.ECGfilename.Text;
    ECGClass{i}.ECGclassification   = ECGfiles{i}.ECGclassification.Text;
    ECGClass{i}.atrialrhythm        = ECGfiles{i}.atrialrhythm.Text;
    ECGClass{i}.atrialwaveform      = ECGfiles{i}.atrialwaveform.Text;    
    ECGClass{i}.ventricularwaveform = ECGfiles{i}.ventricularwaveform.Text;
    ECGClass{i}.ventricularrhythm   = ECGfiles{i}.ventricularrhythm.Text;
    ECGClass{i}.avnode              = ECGfiles{i}.avnode.Text;
    ECGClass{i}.ACS                 = ECGfiles{i}.ACS.Text;
    ECGClass{i}.lvh                 = ECGfiles{i}.lvh.Text;
end
