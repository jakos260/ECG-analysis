
% cases = dir('proCips*');
% cd ('C:\Users\peter\CIPS\proCIPS\proCips\results\');cases = dir('proCips*');
cd('C:\Users\peter\CIPS\medtronic\results_coronly\');cases = dir('pa*');
ECG=[];
for i=1:length(cases)
    if cases(i).isdir
        a = dir(cases(i).name);
%         if isempty(strfind(cases(i).name,'04'))
%             continue;
%         end
        k=1;
        for i1=1:length(a)
            if a(i1).isdir && isempty(strfind(a(i1).name,'.'))
                name = fullfile(cases(i).name,a(i1).name);
                patname  = cases(i).name;
                b=dir(name);
                for i2=1:length(b)
                    if b(i2).isdir && isempty(strfind(b(i2).name,'.'))
                        name2 = fullfile(cases(i).name,a(i1).name,b(i2).name);
                        c=dir(fullfile(name2,'*.ecg'));
                        for j=1:length(c)
                            ECG = loadmat(fullfile(name2,c(j).name));
                            ECGS{i}.ecgs{k}=ECG;
                            ECGS{i}.names{k}=fullfile(name2,c(j).name);
                            k=k+1;
                            figure(1);leadv12(ECG(:,1:min(4000,size(ECG,2))),'paperspeed',25,'max',3);
                            saveas(figure(1),fullfile(name2,[c(j).name '.png']))
                        end                        
                    end
                end                
            end
        end        
    end
end
ecgs=ECGS{end}.ecgs;

vcgs=[];
for i=1:length(ecgs)
    
%     if ~isempty(strfind(patname ,'15')) && i==3
%         ecgs{i} = ecgs{i}(:,605:end-403);        
%     else        
%         ecgs{i} = ecgs{i}(:,600:end-398);
%     end
    vcgs{i} = ecgs{i}(4:end,380:600);
end

if ~isempty(strfind(patname ,'08'))
    ecgs(2)=[]; % probably a pwave
elseif ~isempty(strfind(patname ,'09'))
    ecgs(4:end)=[]; 
elseif ~isempty(strfind(patname ,'15'))
    ecgs(2,4:9)=[]; 
end

        

figure(2);leadv12(ecgs(),'paperspeed',50)
% figure(2);leadv12(ecgs([1 3 10 11]),'paperspeed',50,'max',1.5)

saveas(figure(2),fullfile(name,[patname '.png']))




