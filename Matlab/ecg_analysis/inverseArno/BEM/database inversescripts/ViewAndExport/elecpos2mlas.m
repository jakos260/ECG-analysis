function [mla, standard,leadsys]=elecpos2mlas(SelectChannels,ElecPos)
chanswitch = []; % only use when needed
if ~exist('ElecPos'), ElecPos = 'pigbsm'; end

if ischar(ElecPos),
    switch lower(ElecPos),
        case {'amsbdf'}
            vl=66;
            vr=67;
            vf=65;%Use (bdf) channel numbers
            %vl66,vr67, vf65 should end up as 63, 64, 65
            %amsterdam
            ElecPos = [
                vr,00,00,00,00,00,00,vl,00,00,00,00,00,00,00,00;
                01,00,04,10,00,00,00,36,42,00,00,54,00,00,00,00;
                00,00,00,00,16,00,28,00,43,00,00,00,00,57,00,00;
                00,00,05,00,17,00,29,37,00,00,51,00,56,00,00,60;
                00,00,00,11,00,24,00,38,00,00,00,00,00,00,00,00;
                02,00,06,12,18,00,30,39,44,00,52,00,00,00,00,00;
                00,00,00,13,19,25,31,40,45,00,00,00,00,58,00,61;
                00,00,07,14,20,26,32,00,46,49,00,00,00,00,00,00;
                03,00,00,15,21,00,33,41,00,00,00,55,00,00,00,00;
                00,00,08,00,22,27,00,00,47,50,00,00,00,00,00,00;
                00,00,00,00,00,00,34,00,00,00,00,00,00,59,00,00;
                00,00,09,00,23,00,00,00,48,00,53,00,00,00,00,62;
                00,00,00,00,00,00,35,00,00,00,00,00,00,00,00,00;
                00,00,00,00,00,00,vf,00,00,00,00,00,00,00,00,00
                ];
            standard=[12;18;25;31;40;45;vr;vl;vf];
            if ~all(ismember(standard,SelectChannels))
                standard=[];
            end
            leadsys='ams';
            %             chanswitch=[1:62,66,67,65]; % BSM, VL, VR, VF
            %VL63=67, vr 64=65, vf 65=66, , maar is hier nog mee
            %geschoven?
        case {'amszip'}
            vl=63;
            vr=64;
            vf=65;%Use (bdf) channel numbers
            %vl66,vr67, vf65 should end up as 63, 64, 65
            %amsterdam
            ElecPos = [
                vr,00,00,00,00,00,00,vl,00,00,00,00,00,00,00,00;
                01,00,04,10,00,00,00,36,42,00,00,54,00,00,00,00;
                00,00,00,00,16,00,28,00,43,00,00,00,00,57,00,00;
                00,00,05,00,17,00,29,37,00,00,51,00,56,00,00,60;
                00,00,00,11,00,24,00,38,00,00,00,00,00,00,00,00;
                02,00,06,12,18,00,30,39,44,00,52,00,00,00,00,00;
                00,00,00,13,19,25,31,40,45,00,00,00,00,58,00,61;
                00,00,07,14,20,26,32,00,46,49,00,00,00,00,00,00;
                03,00,00,15,21,00,33,41,00,00,00,55,00,00,00,00;
                00,00,08,00,22,27,00,00,47,50,00,00,00,00,00,00;
                00,00,00,00,00,00,34,00,00,00,00,00,00,59,00,00;
                00,00,09,00,23,00,00,00,48,00,53,00,00,00,00,62;
                00,00,00,00,00,00,35,00,00,00,00,00,00,00,00,00;
                00,00,00,00,00,00,vf,00,00,00,00,00,00,00,00,00
                ];
            standard=[12;18;25;31;40;45;vr;vl;vf];
            
            if ~all(ismember(standard,SelectChannels))
                standard=[];
            end
            leadsys='ams';
            %             chanswitch=[1:62,66,67,65]; % BSM, VL, VR, VF
            %VL63=67, vr 64=65, vf 65=66, , maar is hier nog mee
            %geschoven?
        case {'pigbsm','pbsm','pigbsmegm','pigbsmegm09'}
            %pigbsm
            chanarr=1:60;
            mla(:,2)=mod(chanarr-1,6)+2;%+1 for the mod trick, only limb leads on 1
            mla(:,1)=floor((chanarr-1)/6)+1;
            mla(:,3)=chanarr;
            mla(31:60,3)=2+mla(31:60,3); % 31 and 32 used as limb leads
            mla=[mla;2,1,31;2,8,32;9,1,63;9,8,64];
            if strcmp(lower(ElecPos),'pigbsmegm')
                mla=[mla;4,8,65;5,8,66;6,8,67;7,8,68]; % add 4 egm channels on lower row
            end
            if strcmp(lower(ElecPos),'pigbsmegm09')
                mla=[mla;4,8,76;5,8,75;6,8,65;7,8,66]; % add 4 egm channels on lower row, A, Endo, Epidist, Epi2
            end
            mla=sortrows(mla,3);
            standard=[];
            leadsys='pigbsm';
            %             return % oostep1: cut from ViewandExport. Too lazy to rewrite this
        case '9lead'
%             mla=loadmat('9lds.mla');
%             mla=mla(2:end,:);
%             mla(1,:)=[9,1,0];
            mla(1:9,:)=[ones(1,9);1:9;1:9]';
            standard=[];
            leadsys='9leads';
        otherwise
            error('Unknown electrode configuration');
    end
end
if ~exist('mla')
    mla(1,1:2)=size(ElecPos);
    mla(1,3)=0;
    for i=1:size(ElecPos,1)
        for j=1:size(ElecPos,2)
            if ElecPos(i,j)>0
                mla(end+1,:)=[j,i,ElecPos(i,j)];
            end
        end
    end
    leadsys='';
end