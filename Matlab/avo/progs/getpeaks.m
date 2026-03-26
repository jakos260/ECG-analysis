% getpeaks.m
% function PEAKS=getpeaks(signal,threshold,window)
% Find value PEAKS(:,1) and location PEAKS(:,2) of local peaks in a (scalar) signal.
% dedicated to finding RR-intervals in a record containing several beats.
% If window>0 PEAKS(:,2) documents as events the locations of the extremes above the
% threshold
% If window<0 PEAKS(:,2) takes as events the locations of the upsloping threshold
% crossings,in which case PEAKS(:,1)=threshold
% For all such events within an interval window following any event identified are discarded
% For identifying the downslope crossings:
% apply the function to reverse(signal), followed by reverse(PEAKS);
% modified 2013_01_11; A. van Oosterom

function PEAKS=getpeaks(signal,threshold,window)

nt=size(signal,2);
wind=abs(window);
signal(signal<threshold)=0;
PEAKS=[];
% treat possible extremes at the bounds of the interval
augment=[ 0 signal 0 ];
SEARCH=[(1:nt+1)' augment(1:nt+1)' augment(2:nt+2)'];
starters=SEARCH(SEARCH(:,2)==0 & SEARCH(:,3)~=0);
endings=SEARCH(SEARCH(:,2)~=0 & SEARCH(:,3)==0)-1;
npeaks=size(starters,1);

if window>0,
    if npeaks>0,
        for i=1:npeaks,
            % endspan=min(starters(i)+window-1,nt);
            % [PEAKS(i,1) PEAKS(i,2)]=max(signal(starters(i):endspan));
            [PEAKS(i,1) PEAKS(i,2)]=max(signal(starters(i):endings(i)));
        end
        PEAKS(:,2)=PEAKS(:,2)+ starters-1;
    end
else
    signal(starters);
    PEAKS=[threshold*ones(size(starters,1),1) starters];
end

% pause

SEARCH=[];
while ~isempty(PEAKS),
    SEARCH=[SEARCH;PEAKS(1,:)];
    PEAKS(PEAKS(:,2)<PEAKS(1,2)+wind,:)=[];
end
PEAKS=SEARCH;




