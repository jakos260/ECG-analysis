function varargout = findPeakQRS(rmsSignal)

maxS= max(rmsSignal);
maxSel=find( rmsSignal > maxS*0.8 );
segments = maxSel;
segments(diff(segments )==1)=[];
segments(segments < 30) = [];
if isempty(segments) || segments(1) > size(rmsSignal,2)/2
    tpeak=-1;
elseif length( segments ) == 1 % one peak
    tpeak = find(rmsSignal== maxS);
else
    if segments(2) - segments(1) > 100 % large T wave
        tpeak = find(rmsSignal == max(rmsSignal(maxSel(1):segments(1))));
        segments(2:end)=[];
    else
        tpeak = find(rmsSignal == max(rmsSignal(maxSel(1):segments(1))));
        segments(3:end)=[];
    end
end
varargout{1} = tpeak;
if nargout > 1
    if tpeak > 0
        for t=tpeak-30:-1:max(1,tpeak-200)
            if rmsSignal(t) > rmsSignal(t+1) && rmsSignal(t) < rmsSignal(t+10)
                tQRSOnset = t+1;
                break
            end
        end
        for t=tpeak+30:1:tpeak+120
            if rmsSignal(t) > rmsSignal(t-1) && rmsSignal(t) < rmsSignal(t-10)
                tQRSend = t-1;
                break
            end
        end
        estimatedTend = round(size(rmsSignal,2)*0.8);
        if exist('tQRSend','var') && tQRSend < estimatedTend
            tTpeak = find(rmsSignal==max(rmsSignal(tQRSend:estimatedTend)));
            tTend = find(rmsSignal==min(rmsSignal(tTpeak+30:1:tTpeak+120)));
%             for t=tTpeak+30:1:tTpeak+120
%                 if rmsSignal(t) > rmsSignal(t-1) && abs(rmsSignal(t) - rmsSignal(t+5)) < 0.2 * rmsSignal(t)
%                     tTend = t-1;
%                     break
%                 end
%             end
        end
    end
    if exist('tQRSOnset','var') && exist('tQRSend','var') && exist('tTpeak','var') && exist('tTend','var') 
        varargout{2} = [tQRSOnset tQRSend tTpeak tTend];
    else
        varargout{2} = ones(1,4) * -1;
    end
%      figure(200);clf
%      plot(lowpassma(rmsSignal,10))
end

