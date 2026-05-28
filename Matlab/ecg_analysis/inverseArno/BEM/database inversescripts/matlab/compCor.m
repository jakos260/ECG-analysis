function cor = compCor(MEAS,SIM,mode)


if size(MEAS) ~= size(SIM)
    error('PHI and PSI not the same size');
end
if mode > 10
    meas=bsxfun(@rdivide,MEAS,rms(MEAS) + 0.01);
    sim=bsxfun(@rdivide,SIM,rms(MEAS) + 0.01);
end
if mode == 1
    COR = corrcoef(MEAS,SIM);
    cor = COR(2,1);
%     cor = mean(MEAS.*SIM)./sqrt(std(MEAS).*std(SIM));    
%     cor(isnan(cor))=[];
%     cor = mean(cor);
%     cor = sum(MEAS(:).*SIM(:))/sqrt(sum(MEAS(:).^2) * sum(SIM(:).^2))
elseif mode == 11
    cor = s(MEAS(:).*SIM(:))/(std(MEAS(:)).*std(SIM(:)));
elseif mode == 2
    cor = norm(MEAS-SIM,'fro')/norm(MEAS,'fro');
elseif mode == 22
    cor = norm(meas-sim,'fro')/norm(meas,'fro');
elseif mode == 3
    cor = MEAS./SIM;
    cor(isinf(cor)) = 0;
    cor(SIM<0.0001) = 0;
    cor = 1-sum(cor(:))/size(MEAS,2);
elseif mode == 4
    b = abs(MEAS);
    maxMeas = max(b,[],2);
    limit = maxMeas/1000;
    c=bsxfun(@max,b,limit);
    c= max(b,0.001);
    cor = (MEAS.*SIM)./c;
    cor = 1- sum(abs(cor(:)))/size(MEAS,2);
elseif mode == 5
%     drms=rms(MEAS-SIM);
%     rrms=rms(MEAS);
%     rrms(rrms < max(rrms)/100) = max(rrms)/100;
%     cor = 1 - rms((drms./rrms)');
    cor = (((MEAS-SIM).^2)./(abs(MEAS)+0.01^2));
    cor = ((MEAS-SIM).^2);
    cor=bsxfun(@rdivide,cor,rms(MEAS)+0.01);
    cor = 1-mean((cor(:))');%sqrt(mean(cor(:)));
%     cor = mean(cor(abs(cor)< 10));
end
    
    
        

