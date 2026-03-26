function [ypred,predict,RESULT, thresholds,paramGrad] = logisticregression( x, target )

ytrain = target; % Target variable
xtrain = x;
doNormalize = 2;
if doNormalize == 1
% xtrain = bsxfun(@rdivide,bsxfun(@minus,x,mean(x)), std(x) * 2);% Normalized Predictors
    xtrain = bsxfun(@minus,x,mean(x));% shift to mean Normalized Predictors
    xtrain = bsxfun(@rdivide,xtrain, std(xtrain));% Normalized Predictors
elseif doNormalize == 2
    xtrain = bsxfun(@minus,x,median(x));% shift to mean Normalized Predictors
    XL= quartile(xtrain,0.1);
    XH= quartile(xtrain,0.9);
    xtrain = bsxfun(@rdivide,xtrain, XH-XL);% Normalized Predictors
end
xtrain=[ones(length(xtrain),1) xtrain]; % one is added for calculation of biases.
xtest=xtrain;
ytest=ytrain;

%compute cost and gradient

theta=zeros(size(xtrain,2),1); % Initial weights


[J, grad, thresholds] = costFunction(theta,xtrain,ytrain); % Cost funtion

paramGrad= sortrows([(1:size(x,2))' grad(2:end) abs(grad(2:end)./min(abs(grad(2:end))))],3);

ypred = xtest * thresholds; %target prediction

% probability calculation
predict=sigmoid(ypred); % Hypothesis Function

ypred(predict >= 0.5)=1;
ypred(predict <  0.5)=0;

TP = sum(ypred == 1 & target == 1);
TN = sum(ypred == 0 & target == 0);
FP = sum(ypred == 1 & target == 0);
FN = sum(ypred == 0 & target == 1);
% make confusion matrix here

RESULT.sensitivity = TP / (TP + FN);
RESULT.specificity = TN / (TN + FP);
RESULT.PPV         = TP / (TP + FP);
RESULT.NPV         = TN / (TN + FN);



% Decision Boundary
%%
figure;clf

subplot(1,15,1:7)
hold on
scatter(xtest(ytest==1,2),xtest(ytest==1,3),'b+','linewidth',2.0)
scatter(xtest(ytest==0,2),xtest(ytest==0,3),'r.','linewidth',2.0)
legend('Pos class','Neg. class')
subplot(1,15,8:9); hold on
trueNeg=find(ypred == 0 & target ==0);
falseNeg=find(ypred == 1 & target ==0);
truePos=find(ypred == 1 & target ==1);
falsePos=find(ypred == 0 & target ==1);

scatter(target(truePos),predict(truePos),'+g')
scatter(target(falsePos),predict(falsePos),'xr')
scatter(target(trueNeg),predict(trueNeg),'og')
scatter(target(falseNeg),predict(falseNeg),'or')
axis([-0.2 1.2 0 1])
subplot(1,15,10:15)
plot(J(3:end))


% J

%%
function [J, grad, th] = costFunction(theta, xtrain,ytrain)
% Cost Summary of this function goes here

maxiter     = 1000; % No. of iterations for weight updation
alphaStart  = 0.1; % Learning parameter
alpha       = alphaStart / size(xtrain,2);
th = theta;
m  = size(xtrain,1);
J  = ones(maxiter,1);



for j=1:maxiter
    h = sigmoid( xtrain * th );
    J(j) = -(1 / m) * sum(ytrain .* log(h) + (1-ytrain) .* log(1-h));
    th = th + alpha * xtrain' * (ytrain - h);    
    if j > 1         
        if J(j-1) - J(j) < 1e-5
            break;
        elseif J(j-1) - J(j) < 5e-3
            alpha = max(0.01 / size(xtrain,2),alpha / 2);
        end
    end        
end
J=J(1:j);
grad = zeros(size(theta,1),1);
for i = 1:size(grad)
    grad(i) = (1/m) * sum(( h - ytrain )' * xtrain(:,i));
end

%%
function g = sigmoid( z )
g=1./(1+exp(-z));

%%
function values = quartile(dataIn,quart)


if size(dataIn,1)==1
    data = dataIn';
else
    data = dataIn;
end
values = zeros(1,size(data,2));
for i=1:size(data,2)
    ds = sort(data(:,i));
    i1= floor(quart*length(ds));
    i2= ceil(quart*length(ds));
    di = quart*length(ds) - i1;
    values(i) = ds(i1) - (ds(i2)- ds(i1)) *di;    
%     values(i) = (ds(i1)  + ds(i2))/2;    
end

