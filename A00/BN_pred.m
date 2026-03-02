function [pmf,bin] = BN_pred(x,b)
%% returns binomial pmf with incorrect response being the 1st bin
xPred = log(1+ exp(x - b(:,1))).* b(:,2);
prCorr = 1./(1+exp(-(-log(5) + xPred)));
prInco = 1 - prCorr;
pmf = [prInco,prCorr];
bin = [0,1];
return