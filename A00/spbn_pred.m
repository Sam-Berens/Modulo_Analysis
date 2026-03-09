function [Pmf,bin] = spbn_pred(P,x)
% [Pmf,bin] = spbn_pred(x,b)
%    Predicts the probability mass function (Pmf) using the Softplus
%    binomial model.
x = x(:)';
a = P(:,1);
b = P(:,2);
xPred = log(1 + exp(x-a)) .* b;
prCorr = 1./(1+exp(log(5)-xPred));
prInco = 1 - prCorr;
Pmf = [prInco,prCorr];
bin = [0,1];
return