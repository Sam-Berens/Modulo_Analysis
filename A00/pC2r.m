function [r] = pC2r(pC)
% Returns the resultant vectors (r) corresponding to the model-predicted
% probability of a correct response on the first attempt of a trial in the
% hspvm model. This function is used to compute "memorisation points" -
% i.e., when during training, a given pair has been well-learnt.

% Th function f() takes a value of kappa (the von Mises concentration
% parameter and returns the probability of a correct response.
f = @(k)exp(k)/sum(exp(k*cos((0:5)*pi/3)));

% The following call to fmincon inverts f() for the values in pC. That is,
% it identifies the values of kappa that map to values of pC.
kappa = nan(size(pC));
for ii = 1:numel(pC)
    kappa(ii) = fminunc(@(k)(f(k)-pC(ii)).^2,1,...
        optimoptions(@fminunc,'Display','off'));
end

% Finally, we solve one of four polynomials to invert the r -> kappa
% approximation provided by Fisher (1993). This returns a the values r that
% correspond to the the values of kappa computed above.
[r]= k2r(kappa);
return

function [r] = k2r(kappa)
kappaMag = abs(kappa);
kappaSgn = sign(kappa);
r = nan(size(kappa));
for ii = 1:numel(kappa)
    k = kappaMag(ii);
    if k < 7e-3
        r(ii) = k;
    elseif (k < 3.648) && (k >= 1.2516)
        r(ii) = (-sqrt((10000*k^2 - 19800*k + 33709)) +100*k + 179) / 278;
    elseif k < 1.2516
        temp = roots([5,0,6,0,12,(-6*k)]);
        r(ii) = abs(temp(end));
    else
        temp = roots([1,-4,3,(-1/k)]);
        r(ii) = abs(temp(2));
    end
end
r(r>1) = 1;
r = r.*kappaSgn;
return