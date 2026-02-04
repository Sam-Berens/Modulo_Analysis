function [epsilon] = epsilonFun(M)
%% epsilon is a 1 x 2 for coloc = -1 and +1
% Set mean similarity selectors
persistent S;
if isempty(S) || any(structfun(@isempty,S))
    [A,B] = meshgrid(0:5,0:5);
    Db = min(mod(A-B,6),mod(B-A,6));

    D = kron([0,0;1,0],Db);
    D(triu(true(12))) = NaN;
    S.n1 = D==1;
    S.n2 = D==2;
    S.n3 = D==3;

    D = kron(eye(2),Db);
    D(triu(true(12))) = NaN;
    S.p1 = D==1;
    S.p2 = D==2;
    S.p3 = D==3;
end

% Preallocate the epsilon for both conditions
epsilon = nan(1,2);

% Compute the neural similarity
R = corr(M);

% Extract mean similarities (colocation=-1)
mu1 = mean(R(S.n1));
mu2 = mean(R(S.n2));
mu3 = mean(R(S.n3));

% Return if the target inequality does not hold
if ~(isnan(mu1) || (mu1 < mu2) || (mu2 < mu3))
    rho = [mu1;mu2;mu3];
    epsilon(1) = getEpsilon(rho);
end

% Extract mean similarities (colocation=+1)
mu1 = mean(R(S.p1));
mu2 = mean(R(S.p2));
mu3 = mean(R(S.p3));

% Return if the target inequality does not hold
if ~(isnan(mu1) || (mu1 < mu2) || (mu2 < mu3))
    rho = [mu1;mu2;mu3];
    epsilon(2) = getEpsilon(rho);
end

return


function prop = getEpsilon(rho)
prop = (rho(2) - rho(1))/(rho(3)-rho(1)); % coarse (0) < 0.5 < fine (1)
prop = prop - 0.5; % (-0.5) coarse < 0 < fine (0.5)
prop = prop * 2;  %  (-1) coarse < 0.5 < fine (1)
return
