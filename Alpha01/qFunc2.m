function [q] = qFunc2(M)

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
    % block masks
    lower6 = tril(true(6), -1);
    maskA = [lower6, false(6); false(6), false(6)];
    maskB = [false(6), false(6); false(6), lower6];

    S.p1.a = D==1 & maskA;
    S.p2.a = D==2 & maskA;
    S.p3.a = D==3 & maskA;

    S.p1.b = D==1 & maskB;
    S.p2.b = D==2 & maskB;
    S.p3.b = D==3 & maskB;
end

% Preallocate components
q = nan(1,2);

% Compute the neural similarity
R = corr(M);

% Extract mean similarities (colocation=-1)
mu1 = mean(R(S.n1));
mu2 = mean(R(S.n2));
mu3 = mean(R(S.n3));

%make all similarities 'min centred' (mu3 just becomes 0)
mu1 = mu1 -mu3;
mu2 = mu2 - mu3;

% Get the fit (colocation=-1)
if ~isnan(mu1) && (mu1 > mu2)
    q(1) = solveQ([mu1;mu2]);
end

mu1 = zeros(2,1);
mu2 = zeros(2,1);
mu3 = zeros(2,1);

% Extract mean similarities (colocation=+1)
mu1(1) = mean(R(S.p1.a));
mu2(1) = mean(R(S.p2.a));
mu3(1) = mean(R(S.p3.a));
mu1(2) = mean(R(S.p1.b));
mu2(2) = mean(R(S.p2.b));
mu3(2) = mean(R(S.p3.b));

%make all similarities 'min centred' (mu3 just becomes 0)
mu1(1) = mu1(1) - mu3(1);
mu2(1) = mu2(1) - mu3(1);
mu1(2) = mu1(2) - mu3(2);
mu2(2) = mu2(2) - mu3(2);

%mean across a-a similairties and b-b similarities
mu1 = mean(mu1);
mu2 = mean(mu2);
% Get the fit (colocation=+1)
if ~isnan(mu1) && (mu1 > mu2) 
    q(2) = solveQ([mu1;mu2]);
end

return

function [q] = solveQ(rho)

%formula for model is rho2 = rho1 * ((1/3)^q)
%q is therefore 

%rho2/rho1 = (1/3)^q
%ln(rho2/rho1) = q * ln(1/3)
%q = (ln(rho1) - ln(rho2)) / ln(3)

% Set persistent variables
persistent h 
if isempty(h) 
    h = 3; 
end

q = log((rho(2)/rho(1))) / log(h);

return
