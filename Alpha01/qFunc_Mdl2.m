function [components] = qFunc_Mdl2(M)

% Set mean similarity selectors
persistent S;
if isempty(S1) || isempty(S2) || isempty(S3)
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

% Preallocate components
components = nans(1,8);

% Compute the neural similarity
R = corr(M);

% Extract mean similarities (colocation=-1)
mu1 = mean(R(S.n1));
mu2 = mean(R(S.n2));
mu3 = mean(R(S.n3));

% Get the fit (colocation=-1)
if ~(isnan(mu1) || (mu1 < mu2) || (mu2 < mu3))
    components(1:4) = fitQ([mu1;mu2;mu3]);
end

% Extract mean similarities (colocation=+1)
mu1 = mean(R(S.p1));
mu2 = mean(R(S.p2));
mu3 = mean(R(S.p3));

% Get the fit (colocation=+1)
if ~(isnan(mu1) || (mu1 < mu2) || (mu2 < mu3))
    components(5:8) = fitQ([mu1;mu2;mu3]);
end

return

function [components] = fitQ(rho)

% Set persistent variables
persistent h predFnc A c p0 opts;
if isempty(h) || isempty(predFnc)
    h = [2;1;0]./3; % Assumes [dist1;dist2;dist3]
    predFnc = @(p) p(1) + p(2).*(h.^p(3));
end
if isempty(A) || isempty(c) || isempty(p0) || isempty(opts)
    A = [ 
        0,-1,0;
        0,0,-1;
        0,0,1];
    c = [0;0;7];
    p0 = [0.2;1;1];
    opts = optimoptions('fmincon',...
        'Display','off','Algorithm','sqp');
end

% Run the solver
costFnc = @(p) sum((predFnc(p) - rho).^2);
[pHat,err] = fmincon(costFnc,p0,A,c,[],[],[],[],[],opts);
components = [pHat',err];
return