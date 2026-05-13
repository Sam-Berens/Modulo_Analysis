function [components] = qFunc1(M)
%QFUNC1 Fit q-models separately for across- and within-position similarities.
%
%   COMPONENTS = QFUNC1(M) computes a correlation-based representational
%   similarity matrix from the columns of M, averages similarities between
%   stimulus pairs separated by modular distances 1, 2, and 3, and fits the
%   non-linear distance-similarity model separately for:
%
%       colocation = -1: across-position comparisons.
%       colocation = +1: within-position comparisons.
%
%   The fitted model is:
%
%       rho_i = b0 + b1 * (1 - i/3)^q
%
%   where i is the modular distance between sparks, rho_i is the mean neural
%   similarity at that distance, b0 is an intercept, b1 is a non-negative
%   distance-related scaling term, and q controls the non-linearity of the
%   relationship between modular distance and representational similarity.
%
%   M must be an N-by-12 numeric matrix, where rows are features/voxels and
%   columns are condition-wise activity patterns. Columns 1:6 and 7:12 are
%   assumed to correspond to the same modular-value order, 0:5, in two
%   serial or presentation positions.
%
%   Comparisons between conditions with the same modular value are excluded,
%   so the fitted models are not driven by visual identity or self-similarity
%   effects.
%
%   COMPONENTS is a 1-by-8 vector:
%
%       COMPONENTS(1) = b0 for colocation = -1, across-position comparisons.
%       COMPONENTS(2) = b1 for colocation = -1, across-position comparisons.
%       COMPONENTS(3) = q  for colocation = -1, across-position comparisons.
%       COMPONENTS(4) = SSE for colocation = -1, across-position comparisons.
%
%       COMPONENTS(5) = b0 for colocation = +1, within-position comparisons.
%       COMPONENTS(6) = b1 for colocation = +1, within-position comparisons.
%       COMPONENTS(7) = q  for colocation = +1, within-position comparisons.
%       COMPONENTS(8) = SSE for colocation = +1, within-position comparisons.
%
%   A fit is only attempted for a given colocation level when the mean
%   similarities follow the strict structural ordering rho_1 > rho_2 > rho_3.
%   If this ordering does not hold, the corresponding four output values are
%   left as NaN.
%
%   Interpretation of q:
%       q = 1 gives a linear distance-similarity relationship.
%       q > 1 implies sharper similarity changes at smaller distances.
%       q < 1 implies sharper similarity changes at larger distances.
%
%   The fitting step requires FMINCON from the Optimization Toolbox.
%
%   See also CORR, FMINCON.
%

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

% Preallocate components
components = nan(1,8);

% Compute the neural similarity
R = corr(M);

% Extract mean similarities (colocation=-1)
mu1 = mean(R(S.n1));
mu2 = mean(R(S.n2));
mu3 = mean(R(S.n3));

% Get the fit (colocation=-1)
if ~isnan(mu1) && (mu1 > mu2) && (mu2 > mu3)
    components(1:4) = fitQ([mu1;mu2;mu3]);
end

% Extract mean similarities (colocation=+1)
mu1 = mean(R(S.p1));
mu2 = mean(R(S.p2));
mu3 = mean(R(S.p3));

% Get the fit (colocation=+1)
if ~isnan(mu1) && (mu1 > mu2) && (mu2 > mu3)
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