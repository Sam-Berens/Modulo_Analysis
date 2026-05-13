function [components] = qFunc0(M)
%FITQ Fit the modular-distance q-model to three mean similarity values.
%
%   COMPONENTS = FITQ(RHO) fits the model:
%
%       rho_i = b0 + b1 * (1 - i/3)^q
%
%   to the three mean similarity values in RHO using constrained least
%   squares. RHO must be a 3-element vector ordered as:
%
%       RHO(1) = mean similarity at modular distance 1.
%       RHO(2) = mean similarity at modular distance 2.
%       RHO(3) = mean similarity at modular distance 3.
%
%   COMPONENTS is a 1-by-4 vector:
%
%       COMPONENTS(1) = b0, model intercept.
%       COMPONENTS(2) = b1, non-negative distance-related scaling term.
%       COMPONENTS(3) = q, fitted non-linearity parameter.
%       COMPONENTS(4) = SSE, sum of squared fitting error.
%
%   The optimisation constrains b1 >= 0 and 0 <= q <= 7. The intercept b0
%   is unconstrained. The fitted q parameter describes non-linearities in
%   the relationship between modular distance and representational
%   similarity.
%
%   This helper is called by QFUNC0 and QFUNC1 and requires FMINCON from the
%   Optimization Toolbox.
%
%   See also FMINCON, QFUNC0, QFUNC1.
%

% Set mean similarity selectors
persistent S1 S2 S3;
if isempty(S1) || isempty(S2) || isempty(S3)
    [A,B] = meshgrid(0:5,0:5);
    D = min(mod(A-B,6),mod(B-A,6));
    D = kron(ones(2),D);
    D(triu(true(12))) = NaN;
    S1 = D==1;
    S2 = D==2;
    S3 = D==3;
end

% Compute the neural similarity
R = corr(M);

% Extract mean similarities
mu1 = mean(R(S1));
mu2 = mean(R(S2));
mu3 = mean(R(S3));

% Return if the target inequality does not hold  
if isnan(mu1) || (mu1 < mu2) || (mu2 < mu3)
    components = NaN;
    return
end

% Get the fit
components = fitQ([mu1;mu2;mu3]);
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