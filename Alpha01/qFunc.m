function [components] = qFunc(M)

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