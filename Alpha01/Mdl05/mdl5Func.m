function [z] = mdl5Func(M,imgPerm)
persistent nanzscore VisSim prev_imgPerm S IX;
if isempty(nanzscore) || isempty(VisSim)
    nanzscore = @(v) (v-mean(v,'omitmissing')) ./ std(v,'omitmissing');

    temp = load('VisualSim.mat');
    VisSim = temp.VisualSim;
end

if isempty(prev_imgPerm) || isempty(S) || isempty(IX) || ...
        ~isequal(prev_imgPerm, imgPerm)

    % Construct the basic 6x6 similarity hypothesis
    simFun = @(x,y) 1-((min(mod((x-y),6),mod((y-x),6))/3));
    H_ = nan(6);
    for ii = 1:36
        [x,y] = ind2sub([6,6],ii);
        H_(ii) = simFun(x,y);
    end

    % Hypothesis for colocation=-1
    Hn = kron(ones(2),H_);
    Hn(logical(kron([1,1;0,1],ones(6)))) = NaN;
    Hn(logical(kron([0,0;1,0],eye(6)))) = NaN;

    Haa = kron(ones(2),H_);
    Haa(logical(kron([0,1;1,1],ones(6)))) = NaN;
    Haa(triu(true(12))) = NaN;

    Hbb = kron(ones(2),H_);
    Hbb(logical(kron([1,1;1,0],ones(6)))) = NaN;
    Hbb(triu(true(12))) = NaN;

    Hp = sum(cat(3,Haa,Hbb),3,'omitmissing');
    Hp(isnan(Haa) & isnan(Hbb)) = NaN;

    % Make the selectors for the real data RSM
    S.n = ~isnan(Hn);
    S.aa = ~isnan(Haa);
    S.bb = ~isnan(Hbb);
    S.p = ~isnan(Hp);
    S.all = S.n | S.p;

    % Load in predicted visual similarity from denseNet-169
    V = VisSim(imgPerm,imgPerm);
    V = kron(ones(2),V);

    % Construct the design matrix
    X = [...
        nanzscore(Hn(S.all)),...
        nanzscore(Hp(S.all)),...
        nanzscore(V(S.all))];

    % Replace N/A values with mean
    X(isnan(X)) = 0;

    % Invert for regression
    IX = pinv(X);

    % Set/update prev_imgPerm
    prev_imgPerm = imgPerm;

end

% Get the 12x12 RSM for the searchlight image
R = corr(M);

% Zscore within aa/bb/ba pairs
R(S.n) = nanzscore(R(S.n));
R(S.aa) = nanzscore(R(S.aa));
R(S.bb) = nanzscore(R(S.bb));
y = R(S.all);

% Do regression and Fischer-transform the results
z = atanh(IX*y);
return