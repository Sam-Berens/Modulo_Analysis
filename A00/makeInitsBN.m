
function makeInitsBN(outDir,baseName)
% makeInits
% adapted from ...
% Sam Berens (s.berens@sussex.ac.uk)
% 31/07/2025
%
% Syntax: makeInits([outDir],[baseName])
%
% Description:
%    Generate 8 JSON init files for the hspvm Stan model.
%    The function populates the following parameters:
%       alpha1  (vector[nTypes])          unbounded
%       alpha2  (vector<lower=0>[nTypes]) positive
%       beta1   (vector[nTypes])          unbounded
%       beta2   (vector<lower=0>[nTypes]) positive
%       zau     (vector[nPairs])          unbounded
%       zb      (vector[nPairs])          unbounded
%
% Optional arguments:
%    outDir   - Output directory [default './Inits'];
%    baseName - Prefix for the filenames  [default 'init'];
%
%
% Example
% -------
%    makeInits();
%    %% Then call CmdStan:
%    ./hspvm sample num_samples=1000 num_warmup=1000 \
%       data file=data.json \
%       init=init{{id}}.json \
%       id={{id}} output file=out_chain{{id}}.csv
%
%    (GNU parallel users can replace {{id}} with ::: 1-8.)
%

arguments
    outDir (1,:) char = './Inits/Binomial';
    baseName (1,:) char = 'Init';
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end

nTypes = 4;
nPairs = 36;
nChains = 8;

alpha1 = [
    -1,-0.75,-0.25, 0;
    -1,-0.75,-0.25, 0;
    -1,-0.75,-0.25, 0;
    -1,-0.75,-0.25, 0;
    0, -0.25, -0.75,-1;
    0, -0.25, -0.75,-1;
    0, -0.25, -0.75,-1;
    0, -0.25, -0.75,-1];
alpha2 = [
    0.95,0.975,1.025,1.05; 
    0.95,0.975,1.025,1.05;
    1.05,1.025,0.975,0.95;
    1.05,1.025,0.975,0.95;
    0.95,0.975,1.025,1.05;
    0.95,0.975,1.025,1.05;
    1.05,1.025,0.975,0.95;
    1.05,1.025,0.975,0.95];
beta1 = [zeros(nChains,nTypes)];
beta2 = [
    0.048,0.049, 0.051, 0.052;
    0.052, 0.051, 0.049, 0.048;
    0.048, 0.049, 0.051, 0.052;
    0.052, 0.051, 0.049, 0.048
    0.048,0.049, 0.051, 0.052;
    0.052, 0.051, 0.049, 0.048
    0.048,0.049, 0.051, 0.052;
    0.052, 0.051, 0.049, 0.048];

zau = zeros(nChains,nPairs);
zb = zeros(nChains,nPairs);

for iChain = 1:nChains
    init = struct;
    init.alpha1 = alpha1(iChain,:)';
    init.alpha2 = alpha2(iChain,:)';
    init.beta1 = beta1(iChain,:)';
    init.beta2 = beta2(iChain,:)';
    init.zau = zau(iChain,:)';
    init.zb = zb(iChain,:)';

    % Write JSON
    fname = fullfile(outDir, sprintf('%s_%d.json', baseName, iChain));
    fid = fopen(fname,'w');
    if fid < 0
        error('Could not open %s for writing.', fname);
    end
    fprintf(fid,'%s',jsonencode(init,'PrettyPrint',true));
    fclose(fid);
end

return
