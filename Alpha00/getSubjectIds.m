function [subjectIds] = getSubjectIds(G)

if nargin < 1
    G = 'G0';
    warning('Group ID not specified. Defaulting to all subjects in the /Data directory (G0).');
end

% Get the full list of subject IDS:
dirLits = dir('../../Data');
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';

switch G
    case 'G0'
        return
    case 'G1'
        toKeep = ~ismember(subjectIds,{'eade18a5'});
        subjectIds = subjectIds(toKeep);
    case '25bc1838'
        subjectIds = {'25bc1838'};
    case 'special'
        subjectIds = {'7ca435e9';
            '7db10ad2';
            '802dbd4b';
            '8b6a0a9d';
            '8c5da281';
            '8e86fbfd';
            '9224122e';
            '93bb6a53';
            'a3b38616';
            'aad7f417';
            'b1a7739a';
            'b777a363';
            'b8282457';
            'bc74ee19';
            'c67b7f77';
            'c985b5bd';
            'cc7a3caa';
            'd0d844d6';
            'd90d8894';
            'e64c5393';
            'eec99b44';
            'efcb7c45'};
    otherwise
        error('Unrecognised group ID.');
end

return