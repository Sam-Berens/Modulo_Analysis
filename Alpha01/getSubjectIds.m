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
    case '7ca435e9'
        subjectIds = {'7ca435e9'};
    case '~25bc1838'
        toKeep = ~ismember(subjectIds,{'eade18a5'});
        subjectIds = subjectIds(toKeep);
        toKeep = ~ismember(subjectIds,{'25bc1838'});
        subjectIds = subjectIds(toKeep);
    otherwise
        error('Unrecognised group ID.');
end

return