function [subjectId] = getSubjectIds(G)

if nargin < 1
    G = 'G0';
    warning('Group ID not specified. Defaulting to all subjects in the /Data directory (G0).');
end

% Get the full list of subject IDS:
dirLits = dir('../../Data');
subjectId = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectId = sort(subjectId); %#ok<TRSRT>

switch G
    case 'G0'
        return
    case 'G1'
        toKeep = ~ismember(subjectId,{'eade18a5'});
        subjectId = subjectId(toKeep);
    otherwise
        error('Unrecognised group ID.');
end

subjectId = categorical(subjectId);
return