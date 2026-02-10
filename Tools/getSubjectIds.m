function [subjectId] = getSubjectIds(G)
%GETSUBJECTIDS Return subject IDs for a specified group.
%
%   subjectId = GETSUBJECTIDS(G) returns a categorical array of subject IDs
%   found in the ../../Data directory, filtered according to the specified
%   group ID G.
%
%   Group IDs:
%       'G0' - All subjects found in the Data directory (default).
%       'G1' - Subset of subjects included for analysis.
%       'G2' - Subset of subjects classified as strong generalisers.
%
%   If G is not provided, the function defaults to 'G0'.
%
%   Subject IDs are defined as directory names of length 8 characters and
%   are returned in sorted order.
%
%   Input:
%       G         - Character vector specifying the group ID.
%
%   Output:
%       subjectId - Categorical array of subject IDs.
%
%   Example:
%       subjectId = getSubjectIds('G1');

% Input check
if nargin < 1
    G = 'G0';
    warning([...
    'Group ID not specified.%c',...
    'Defaulting to all subjects in the /Data directory (G0).'],10);
end

% Get the full list of subject IDS:
dirLits = dir('../../Data');
subjectId = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectId = sort(subjectId); %#ok<TRSRT>

% G-switch
switch G
    case 'G0'
        subjectId = categorical(subjectId);
        return
    case 'G1'
        toKeep = ~ismember(subjectId,{'eade18a5'});
        subjectId = subjectId(toKeep);
    case 'G2'
        pNonc = get_pNonc('G1');
        pNonc = pNonc(pNonc.pNonc>5,:);
        toKeep = ismember(subjectId,cellstr(pNonc.subjectId));
        subjectId = subjectId(toKeep);
    case 'G3'
        toKeep = ~ismember(subjectId,{'eade18a5';'efcb7c45'});
        subjectId = subjectId(toKeep);
    otherwise
        error('Unrecognised group ID.');
end
subjectId = categorical(subjectId);
return