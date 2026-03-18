function [imgPerm] = getImgPerm(subjectId)

% Check inputs
requestedSingleSubject = false;
requestedG0 = false;
if nargin < 1
    subjectId = getSubjectIds('G0');
    subjectId = cellstr(subjectId);
    requestedG0 = true;
elseif ischar(subjectId) || isstring(subjectId)
    requestedSingleSubject = true;
    subjectId = {subjectId};
elseif iscategorical(subjectId)
    subjectId = cellstr(subjectId);
end

% Check whether the data have been saved
if exist('imgPerm.mat','file')
    temp = load('imgPerm.mat');
    imgPerm = temp.imgPerm;
    clear temp;
    s = contains(cellstr(imgPerm.subjectId),subjectId);
    imgPerm = imgPerm(s,:);
    if requestedSingleSubject
        imgPerm = imgPerm.imgPerm{1};
    end
    return
end

% Query the database and populate imgPerm
imgPerm = cell(size(subjectId));
for iSubject = 1:numel(subjectId)
    serverResponse = webwrite(...
        'https://b01.learningandinference.org/GetSessionDef.php',...
        'SubjectId',subjectId{iSubject});
    imgPerm{iSubject} = serverResponse.ImgPerm;
end

% If a single subject was requested, return a vector ...
if requestedSingleSubject
    imgPerm = imgPerm{1};
    return
end

% ... otherwise, return a table
subjectId = categorical(subjectId);
imgPerm = table(subjectId,imgPerm);

% Save the future
if requestedG0
    save('imgPerm.mat','imgPerm');
end

return