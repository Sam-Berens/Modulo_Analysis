function [count] = checkVDMStatus()

checksDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);

%% Get list of Subject dirs
dirList = dir();
dirList = dirList(...
    [dirList.isdir]' & ...
    cellfun(@(s)numel(s)==8,{dirList.name}'));

%% Loopy loopy
count = zeros(size(dirList,1),1);
for iSubject = 1:size(dirList,1)
    subjectId = dirList(iSubject).name;
    cd(subjectId);
    cd Fieldmap;
    fileList = dir('vdm5__*_FmORI_R*.nii');
    count(iSubject) = numel(fileList);
    cd ..;
    cd ..;
end

%% Cd back
cd (checksDir);

%% Display
count = table(categorical({dirList.name}'),count,...
    'VariableNames',{'SubjectId','VDMcount'});

return