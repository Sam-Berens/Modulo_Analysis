function [todo] = checkFmMagStatus()

checksDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);

%% Get list of Subject dirs
dirList = dir();
dirList = dirList(...
    [dirList.isdir]' & ...
    cellfun(@(s)numel(s)==8,{dirList.name}'));

%% Loopy loopy
todo = true(size(dirList,1),1);
for iSubject = 1:size(dirList,1)
    subjectId = dirList(iSubject).name;
    cd(subjectId);
    cd Fieldmap;
    fmList = dir('*_FmMag.nii');
    if ~isempty(fmList)
        todo(iSubject) = false;
    end
    cd ..;
    cd ..;
end

%% Cd back
cd (checksDir);

%% Display
todo = {dirList(todo).name}';

return