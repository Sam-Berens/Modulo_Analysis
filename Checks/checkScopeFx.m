function [] = checkScopeFx()

checksDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);

%% Get list of Subject dirs
dirList = dir();
dirList = dirList(...
    [dirList.isdir]' & ...
    cellfun(@(s)numel(s)==8,{dirList.name}'));

%% Loopy loopy
problem = true(size(dirList,1),1);
for iSubject = 1:size(dirList,1)
    subjectId = dirList(iSubject).name;
    cd(subjectId);
    cd Fieldmap;
    movPar = dir('*_FmTopupCoefs-Movpar.txt');
    rawImg = dir('_*_FmPA-*-001.nii');
    try
        problem(iSubject) = rawImg.datenum > movPar.datenum;
    catch
        problem(iSubject) = false; % Not everyone has been done yet
        disp(subjectId);
    end
    cd ..;
    cd ..;
end

%% Cd back
cd (checksDir);

%% Display
disp({dirList(problem).name}');

return