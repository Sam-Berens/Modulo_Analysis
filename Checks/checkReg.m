function [] = checkReg(G)
if nargin == 0
    [subjectIds] = getSubjectIds();
elseif contains(G,'G')
    [subjectIds] = getSubjectIds(G);
else 
    error('Missing input');
end 
cd ../..;
cd Data/;

for iSubject = 1:1:size(subjectIds,1)
    getAndDisplay(subjectIds{iSubject,1});
    pause();    
end
cd ../Scripts/Checks;
return

function [] = getAndDisplay(subjectId)

imgsToDisplay = cell(2,1);

subjectDir = subjectId;
cd(subjectDir);

%% Get bias corrected in FilesNames{1,1}:
cd Structural;
imgsToDisplay{1,1} = dir('m_*.nii');
imgsToDisplay{1,1} = [pwd,filesep,imgsToDisplay{1,1}.name];
cd ..; %now youre in the subject folder

%% Get mean EPI in FilesNames{2,1}:
cd EPI;
imgsToDisplay{2,1} = dir('meanu_*.nii');
imgsToDisplay{2,1} = [pwd,filesep,imgsToDisplay{2,1}.name];
cd ..;%now youre in the subject folder

%% Cd out:
cd ..;%now youre in the data folder 

%% Hand over to SPM Job man:
SpmBatch = {};
SpmBatch{1}.spm.util.checkreg.data = imgsToDisplay;
spm_jobman('initcfg');
spm_jobman('run',SpmBatch);
return