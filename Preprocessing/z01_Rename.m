function [] = z01_Rename(subjectId)

scriptsDir = pwd;
cd(['..',filesep,'..',filesep,'Data']);
% Cd into subject dir
cd(subjectId);

%% EPIs
cd(['EPI',filesep,'0_Raw']);
runDirs = dir('R*');
for rr = 1:numel(runDirs)
    cd(runDirs(rr).name);
    renameEPIs(subjectId);
    mkdir Json;
    movefile('*.json','Json');
    cd ..;
end
cd(['..',filesep,'..']);

%% FMs
cd('Fieldmap');
apList = {dir('FieldMap_AP_0*').name}';
paList = {dir('FieldMap_PA_0*').name}';
nFms = numel(apList);
if nFms == 1
    cd(dir('FieldMap_AP_0*').name);
    renameFMs(subjectId,'AP');
    movefile('*.nii','..');
    cd('..');
    cd(dir('FieldMap_PA_0*').name);
    renameFMs(subjectId,'PA');
    movefile('*.nii','..');
    cd('..');
    mkdir('SBref+Json')
    cellfun(@(s) movefile(s,'SBref+Json'), {dir('FieldMap*').name}'); 
    cd('..');
elseif nFms > 1 && size(paList,1) == nFms
    for ii = 1:nFms
        cd(apList{ii,1});
        renameFMs(subjectId,'AP');
        movefile('*.nii','..');
        cd('..');
        cd(paList{ii,1});
        renameFMs(subjectId,'PA');
        movefile('*.nii','..');
        cd('..');
    end
    mkdir('SBref+Json')
    cellfun(@(s) movefile(s,'SBref+Json'), {dir('FieldMap*').name}');
    cd('..');
end   
%% Structural
cd('Structural')
renameT1(subjectId)
mkdir Json;
movefile('*.json','Json');
movefile('*.mat','Json');
cd('..');

%% Cd out of subject dir
cd ..;
cd(scriptsDir);

return

function [] = renameEPIs(subjectId)
%% Get SeriesNumber
seriesN = getSeriesNumber();
%% Get the image numbers
epiList = dir('*.nii');
imgNum = cellfun(@(s)str2double(s((end-10):(end-6))),{epiList.name}');
imgCount = numel(imgNum);
%% Construct the target file names
fnc = @(sId,sN,iN) sprintf('_%s_EPI-%02d-%03d.nii',sId,sN,iN);
targetFns = cellfun(fnc,...
    repmat({subjectId},imgCount,1),...
    repmat({seriesN},imgCount,1),...
    num2cell(imgNum),...
    'UniformOutput',false);
%% Rename the files
for iImg = 1:numel(epiList)
    movefile(epiList(iImg).name, targetFns{iImg});
end
return

function [] = renameFMs(subjectId,direction)
%% Get SeriesNumber
seriesN = getSeriesNumber();
%% Get the image numbers
epiList = dir('*.nii');
imgNum = cellfun(@(s)str2double(s((end-10):(end-6))),{epiList.name}');
imgCount = numel(imgNum);
%% Construct the target file names
fnc = @(sId,sN,iN) sprintf('_%s_Fm%s-%02d-%03d.nii',sId,direction,sN,iN);
targetFns = cellfun(fnc,...
    repmat({subjectId},imgCount,1),...
    repmat({seriesN},imgCount,1),...
    num2cell(imgNum),...
    'UniformOutput',false);
%% Rename the files
for iImg = 1:numel(epiList)
    movefile(epiList(iImg).name, targetFns{iImg});
end
return

function [seriesN] = getSeriesNumber()
fileList = dir('*.json');
text = fileread(fileList(1).name);
data = jsondecode(text);
seriesN = data.acqpar.SeriesNumber;
return

function [] = renameT1(subjectId)
%% Get SeriesNumber
seriesN = getSeriesNumber();
%% Construct the target file name and move
targetFn = sprintf('_%s_T1-%02d.nii',subjectId,seriesN);
movefile(dir('*.nii').name,targetFn)
return