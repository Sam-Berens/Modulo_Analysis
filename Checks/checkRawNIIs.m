function [countsTable,discountsTable] = checkRawNIIs()

cd(['..',filesep,'..',filesep,'Data']);

%% Get list of Subject dirs
dirList = dir();
dirList = dirList(...
    [dirList.isdir]' & ...
    cellfun(@(s)numel(s)==8,{dirList.name}'));

%% Loopy loopy
for iSubject = 1:size(dirList,1)
    
    subjectId = dirList(iSubject).name;
    cd(subjectId);
    
    niiCounts = struct();
    niiDiscounts = struct();
    niiCounts.SubjectId = subjectId;
    niiDiscounts.SubjectId = subjectId;
    
    cd(['EPI',filesep,'0_Raw']);
    [niiCounts,niiDiscounts] = getNiiList(niiCounts,niiDiscounts);
    cd(['..',filesep,'..']);
    
    cd Fieldmap;
    [niiCounts,niiDiscounts] = getNiiList(niiCounts,niiDiscounts);
    cd ..;
    
    cd Structural;
    [niiCounts,niiDiscounts] = getNiiList(niiCounts,niiDiscounts);
    cd ..;
    
    cd ..;
    if iSubject == 1
        countsTable = struct2table(niiCounts);
        discountsTable = struct2table(niiDiscounts);
    else
        countsTable = outerjoin(countsTable,struct2table(niiCounts),...
            'MergeKeys',true);
        discountsTable = outerjoin(discountsTable,struct2table(niiDiscounts),...
            'MergeKeys',true);
    end
end

return

function [niiCounts,niiDiscounts] = getNiiList(niiCounts,niiDiscounts)
dirList = dir();
dirList = dirList(cellfun(@(s)~strcmp('.',s(1)),{dirList.name}'));
if size(dirList,1)==0
    return
end
[~,seriesName] = fileparts(dirList(1).folder);
if contains(seriesName,'SBref') || contains(seriesName,'SBRef')
    return
end
imgNums = [];
for iDir = 1:size(dirList,1)
    if dirList(iDir).isdir
        cd(dirList(iDir).name);
        [niiCounts,niiDiscounts] = getNiiList(niiCounts,niiDiscounts);
        cd ..;
    elseif strcmp(dirList(iDir).name(end-3:end),'.nii')
        if strcmp(dirList(iDir).name(1),'_')
            % Converted
            imgNums = [imgNums;
                str2double(dirList(iDir).name(end-6:end-4))]; %#ok<*AGROW>
        else
            % Not converted
            imgNums = [imgNums;
                str2double(dirList(iDir).name(end-10:end-6))];
        end
    end
end
n = numel(imgNums);
if n > 0
    niiCounts.(seriesName) = n;
    if n > 1
        niiDiscounts.(seriesName) = var(diff(imgNums))~=0;
    else
        niiDiscounts.(seriesName) = false;
    end
end
return