function [] = checkLastMod()
subjectIds = dir();
subjectIds = subjectIds(...
    [subjectIds.isdir]' & ...
    cellfun(@(s)numel(s)==8,{subjectIds.name}'));
subjectIds = {subjectIds.name}';
LastMod = nan(size(subjectIds,1),5);
for iSubject = 1:numel(subjectIds)
    cd(subjectIds{iSubject});
    cd EPI/0_Raw;
    runList = dir('R*');
    for ii = 1:size(runList,1)
        iRun = str2double(runList(ii).name(2));
        cd(runList(ii).name);
        fileList = dir('*EPI*');
        LastMod(iSubject,iRun) = max([fileList.datenum]);
        cd ..;
    end
    cd ../..;
    cd ..;
end
LastMod = datetime(LastMod,'ConvertFrom','datenum');
LastMod = array2table(LastMod,...
    'VariableNames',cellfun(@(i)sprintf('R%i',i),num2cell(1:5),...
    'UniformOutput',false),...
    'RowNames',subjectIds);
return