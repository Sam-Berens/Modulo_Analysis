function [] = revertAnalysis()
confirm = questdlg('Are you sure you wish to revert?', ...
	'Confirm Revert', ...
	'YES (delete)','NO (ESCAPE)','NO (ESCAPE)');
if ~strcmp(confirm,'YES (delete)')
    return
end
%%
dirLits = dir('../../Data');
subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
subjectIds = categorical(subjectIds);
for iSubject = 1:numel(subjectIds)
    sId = subjectIds(iSubject);
    pathToTaskIO = sprintf('..%s..%sData%s%s%sBehavioural%s',...
        filesep,filesep,filesep,char(sId),filesep,filesep);

    %% A00/
    system(['rm -R ',pathToTaskIO,'A00']);

    %% TrainTaskIO.mat
    % [~,lsRes] = system(sprintf('ls -l %s | grep TrainTaskIO.mat',...
    %     pathToTaskIO));
    % if contains(lsRes,'TrainTaskIO.mat')
    %     delete([pathToTaskIO,'TrainTaskIO.mat']);
    % end
    % 
    %% A00.mat
    % [~,lsRes] = system(sprintf('ls -l %s | grep A00.mat',...
    %     pathToTaskIO));
    % if contains(lsRes,'A00.mat')
    %     delete([pathToTaskIO,'A00.mat']);
    % end
end
return