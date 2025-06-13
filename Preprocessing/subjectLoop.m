function [] = subjectLoop(func,subjectIds,numWorkers)

% Check inputs
if nargin == 0
    error('subjectLoop requires at least one input argument.');
elseif nargin == 1
    dirLits = dir('../../Data');
    subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
elseif nargin < 3
    numWorkers = 1;
elseif nargin > 3
    error('Too many input arguments.');
end
if numWorkers > 8
    result = questdlg(sprintf([...
        'Are you sure that you want to run this job with %i workers?%c',...
        'That is quite a lot of workers...'],...
        numWorkers,10),...
        'Are you sure?',...
        'Yes, I am sure!','No, escape.',...
        'No, escape.');
    if ~sprintf(result,'Yes, I am sure!')
        return
    end
end

%% Apply the function
if numWorkers==1
    for iSubject = 1:numel(subjectIds)
        func(subjectIds{iSubject});
    end
else
    parpool(numWorkers);
    parfor iSubject = 1:numel(subjectIds)
        func(subjectIds{iSubject}); %#ok<PFBNS>
    end
end
return