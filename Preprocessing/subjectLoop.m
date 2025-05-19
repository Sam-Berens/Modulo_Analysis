function [] = subjectLoop(func,subjectIds)

%% Check inputs
if nargin == 0
    error('subjectLoop requires at least one input argument.');
elseif nargin == 1
    dirLits = dir('../../Data');
    subjectIds = {dirLits(cellfun(@(s)numel(s)==8,{dirLits.name}')).name}';
elseif nargin > 2
    error('Too many input arguments.');
end

%% Apply the function
for iSubject = 1:numel(subjectIds)
    func(subjectIds{iSubject});
end
return