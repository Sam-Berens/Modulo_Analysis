function [] = z06_testIntercept()
% Cd out
wd = pwd;
cd ..;

dirs.Data = '../../Data';
G = 'G1';
subjectIds = getSubjectIds(G);
%make Y a cell array of NIfTI file names
Y = arrayfun(@(x) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Searchlight',G,'wzTemplate_colocation=+1.nii'),...
    subjectIds,'UniformOutput', false);
%call tfce_nullBoot_oneSample() to test intercept against null across
%people using tfce

%read in y as 4d and then mask it out

%% One-sample test (H0: mean = 0), false just means no parfor
[pVal, tfceStat] = tfce_nullBoot_oneSample(Y,0,false);

%read the group epi mask normed to mni so you have an appropriate header
%for your pVal map
src = fullfile(dirs.Data,'_Group',G,'Structural',...
    'GrpEpiMask00','G1_GrpEpiMask00.nii');
pVname = fullfile(dirs.Data,'_Group',G,'Analysis',...
    'Alpha01','colocation=+1','tfce.nii');
V0 = spm_vol(src);
V1 = V0;
V1.dt = [16,0];
V1.n = [1,1];
V1.fname = pVname;
V1.descrip = 'TFCE One-Sample t-test pVals and tfceStat';
spm_write_vol(V1,pVal);
V1.n = [1,2];
spm_write_vol(V1,tfceStat);

% Cd back in
cd(wd);
save('tfceResult.mat',"tfceStat","pVal");
return
