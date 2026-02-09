function [] = z07_testPnonc()
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
%mean center the nonCom predictor 
t = get_pNonc(G);
x = t{:,'pNonc'} - mean(t.pNonc);

% call tfce_nullBoot_corr() to test the correlation between 
% the mean-centred pNonc predictor and zTemplate score for each voxel
% across subjects, using tfce

%% Voxel-wise correlation, no parfor
[pVal, tfceStat] = tfce_nullBoot_corr(Y,x,false);

%read the group epi mask normed to mni so you have an appropriate header
%for your pVal map
src = fullfile(dirs.Data,'_Group',G,'Structural',...
    'GrpEpiMask00','G1_GrpEpiMask00.nii');
pVname = fullfile(dirs.Data,'_Group',G,'Analysis',...
    'Alpha01','colocation=+1','tfce_pNonc.nii');
V0 = spm_vol(src);
V1 = V0;
V1.dt = [16,0];
V1.n = [1,1];
V1.fname = pVname;
V1.descrip = 'TFCE Voxel-wise correlation pVals and tfceStat';
spm_write_vol(V1,pVal);
V1.n = [1,2];
spm_write_vol(V1,tfceStat);

% Cd back in
cd(wd);
save('tfceResult_pNonc.mat',"tfceStat","pVal");
return
