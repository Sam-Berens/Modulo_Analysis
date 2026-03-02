function [] = z04_normEpsilon()
% Cd out
wd = pwd;
cd ..;
G = 'G1';
imagePath = ['Analysis',filesep,'Alpha01',...
    filesep,'Mdl03'];
norm2MNI(G,imagePath,[0,0,0],[nan,nan,nan]);
% Cd back in
cd(wd);
return
