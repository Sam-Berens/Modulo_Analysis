function [] = zX5_normE_mdl3()
G = 'G1';
imagePath = ['Analysis',filesep,'Alpha01',...
    filesep,'Mdl03'];
norm2MNI(G,imagePath,[0,0,0],[nan,nan,nan]);
return
