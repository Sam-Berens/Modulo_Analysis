function [] = zX4_normQcoLoc1()
G = 'G1';
imagePath = ['Analysis',filesep,'Alpha01',...
    filesep,'Mdl02'];
norm2MNI(G,imagePath,[0,0,0],[nan,nan,nan]);
return
