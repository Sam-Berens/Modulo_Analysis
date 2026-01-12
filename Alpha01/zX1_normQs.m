function [] = zX1_normQs()
G = 'G1';
imagePath = ['Analysis',filesep,'Alpha01',...
    filesep,'Mdl01'];
norm2MNI(G,imagePath);
return
