function [] = z04_normSearchlight(G)
wd = pwd;
cd ..
norm2MNI(G,'Analysis/Alpha01/Searchlight',[0,0,0],[2,2,2]); 
cd(wd);
return