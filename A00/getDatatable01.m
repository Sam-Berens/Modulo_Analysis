function [datatable01] = getDatatable01()
if exist('Datatable01.mat','file')
    datatable01 = load('Datatable01.mat');
    return
end

bNtb = getWAIC('Binomial');
vMtb = getWAIC('vonMises');
datatable01 = join(vMtb,bNtb);
save('Datatable01.mat','datatable01');
return