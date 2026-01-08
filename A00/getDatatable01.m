function [datatable01] = getDatatable01()
if exist('Datatable01.mat','file')
    datatable01 = load('Datatable01.mat');
    return
end

vMtb = getWAIC('vonMises');
bNtb = getWAIC('Binomial');
datatable01 = innerjoin(vMtb,bNtb);
save('Datatable01.mat','datatable01');
return