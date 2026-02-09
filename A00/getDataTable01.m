function [DataTable01] = getDataTable01()
if exist('DataTable01.mat','file')
    DataTable01 = load('DataTable01.mat');
    DataTable01 = DataTable01.DataTable01;
    return
end
G = 'G1';
DataTable01 = get_pNonc(G);
DataTable01 = join(DataTable01,getWAIC(G,'Binomial'));
DataTable01 = join(DataTable01,getWAIC(G,'vonMises'));
save('DataTable01.mat','DataTable01');
return