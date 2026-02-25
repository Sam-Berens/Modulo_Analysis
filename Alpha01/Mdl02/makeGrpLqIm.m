function [] = makeGrpLqIm(G)

%% Cd out
wd = pwd;
cd ..;

%get initial folder locations
dirs.Data = '../../Data' ;
dirs.Group =  fullfile(dirs.Data,'_Group','G1');
dirs.Mdl2 = fullfile(dirs.Group,'Analysis','Alpha01','Mdl02');
if ~exist(dirs.Mdl2,'dir')
    mkdir(dirs.Mdl2);
end 

%TO DO make paths not absolute 
%get hipp mask
src = '/mnt/Hestia/Modulo/Data/_Group/MniRois/';
lFn = '_RitcheyProb_lHippComb.nii';
rFn = '_RitcheyProb_rHippComb.nii';
lV = spm_vol(fullfile(src,lFn));
rV = spm_vol(fullfile(src,rFn));
lM = spm_read_vols(lV);
rM = spm_read_vols(rV);
maskM = lM + rM;

%get 'scan' names 
subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
srchLghtFn = categorical({'wlQ.nii'});
srchLghtFn = repmat(srchLghtFn,[nSubs,1]);
scans = arrayfun(@(x,y) fullfile(dirs.Data,char(x),'Analysis',...
    'Alpha01','Mdl02','G1',char(y)),subjectIds,srchLghtFn,...
    'UniformOutput',false);
%load in both colocs
v = spm_vol(scans);


% [~,~,z] = ind2sub(size(mp),s);
% slice = squeeze(lqRGB(:,:,z(1000),:));
% imshow(slice);
coloc = 'np';
for ii=1:2
    cCl = coloc(ii);
    cV = cellfun(@(x) x(ii),v);
    cM = spm_read_vols(cV);
   
    %add hipp mask
    %removed the masked values of q (some of which come from being invalid q
    %not from hipp mask
    cM = cM.*maskM;
    cM(cM==0) = nan;
    
    %make the group image and save it 
    gM = mean(cM,4,'omitmissing');
    cV = cV(1);
    fname = fullfile(wd,G,'_coloc',upper(cCl),'1.nii');
    cV.fname = fname;
    spm_write_vol(cV,gM);

    s = find(~isnan(cM));
    data_norm =(cM - min(cM(:))) / (max(cM(:)) - min(cM(:)));
    data_norm(~s) = nan;
    cmap = colormap(parula(numel(s)));
    [~,order] = sort(data_norm(s));
    rgb = cmap(order,:);

    lqR = nan(181,217,181);
    lqG = nan(181,217,181);
    lqB = nan(181,217,181);
    lqR(s) = rgb(:,1);
    lqG(s) = rgb(:,2);
    lqB(s) = rgb(:,3);
    lqRGB = cat(4,lqR,lqG,lqB);

    lqRGB = permute(lqRGB,[1,2,4,3]);
    lqRGB(isnan(lqRGB)) = 0;
    fname = fullfile(wd,G,'_coloc',cCl,'1.avi');
    v = VideoWriter(fname, 'Motion JPEG AVI');
    v.FrameRate = 30;   % Set frame rate
    open(v)
    for k = 1:size(lqRGB,4)
        frame = lqRGB(:,:,:,k);   % Extract kth RGB frame
        writeVideo(v, frame);
    end
    close(v)
end
% Cd back in
cd(wd);
return