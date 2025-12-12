function [] = Searchlight(G)
dirs.Data = '../../Data';

subjectIds = getSubjectIds(G);
nSubs = numel(subjectIds);
roiId = 'lHi'; %(l?/r? one at a time?)

%% Loop through subjects
for iSub=1:nSubs
    cId = subjectIds{iSub};
    dirs.sub = [dirs.Data,filesep, cId];
    dirs.EPI = [dirs.sub,filesep,'EPI'];
    dirs.maskPath = [dirs.EPI,filesep,'_',cId,'_epiMask00.nii'];
    dirs.Alpha01 = [dirs.sub,filesep,'Analysis',filesep,'Alpha01'];
    %% to do check if this is feasible
    % load in EPI mask for each person to get the neighbourhood for all
    % possible searchlights:
    
    ds = cosmo_fmri_dataset(dirs.maskPath);
    %in ds, ds.fa  contains .i .j and .k which are the xyz coords in voxel
    %space of each voxel

    %... so you end up with samples having a row for each #stim

    %maybe we like call this function on one row of the dataset just to
    %get the nbrhood coords
    nbrhood = cosmo_spherical_neighborhood(ds,'radius',3);

    %nbrhood.neighbors{5}
    %Might be: [3, 4, 5, 6, 7] %these are the lin indicies of voxels inside the
    %search light.

    %so in order to mask out the nierghboods with what is inside the roi, we
    %need the linear indices of what is in the roi - that is exactly what
    %isEPIinROI is!!

    %Data coming out of this is [12,nVox,nPos]
    [Data,isEpiInRoi] = getTpatterns_EpiRes(G,cId,roiId,dirs);

    %Select only the searchlight locations inside the Hippocampus
    %remember they are specified in the voxel linear indices
    idxs = cellfun(@(x) any(ismember(x,isEpiInRoi)), nbrhood.neighbors);
    NB.neighbors = nbrhood.neighbors(idxs);
    NB.origin =  nbrhood.origin.fa.j(idxs); %we only need the y coords

    nSearchlights = numel(NB.neighbors);
    Q = nan(nSearchlights,1);
    Y = nan(nSearchlights,1);
    %% LOOP through searchlights
    for iSL = 1:numel(nSearchlights)

        %then if yes, then we have to mask out any of the neirbhood that might be on the edge
        cNB = NB.neighbors{iSL};
        %get the centre coords of this NB
        mmY = NB.origin(iSL);

        %remove any members of a specific search light which are outside the roi
        trimedNB = cNB(ismember(cNB,isEPIinROI));

        %% TODO: collect up Tstats for voxels in trimedNB
        %loop for all #stim
        cM = spm_read_vols(t.V);
        cM = cM(trimedNB); %select searchlight voxels
        cM = cM(:); %flatten
        %make M such that its [nvoxels,nStims]
        %D = pdist(M,'euclidean'); %this function aslways compares rows
        D = corr(M);
        lower = tril(true(6),-1);
        D = D(lower); %these are ordered the same as nchoose im pretty sure

        %we will compute the arithmetic mean of all similarity scores (𝑟𝑥,𝑦) that have the same predicted value according to our structural hypothesis
        nrows =6;
        f = @(x,y)min(mod((x-y),6),mod((y-x),6));
        distPairs = nchoosek(1:nrows,2);
        x = distPairs(:,1);
        y = distPairs(:,2);
        distBin = f(x,y); %now we have a value of the expected modular distance for each value of D,
        %these will be the i term in the model:

        %𝜌𝑖=𝑏0+𝑏1(1−𝑖3)𝑞
        %pi is the mean simiarity scores for 𝜌1, 𝜌2, and 𝜌3 where i is the distance (so we have to collect up each bin's D values and mean them
        P = nans(3,1);
        I = nans(3,1);
        for ii=1:3
            idxs = find(distBin==ii);
            %PUT p into the design Matrix for searchlight location jj
            P(ii,1) = mean(D(idxs));
            I(ii,1) = ii;
        end

        %transform I to get our X predictor value
        X = 1-(I/3);
        designMat = table(P,X);
        modelfun = @(b,x) b(1) + b(2)*(x^b(3));
        beta0 = [0,0,0]; %these are the initial starting values for all the free parameters estimated by the model
        mdlObj = fitnlm(designMat,modelfun,beta0);
        Q(iSL,1) = mdlObj.Coefficients{3,2};
        Y(iSL,1) = mmY;

        %so after you have estimated the 3 betas for predicting similartiy scores
        %for each searchlight location you then basically correlate all q values estimated by the model i.e. b(3), with the A-P MNI coordinate y (basically the y coordinate of the searchlight's centre, as transformed into mm space)


        %𝜔= tanh−1(corr(𝑦𝑖,ln(𝑞𝑖)))


    end

    % and then test 𝜔 against 0 across all participants;

end
return


