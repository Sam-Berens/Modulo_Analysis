function [vPreds] = getVisCorr()
wd = pwd;
fName = fullfile(wd,'Mdl05','vPreds.mat');
if exist(fName,'file')
    strct = load(fName);
    vPreds = strct.vPreds;
    return
end

% Loads the DenseNet-169 model
projDir = '../..';
net = importNetworkFromONNX(fullfile(projDir, 'densenet169s.onnx'));
inputSize = net.Layers(1).InputSize;
%remove 'functional' layers
%toExlude = {'images','images_BatchSizeVerifier','x_Flatten','x_classifier_Gemm','outputOutput'};
%toInclude = {'x_features_dense_835'};
%% Code to verify we're getting the right layer:
layers  = net.Layers;
isConcat = arrayfun(@(l) isa(l,'nnet.cnn.layer.ConcatenationLayer'), layers);
cLayers = layers(isConcat);
% We are looking for conv5_block16_concat, conv5 is in dense block 4 (
% first conv stage is before denseblock 1)
idx = 6 +12 + 32 + 16;% this is n dense layers in blocks 1:3 + 16 
layers = cLayers(idx);
%layers = layers(~ismember({layers.Name},toExlude));
%layers = layers(ismember({layers.Name},toInclude));

%% Loop through the images and resize them to be the correct input size
imgFldr = fullfile(projDir,'stims');
for ii = 0:5
    imgName = sprintf('i%02d.png',ii);
    imgFn = fullfile(imgFldr,imgName);
    img = imread(imgFn);
    if  any(size(img)~= inputSize)
        img = imresize(img,inputSize(1:2));
        [~,imgName] = fileparts(imgName);
        %move them to resized folder so if we run a different network with
        %a different input dim we're not interpolating back and forth
        imwrite(img,fullfile(imgFldr,'resized',imgName));
    end
end

%% Get layer represenations
imgs = arrayfun(@(x) [imgFldr,filesep,sprintf('i%02d.png',x)],(0:5)','UniformOutput',false);
X = cell(numel(layers),1);
for ii = 1:numel(imgs)
    imgName = imgs{ii};
    cImg = imread(imgName);
    [x,y,z] = size(cImg);
    % prep for dlarray reqs (add batch dimension etc)
    cImg = reshape(cImg,x,y,z,1);
    dlImg = dlarray(single(cImg),'SSCB');
    for iL = 1:numel(layers)
        layerName = layers(iL).Name;
        x = forward(net, dlImg, Outputs=layerName);
        x = extractdata(x);
        x = squeeze(x);
        % vectorize for RSM
        x = x(:);
        if ii == 1
            X{iL} = nan(numel(x),numel(imgs));
        end
        X{iL}(:,ii) = x;
    end
end
simMats = cellfun(@(m) atanh(corr(m)),X,'UniformOutput',false);

%vPreds = firstComp(simMats);
vPreds = tanh(simMats{:,:});
save(fName,'vPreds');
return


% function [vPreds] = firstComp(simMats)
% sele = tril(true(6),-1);
% A = cellfun(@(x) x(sele)',simMats,'UniformOutput',false);
% %note that the rows are layers of the network, elements are stim-wise
% %comparisons
% A = cell2mat(A);
% %extract the layer-wise and element-wise weights, along with the first
% %singular value
% [U,S,V] = svds(A,1);
% %reconstruct the rsm from first component across layers
% Ahat = U*S*V';
% Ahat = mean(Ahat,1);
% Ahat = squareform(Ahat);
% %technically I don't think I need to do this
% Ahat(logical(eye(6))) = inf;
% vPreds = tanh(Ahat);
% return

