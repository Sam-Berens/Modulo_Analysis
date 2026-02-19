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
toExlude = {'images','images_BatchSizeVerifier','x_Flatten','x_classifier_Gemm','outputOutput'};
layers  = net.Layers;
layers = layers(~ismember({layers.Name},toExlude));

%% Loop through the images and resize them to be the correct input size
imgFldr = fullfile(projDir,'stims');
for ii = 0:5
    imgName = fullfile(imgFldr,sprintf('i%02d.png',ii));
    img = imread(imgName);
    if  any(size(img)~= inputSize)
        img = imresize(img,inputSize(1:2));
        [~,imgName] = fileparts(imgName);
        %move them to resized folder so if we run a different network with
        %a different input dim we're not interpolating back and forth
        imwrite(img,fullfile(imgFldr,'resized',imgName));
    end
end

%% Get layer represenations
imgs = arrayfun(@(x) [imgFldr,sprintf('i%02d.png',x)],(0:5)','UniformOutput',false);
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
vPreds = firstComp(simMats);
save(fName,'vPreds');
return


function [vPreds] = firstComp(simMats)
sele = tril(true(6),-1);
A = cellfun(@(x) x(sele)',simMats,'UniformOutput',false);
%note that the rows are layers of the network, elements are stim-wise
%comparisons
A = cell2mat(A);
%extract the layer-wise and element-wise weights, along with the first
%singular value
[U,S,V] = svds(A,1);
%reconstruct the rsm from first component across layers
Ahat = U*S*V';
Ahat = mean(Ahat,1);
Ahat = squareform(Ahat);
%technically I don't think I need to do this
Ahat(logical(eye(6))) = inf;
vPreds = tanh(Ahat);
return

