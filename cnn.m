clc
clear

addpath IGA\
imds_dataset = imageDatastore("data/train/","LabelSource","foldernames");
imds_solution = imageDatastore("data/solution/","LabelSource","foldernames");

%imcrop(imds_solution,[2 2 6 6])

[train_dataset, validate_dataset] = splitEachLabel(imds_dataset,0.8);
[train_solution, validate_solution] = splitEachLabel(imds_solution,0.8);

% train_dataset = transform(train_dataset,@im2double);
% validate_dataset = transform(validate_dataset,@im2double);
% train_solution = transform(train_solution,@im2double);
% validate_solution = transform(validate_solution,@im2double);

imds_train = combine(train_dataset,train_solution);
imds_validate = combine(validate_dataset,validate_solution);

imds_train = shuffle(imds_train);
imds_validate = shuffle(imds_validate);


% imageLayer = imageInputLayer([8,8,3]);

layers = [
    imageInputLayer([8,8,3],Name="Input")

    convolution2dLayer(3,8,Padding=[0 0 0 0],WeightsInitializer="he",Name="Conv 1")
    reluLayer(Name="ReLU 1")
    maxPooling2dLayer(2,Padding="same",Stride=2, Name="MaxPool 1")

    convolution2dLayer(3,16,Padding="same",WeightsInitializer="he",Name="Conv 2")
    reluLayer(Name="ReLU 2")
    maxPooling2dLayer(3,Padding="same",Stride=3, Name="MaxPool 2")

    transposedConv2dLayer(3,8,Stride=3, Name="TransConv 1")
    reluLayer(Name="ReLU 3")

    transposedConv2dLayer(2,4,Stride=2, Name="TransConv 2")
    reluLayer(Name="ReLU 4")

    convolution2dLayer(1,1,Padding="same",WeightsInitializer="he",Name="Conv 3")
    clippedReluLayer(255.0,Name="Clipped ReLU 1")

    regressionLayer(Name="Output")
    ];

% decodingLayers = [
%     transposedConv2dLayer(2,32,Stride=2), ...
%     reluLayer, ...
%     transposedConv2dLayer(2,16,Stride=2), ...
%     reluLayer, ...
%     transposedConv2dLayer(2,8,Stride=2), ...
%     reluLayer, ...
%     convolution2dLayer(1,1,Padding="same")
%     clippedReluLayer(1.0)
%     regressionLayer
%     ]; 

%layers = encodingLayers;%[imageLayer, encodingLayers];

options = trainingOptions('adam', ...
    InitialLearnRate=0.01,...
    MaxEpochs=100, ...
    ValidationData=imds_validate, ...
    Plots="training-progress", ...
    OutputNetwork="best-validation-loss",...
    Shuffle="every-epoch",...
    Verbose=false);

[net, net_stats] = trainNetwork(imds_train,layers,options);