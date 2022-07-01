clc
clear

addpath IGA\
imds_dataset = imageDatastore("data/train_23x23/","LabelSource","foldernames");
imds_solution = imageDatastore("data/solution_23x23/","LabelSource","foldernames");


[train_dataset, validate_dataset] = splitEachLabel(imds_dataset,0.8);
[train_solution, validate_solution] = splitEachLabel(imds_solution,0.8);

imds_train = combine(train_dataset,train_solution);
imds_validate = combine(validate_dataset,validate_solution);

imds_train = shuffle(imds_train);
imds_validate = shuffle(imds_validate);

layers = [
    imageInputLayer([23,23,3],Name="Input")

    convolution2dLayer(5,8,Padding=[1 1 1 1],WeightsInitializer="he",Name="Conv 1")
    reluLayer(Name="ReLU 1")
    maxPooling2dLayer(2,Padding="same",Stride=2, Name="MaxPool 1")

    convolution2dLayer(5,16,Padding="same",WeightsInitializer="he",Name="Conv 2")
    reluLayer(Name="ReLU 2")
    maxPooling2dLayer(2,Padding="same",Stride=2, Name="MaxPool 2")

    convolution2dLayer(5,32,Padding="same",WeightsInitializer="he", Name="Conv 3")
    reluLayer(Name="ReLU 3")
    maxPooling2dLayer(2,Padding="same",Stride=2, Name="MaxPool 3")

    transposedConv2dLayer(2,16,Stride=2, Name="TransConv 1")
    reluLayer(Name="ReLU 4")    

    transposedConv2dLayer(2,16,Stride=2, Name="TransConv 2")
    reluLayer(Name="ReLU 5")

    transposedConv2dLayer(2,8,Stride=2, Name="TransConv 3")
    reluLayer(Name="ReLU 6")

    crop2dLayer("centercrop",Name="Crop 1")

    convolution2dLayer(1,1,Padding="same",Name="Conv 4")

    clippedReluLayer(255.0,Name="Clipped ReLU 1")
    

    regressionLayer(Name="Output")
    ];

lgraph = layerGraph(layers);

lgraph = connectLayers(lgraph, "Conv 1", "Crop 1/ref");


options = trainingOptions('adam', ...
    InitialLearnRate=0.01,...
    MaxEpochs=100, ...
    ValidationData=imds_validate, ...
    Plots="training-progress", ...
    OutputNetwork="best-validation-loss",...
    Shuffle="every-epoch",...
    Verbose=false);

[net, net_stats] = trainNetwork(imds_train,lgraph,options);