function [CC, CCjack] = Expt8calculatePCACCvector(dynamics)
%EXPT8CALCULATEPCACCVECTOR Calculates the channel capacities for dynamics
%corresponding to each model obtained from Experiments 2 and 3, but
%preprocesses the dataset using principal components analysis, and then
%subsamples according to vectorial dimension
%
%  EXPT8CALCULATEPCACCVECTOR(dynamics) computes the channel capacity for
%  given dynamics across different values of vectorial dimension 
%
k = 5; % Number of neighbors to use for k nearest neighbors density estimation
nrepeats = 10; % Number of replicates
CC = cell(1,10);
CCjack = cell(1,10);
jackProb = linspace(0.6,0.95,20);
njack = length(jackProb);
PCAdata = getEsbPCA(dynamics);
% Cycle over dimensions of multivariate vector
for d = 1:10 
    d
    %% Obtain channel capacity for given multivariate vector
    data = getEsb(PCAdata,d);
    I = zeros(nrepeats,1);
    Ijack = zeros(nrepeats,njack);
    for repeat = 1:nrepeats
        tic
        [I(repeat),~,Ijack(repeat,:)] = calculateCC(data,k);
        toc
    end
    I
    CC{1,d} = I;
    CCjack{1,d} = Ijack;
end
end

function Esb = getEsbPCA(data)
%GETESB Obtain entries from dynamics matrix which correspond to d principal
%   components
Esb = cell(size(data));
for i = 1:numel(data)
    % First perform z-standardization to properly prepare dataset for PCA
    tmp = zscore(data{i}');
    [~,transformData,~,~,explained] = pca(tmp,'NumComponents',10);
    Esb{i} = transformData';
end
end

function Esb = getEsb(data,t)
Esb = cell(size(data));
for i = 1:numel(data)
    tmp = data{i};
    Esb{i} = tmp(1:t,:);
end
end