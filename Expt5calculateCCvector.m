function [CC, CCjack] = Expt5calculateCCvector(dynamics)
%EXPT5CALCULATECCVECTOR Calculates the channel capacities for dynamics
%corresponding to each model obtained from Experiments 2 and 3
%
%  EXPT5CALCULATECCVECTOR(dynamics) computes the channel capacity for
%  given dynamics across different values of vectorial dimension 
%
k = 5; % Number of neighbors to use for k nearest neighbors density estimation
nrepeats = 10; % Number of replicates
CC = cell(1,10);
CCjack = cell(1,10);
jackProb = linspace(0.6,0.95,20);
njack = length(jackProb);
% Cycle over dimensions of multivariate vector
for d = 1:10 
    d
    %% Select d timepoints for inclusion in multivariate vector
    lb = 1;
    ub = 48;
    seg = round(linspace(lb,ub,d+1));
    % Choose the midpoint of each segment
    midpt = zeros(1,d);
    for tmp = 1:d
        midpt(tmp) = round((seg(tmp)+seg(tmp+1))/2);
    end
    frames = midpt;
    %% Obtain channel capacity for given multivariate vector
    data = getEsb(dynamics,frames);
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

function Esb = getEsb(data,t)
Esb = cell(size(data));
for i = 1:numel(data)
    tmp = data{i};
    Esb{i} = tmp(t,:);
end
end