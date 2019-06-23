function [CC, CCjack] = Expt6calculateAsymCCvector(dynamics)
%EXPT5CALCULATEASYMCCVECTOR Calculates the channel capacities for dynamics
%corresponding to each model obtained from Experiments 2 and 3, but using
%asymmetric vector sampling from time points of maximal variance
%
%  EXPT6CALCULATEASYMCCVECTOR(dynamics) computes the channel capacity for
%  given dynamics across different values of vectorial dimension 
%
k = 5; % Number of neighbors to use for k nearest neighbors density estimation
nrepeats = 10; % Number of replicates
CC = cell(1,10);
CCjack = cell(1,10);
jackProb = linspace(0.6,0.95,20);
njack = length(jackProb);
% Precompute variance of dynamics and then sort by time points
[~,dynamicsVariance] = cellfun(@(X) sort(var(X,[],2),'descend'),dynamics,'UniformOutput',false);
% Cycle over dimensions of multivariate vector
for d = 1:10 
    d
    %% Select d timepoints for inclusion in multivariate vector
    % Subdivide total time interval into d subintervals
    lb = 1;
    ub = 48;
    seg = round(linspace(lb,ub,d));
    % Choose the midpoint of each segment
    midpt = [];
    for tmp = 1:d-1
        midpt(tmp) = round((seg(tmp)+seg(tmp+1))/2);
    end
    frames = [midpt ub];
    %% Obtain channel capacity for given multivariate vector
    data = getEsb(dynamics,dynamicsVariance,frames);
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

function Esb = getEsb(data,dataVariance,frames)
%GETESB Obtain entries from dynamics matrix which correspond to maximum
%   variance in each of d symmetric intervals
Esb = cell(size(data));
for i = 1:numel(data)
    tmp = data{i};
    for j = 1:numel(frames)
        % Option 1: find best time points with maximum variance regardless
        % of which subinterval it belongs to (GREEDY)
%         idx = dataVariance{i}(j); 
        % Option 2: find best time points within each of d subintervals
        % (BALANCED)
        idx = dataVariance{i}(find(dataVariance{i} <= frames(j),1,'first'));
        Esb{i}(j,:) = tmp(idx,:);
        dataVariance{i}(dataVariance{i} <= frames(j)) = [];
    end
end
end