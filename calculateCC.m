function [Imax,qopt,Ijack,b] = calculateCC(Esb,k,varargin)
%CALCULATECC Calculates the channel capacities for a given dataset and k
%value, in the k-nearest neighbors density estimation algorithm
%
arg.jackknife = true;
arg.jackreplicas = 1; 
arg.festim = {}; 
arg.qopt = []; 
arg.k=1; 
arg = parseVarargin(varargin,arg);
%% Jackknife estimates of channel capacity if needed
if arg.jackknife && isempty(arg.qopt)
    jackProb = repmat(linspace(0.6,0.95,20),1,arg.jackreplicas);
    Ijack = nan(size(jackProb));
    parfor i=1:numel(jackProb)
        IX = cellfun(@(e) randsample(size(e,2),ceil(jackProb(i)*size(e,2))),Esb,'uniformoutput',0);
        Ejack = cellfun(@(e,ix) e(:,ix),Esb,IX,'uniformoutput',0);
        Ijack(i) = calculateCC(Ejack, k,'jackknife', false);
    end
    b = regress(Ijack(:),[ones(numel(jackProb),1) 1./jackProb(:)]);
    Imax = b(1); 
    qopt = []; 
    return
end
%% Obtain probability densities if they were not supplied
if isempty(arg.festim)
    Nq = numel(Esb);
    Festim = cell(Nq);
    for i=1:Nq
        Xi = Esb{i};
        parfor j=1:Nq
            if j~=i
                Festim{i,j} = kNNdensityEstimation(Esb{j},Xi,k);
            else
                Festim{i,j} = kNNdensityEstimation(Esb{j},Xi,k);
            end
        end
    end
else
    % arg.festim could be matrix Nq x Nq of estimation or a row vectpr of
    % prob distribution (or gmm) objects with the method pdf
    Nq=size(arg.festim,2); 
    if isa(arg.festim{1},'ProbDist') || isa(arg.festim{1},'gmdistribution')
        Festim = cell(Nq); 
        for i=1:Nq
            for j=1:Nq
                Festim{i,j} = arg.festim{j}.pdf(Esb{i}); 
            end
        end
        
    else % assume that arg.festim are actual densities
        Festim = arg.festim;
    end
end
%% Optimize channel capacity using fmincon
fun = @(q) -fMIfromDensities(Festim,q);
opts = optimoptions(@fmincon,'Display','Off');
q0 = ones(Nq,1)/Nq;
if ~isempty(arg.qopt)
    qopt = arg.qopt(:); 
else
    qopt= fmincon(fun,q0,[],[],ones(1,Nq),1,zeros(Nq,1),ones(Nq,1),[],opts);
end
Imax = -fun(qopt);
