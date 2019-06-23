function f = kNNdensityEstimation(Xw,Xq,k,varargin)
%KNNDENSITYESTIMATION Computes probability densities using k-nearest
%neighbors
%
%% set input arguments; 
arg.k = k; 
arg.kdtree = []; 
arg = parseVarargin(varargin,arg); 
k = arg.k;
if size(Xq,1) == 1
    Xw = unique(Xw); % remove duplicate entries
    Xw = reshape(Xw,size(Xq,1),[]); % reshape back into 1-dimensional vector
end
[d,ni]=size(Xw); % d is dimension of trajectories, ni is number of traj
%% deal with the case that Xq is empty
if ~exist('Xq','var') || isempty(Xq) || isequal(Xw,Xq)
    Xq = Xw;
    k = k + 1; % to remove the selfie case
end

if isempty(arg.kdtree)
    % to speed things a bit can use annquery => but mex file can cause Matlab
    % to crash
    [~,kDist] = annquery(Xw,Xq,k);
elseif isa(arg.kdtree,'KDTreeSearcher')
    [~,kDist] = knnsearch(arg.kdtree,Xq','K',k); 
else
    error('If you suppply a kdtree it must be of KDTreeSearcher class!'); 
end

% Fixes any zero distance that converts to inf as reciprocal (for repeat
% pts)
if any(kDist(k,:) == 0)
    kDist(kDist==0) = min(kDist(kDist~=0));
end

kDist = kDist(k,:);

% [~,D] = knnsearch(Xw',Xq','K',k); 
% D=D(:,k); 

NDsphere = pi.^(d/2)./gamma(d/2+1); 
f = (k-1)./kDist.^d/ni/NDsphere; 
f=f(:); 

