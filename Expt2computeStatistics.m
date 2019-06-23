function [AUC, maxFoldChange, adapt, maxResponse, maxResponseTime] = Expt2computeStatistics(dynamics)
%EXPT2COMPUTESTATISTICS Computes relevant descriptors for
%  Experiment 2, in which various dynamics produced across noise levels and
%  replications
%
%  EXPT2COMPUTESTATISTICS(dynamics) computes the summary descriptors for 
%  the cell dynamics produced across noise levels and replications
%
tspan = 1:48;
GFAPinit = 0.1834;
[AUC, maxFoldChange, adapt, maxResponse, maxResponseTime] = deal(cell(4,4));
for i = 1:size(dynamics,1)
   for j = 1:size(dynamics,2)
      for k = 1:size(dynamics{i,j},2)
            AUC{i,j}(1,k) = trapz(tspan,dynamics{i,j}(:,k));
            [maxResponse{i,j}(1,k),maxResponseTime{i,j}(1,k)] = max(dynamics{i,j}(:,k));
            maxFoldChange{i,j}(1,k) = maxResponse{i,j}(1,k)/GFAPinit;
            adapt{i,j}(1,k) = abs(dynamics{i,j}(end,k) - GFAPinit)/GFAPinit;
      end
   end    
end