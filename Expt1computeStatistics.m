function [diffPotential,entropy] = Expt1computeStatistics(dynamics)
%EXPT1COMPUTESTATISTICS Computes relevant statistics to report from 
%  Experiment 1, in which entropy and differentiation potential are found
%  for various noise levels given CT dose of 10 ng/ml to determine if there
%  is a relationship between the two quantities
%
%  EXPT1COMPUTESTATISTICS(dynamics) computes the entropy and
%  differentiation potential for the cell dynamics produced across noise
%  levels and replications
%
diffPotential  = cell(1,size(dynamics,2));
entropy = cell(1,size(dynamics,2));
bins = -0.05:0.1:1.25;
for j = 1:size(dynamics,2)
    for i = 1:size(dynamics,1)
        for k = 1:size(dynamics{i,j},1)
            tmp = dynamics{i,j}(k,:);
            H = histcounts(tmp,bins,'Normalization','probability');
            J = -H.*log2(H);
            entropy{j}(k,i) = sum(J(~isnan(J)));
            diffPotential{j}(k,i) = sum(tmp > 0.8)/length(tmp);
        end
    end
end