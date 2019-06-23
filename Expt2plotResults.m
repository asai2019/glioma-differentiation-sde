function Expt2plotResults(descriptor,metric)
%EXPT2PLOTRESULTS Plots heatmaps corresponding to target metric from
%  Experiment 2
%
%  EXPT2PLOTRESULTS(descriptor,metric) plots results corresponding to each
%  dose and noise level for a specific summary descriptor calculated from
%  data in Experiment 2
%
meanDescriptor = cellfun(@mean,descriptor);
figure;
colormap('jet');
imagesc(meanDescriptor);
colorbar;
xticks([1 2 3 4]);
xticklabels({'0','5','7.5','10'});
xlabel('CT Dose (ng/ml)');
yticks([1 2 3 4]);
yticklabels({'0.1','1','5','10'});
ylabel('Noise Intensity (%)');
title(metric);