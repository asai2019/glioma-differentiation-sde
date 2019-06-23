function [C1,C2,C3] = Expt7clusterDynamics(dynamics)
% Produces 3 distinct clusters grouped by response dynamics
% and then sorted by mean response value
% C1 has highest mean response value, followed by C2 and C3
    [C1,C2,C3] = deal(cell(size(dynamics)));
    for i = 1:numel(dynamics)
        tmp = dynamics{i};
        idx = kmeans(zscore(tmp'),3,'replicates',100);
        % Compute mean response activation of each cluster
        clusterActivation = [mean(mean(tmp(:,idx == 1),2)),...
            mean(mean(tmp(:,idx == 2),2)),...
            mean(mean(tmp(:,idx == 3),2))];
        % Sort clusters by mean response activation
        [~,cind] = sort(clusterActivation,'descend');
        C1{i} = tmp(:,idx == cind(1));
        C2{i} = tmp(:,idx == cind(2));
        C3{i} = tmp(:,idx == cind(3));
    end
end