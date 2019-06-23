function [CC, CCjack] = Expt4calculateCCdescriptor(descriptor)
%EXPT4CALCULATECCDESCRIPTOR Calculates the channel capacities for each
%summary descriptor obtained from Experiments 2 and 3
%
%  EXPT4CALCULATECCDESCRIPTOR(descriptor) computes the channel capacity for
%  a given dataset of descriptors across different values of k in the
%  k-nearest neighbors density estimation algorithm
%
warning('off','all')
cnt = 1;
[CC, CCjack] = deal(cell(1,1));
for k = [3:20]
    k
    nrepeats = 10;
    jackProb = linspace(0.6,0.95,20);
    njack = length(jackProb);
    %% Calculate information for each descriptor
    I = zeros(nrepeats,1);
    Ijack = zeros(nrepeats,njack);
    for repeat = 1:nrepeats
        tic
        [I(repeat),~,Ijack(repeat,:)] = calculateCC(descriptor,k);
        toc
    end
    I
    CC{1,cnt} = I;
    CCjack{1,cnt} = Ijack;
    cnt = cnt + 1;
end
