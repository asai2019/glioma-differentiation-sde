function X = simulateCLE(numRealizations,dose,intrinsicNoise,extrinsicNoise,model)
%SIMULATECLE Simulates Chemical Langevin Equation model obtained from 
%   Sun et al. (2016) [BMC Systems Biology].
%
%   SIMULATECLE(numRealizations,dose,noise,model) produces dynamics 
%   pertaining to the glioma differentiation model. numRealizations 
%   indicates the number of cells/dynamics to expect. dose is the CT dose 
%   expected. intrinsic and extrinsic noise are the noise levels expected.
% 
K = numRealizations;
output = model.States.Outputs;
Nx = numel(output);
h = model.Time.Step;
tspan = model.Time.Duration(1):h:model.Time.Duration(2);
params = model.Parameters.Values;
X0 = model.States.X0;
X = zeros(Nx*numel(tspan),K);
% Iterate over individual cells
for k = 1:K
    y = X0;
    X(1,k) = y(output);
    step = 0;
    GWnoise=randn((1+max(tspan))*100,22);
    paramsk = params +randn(1,39)*extrinsicNoise;
    for p=h:h:max(tspan)
       step = step + 1;
       y = y + CLEdet(y,paramsk,dose)*h+sqrt(h)*intrinsicNoise*CLEsto(y,paramsk,dose,step,GWnoise);
       y(y<0) = 0;
       y(y>3) = 3;
       X(step+1,k) = real(y(output));
    end
end