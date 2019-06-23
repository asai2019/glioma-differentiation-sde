function X = simulateANM(numRealizations,dose,noise,model)
%SIMULATEANM Simulates Additive Noise Model obtained from Sun et al. (2016)
%   [BMC Systems Biology].
%
%   SIMULATEANM(numRealizations,dose,noise,model) produces dynamics pertaining
%   to the glioma differentiation model. numRealizations indicates the
%   number of cells/dynamics to expect. dose is the CT dose expected. 
%   noise is the noise level expected.
% 
% Please note that usage of sde_euler function requires SDETools Matlab
% package
K = numRealizations;
output = model.States.Outputs;
Nx = numel(output);
tspan = model.Time.Duration;
index = model.Parameters.Indices;
X0 = model.States.X0;
X = zeros(Nx*numel(tspan),K);
for k = 1:K
    % Set SDE solver settings
    opt = sdeset('CONSTGFUN','yes','SDEType','Ito');
    g = noise*ones(size(X0));
    % Simulate realization k of K
    x = sde_euler(model.Name,g,tspan,X0,opt,dose,index);
    % Extract selected outputs
    x = abs(x(:,output));
    % Collect individual cell trajectories
    X(:,k) = x;
end