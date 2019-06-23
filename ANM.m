function [dx,p] = ANM(varargin)
%ANM Simulate ODEs for ANM model by Sun et al. 
% "Stochastic Modeling Suggests that Noise Reduces Differentiation
% Efficiency by Inducing a Heterogeneous Drug Response in Glioma
% Differentiation Therapy". BMC Systems Biology (2016).
%
%   Inputs:
%   t = time point at current step
%   x = model states at current time point
%   p = array of model parameters
%   I = indices of unknown parameters
%
%   Outputs:
%   dx = vector of parameters used to compute model states
%       OR
%        vector of initial conditions (if nargin == 0)
%       OR
%        names of all states (if nargin == 1)
%       
%    p = structure containing all parameters of the system
%       OR
%        names of all applicable inputs of the system (if nargin == 1)
%

% Extract input arguments
if(nargin >= 4)
    t = varargin{1};
    x = varargin{2};
    P = varargin{3};
    I = varargin{4};
end
%% Model Parameters
% Define model parameters
p.V1 = 1.6245; % 1
p.K1 = 9.5398; % 2
p.V2 = 0.6602; % 3
p.K2 = 0.5351; % 4
p.K3 = 0.7018; % 5
p.V4 = 0.1136; % 6
p.K4 = 0.4037; % 7
p.V5 = 0.4449; % 8
p.K5 = 0.2951; % 9
p.V6 = 0.5448; % 10
p.K6a = 0.9413; % 11
p.K6b = 0.8489; % 12
p.V7 = 0.639; % 13
p.K7 = 0.2324; % 14
p.V8 = 0.5359; % 15
p.K8 = 0.2558; % 16
p.V9 = 0.7861; % 17
p.K9 = 0.6108; % 18
p.V10ab = 0.0841; % 19
p.K10a = 0.1316; % 20
p.K10b = 0.1692; % 21
p.V10c = 0.5249; % 22
p.K10c = 0.7371; % 23
p.K10d = 0.8383; % 24
p.d1 = 0.4668; % 25 
p.d2 = 0.5774; % 26
p.d3 = 0.4602; % 27
p.d4 = 0.77; % 28
p.d5 = 0.7147; % 29
p.d6 = 0.8363; % 30
p.d7 = 0.7163; % 31
p.d8 = 0.9704; % 32
p.d9 = 0.7398; % 33
p.d10 = 0.4881; % 34
p.C = 1.25; % 35
p.CT = 10; % 36
p.n1 = 3; % 37
p.n2 = 10; % 38
p.n3 = 8; % 39
% Set initial conditions, define parameters and evaluate system of ODEs
if nargin == 0
    %% Set Initial Conditions
    % Set initial conditions if called with no current states
    x0 = zeros(10,1); % If not set below, intial condition = 0
    x0(1)=1; % PKA
    x0(2)=0.26; % CREB
    x0(3)=1; % P13K
    x0(4)=1; % AKT
    x0(5)=0.3; % GSK3B
    x0(6)=0.7239; % IL6
    x0(7)=0.408; % JAK2
    x0(8)=0.6818; % STAT3
    x0(9)=1.05; % Cyclin D1
    x0(10)=0.1834; % GFAP
    % Return initial conditions
    dx = x0;
else if nargin == 1
        %% Extract names of model states and inputs
        dx = {'PKA','CREB','P13K','AKT','pGSK3B','IL6','JAK2','STAT3','CyclinD1','GFAP'};
        p = {};      
        %% Evaluate System of ODEs
    else
        % Reassign parameters specified in I to new values specified in P
        if ~isempty(I) && ~isempty(P)
            names = fieldnames(p);
            for i = 1:length(I)
                ind = I(i);
                nam_temp = names{ind};
                eval(['p.',nam_temp,'=',num2str(P(i)),';'])
            end
        end
        % System of ODEs
        dx = zeros(size(x));% Initialize output vector
        dx(1,:) = p.d1 + (p.V1*p.CT^p.n1) / (p.K1^p.n1 + p.CT^p.n1) - p.d1*x(1,:); % PKA
        dx(2,:) = (p.V2*x(1,:)) ./ (p.K2 + x(1,:)) - p.d2*x(2,:); % CREB
        dx(3,:) = p.K3 ./ (p.K3 + x(1,:)) - p.d3*x(3,:); % P13K
        dx(4,:) = (p.V4 * x(3,:)) ./ (p.K4 + x(3,:)) - p.d4*x(4,:); % AKT
%         dx(5,:) = (p.V5 * (1-x(5,:)) .* x(4,:)) ./ (p.K5 + 1 - x(5,:)) - p.d5*x(5,:); % pGSK3B
        dx(5,:) = (p.V5 * x(4,:)) ./ (p.K5 + x(4,:)) - p.d5*x(5,:); % pGSK3B
        dx(6,:) = (p.V7 * x(1,:)) ./ (p.K7 + x(1,:)) - p.d7*x(6,:); % IL6
        dx(7,:) = (p.V8 * x(6,:)) ./ (p.K8 + x(6,:)) - p.d8*x(7,:); % JAK2
        dx(8,:) = (p.V9 * x(7,:)) ./ (p.K9^p.n1 + x(7,:)) - p.d9*x(8,:); % STAT3
        dx(9,:) = (p.V6 * x(9,:).^p.n2) ./ (p.K6a^p.n2 + x(9,:).^p.n2) - (p.d6 * (1-x(5,:)) .* x(9,:)) ./ (p.K6b + (1 - x(5,:))); % CyclinD1
%         dx(10,:) = max((p.C - x(9,:)) ./ p.C, 0) .* ((p.V10ab*x(2,:).*x(8,:)) ./ ((p.K10a+x(2,:)) .* (p.K10b+x(8,:))) + (p.V10c * (1 - x(5,:)).^p.n3) ./ (p.K10c + (1 - x(5,:)).^p.n3)) - p.d10*x(10,:); % GFAP
        dx(10,:) = p.V10ab*x(2,:)./(p.K10a+x(2,:)).*((x(8,:)./(p.K10b+x(8,:)))+(x(9,:)./(p.K10d+x(9,:))))+p.V10c*(((p.C-x(9,:))./p.C.*(1-x(5,:))).^p.n3)./(p.K10c^p.n3+(1-x(5,:)).^p.n3)-p.d10*x(10,:);
    end
end
end




