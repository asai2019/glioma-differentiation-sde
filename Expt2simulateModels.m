% Chooses experiment to perform in experiments 2 and 3
experiment = 'Experiment 2-CLE';
addpath(genpath('SDETools'));
switch experiment
    case 'Experiment 1'
        %% INITIALIZE MODEL SPECIFICATIONS
        % EXPERIMENT 1: Determine nature of relationship between differentiation
        % potential and conditional entropy of GFAP response given  
        % cholera toxin (CT) dose (10 ng/ml)
        % Uses 20 replications of 100 realizations for each noise level
        % Outputs: GFAP dynamics, which can be used to compute conditional
        % entropy of GFAP given CT = 10 ng/ml
        model.Name = @ANM; % Model function handle
        [x0,p0] = ANM; % Obtain initial conditions and nominal parameters
        model.States.X0 = x0; % Initial conditions
        model.States.Outputs = 10; % Outputs of interest (GFAP)
        model.Time.Duration = 0:0.1:48; % Timespan of interest (48 hrs)
        model.Parameters.Indices = 36; % Parameter of interest (CT)
        numRealizations = 100; % Number of realizations (i.e. cells) 
        dose = 10; % CT dose (ng/ml)
        noise = 0.001; % noise intensity (%)
        dynamics = simulateANM(numRealizations,dose,noise,model);
    case 'Experiment 2-ANM'
        %% INITIALIZE MODEL SPECIFICATIONS
        % Experiment 2A: Compute summary descriptors of GFAP dynamics given
        % ANM model. Use combinations of CT doses (0, 5, 7.5, 10 ng/ml) and
        % noise levels (0.1, 1, 5, 10%)
        % 500 realizations of each combination will be generated
        % Outputs: GFAP dynamics, which can be used to compute summary
        % descriptors (AUC, max fold change, change in adaptation, max GFAP
        % response, max GFAP response time)
        model.Name = @ANM; % Model function handle
        [x0,p0] = ANM; % Obtain initial conditions and nominal parameters
        model.States.X0 = x0; % Initial conditions
        model.States.Outputs = 10; % Outputs of interest (GFAP)
        model.Time.Duration = 0:0.01:48; % Timespan of interest (48 hrs)
        model.Parameters.Indices = 36; % Parameter of interest (CT)
        numRealizations = 500; % Number of realizations (i.e. cells)
        dose = 0; % CT dose (ng/ml)
        noise = 0.001; % noise intensity (%)
        dynamics = simulateANM(numRealizations,dose,noise,model);
        noise_levels = [0.001, 0.01, 0.05, 0.1];
        dose_levels = [0, 5, 7.5, 10];
        for i = 2:4
            for j = 1:4
                tmp = simulateANM(numRealizations,dose_levels(j),noise_levels(i),model);
                tmp = tmp(2:100:end,:);
                Expt2ANMdynamics{i,j} = tmp;
            end
        end
    case 'Experiment 2-CLE'
        %% INITIALIZE MODEL SPECIFICATIONS
        % Experiment 2B: Compute summary descriptors of GFAP dynamics given
        % CLE model. Use combinations of CT doses (0, 5, 7.5, 10 ng/ml) and
        % noise levels (LL, HL, LH, HH) intrinsic/extrinsic 0.001, 0.1
        % 500 realizations of each combination will be generated
        % Outputs: GFAP dynamics, which can be used to compute summary
        % descriptors (AUC< max fold change, change in adaptation, max GFAP
        % response, max GFAP response time)
        a = load('CLEparams.mat');
        model.Parameters.Values = a.params;
        model.States.X0 = a.x0;
        model.States.Outputs = 10;
        model.Time.Duration = [0 48];
        model.Time.Step = 0.01;
        numRealizations = 500;
        dose = 10;
        intrinsicNoise = 0.1;
        extrinsicNoise = 0.1;
        noiseLevels = [0.001 0.001; 0.1 0.001; 0.001 0.1; 0.1 0.1];
        doseLevels = [0 5 7.5 10];
        for i = 1:4
            for j = 1:4
                tmp = simulateCLE(numRealizations,doseLevels(j),noiseLevels(i,1),noiseLevels(i,2),model);
                Expt3CLEdynamics{i,j} = tmp(2:100:end,:);
            end 
        end
end

