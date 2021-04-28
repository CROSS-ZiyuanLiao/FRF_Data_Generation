fclose('all');
clc, clear;

%% config
% number of the FRFs to be generated in each dataset
n = 6E4;

% standard deviation for roughness factor error
sigmaRough = [5e-3, 6.75e-3, 8.5e-3, 10.25e-3];

% standard deviation for wavespeed error
sigmaSpeed = [50, 70, 90, 110];

% leak scale multiplier (lowerbound and upperbound in each row)
scale = [2.0, 2.5];


%% generate data
% initialize
leakPipeID = zeros(n, 1);  % ID of the leak pipe
lFromStart = zeros(n, 1);  % leak position from pipe start
uRoughness = zeros(n, 12);  % roughness factor errors
uWavespeed = zeros(n, 12);  % wavespeed errors
multiplier = zeros(n, 1);  % leak scale multiplier
normalizedFRF = zeros(n, 640);  % FRF after min-max normalizaiton

fprintf('Generating data. Please wait...\n');

% start timing
timerVal = tic;

% generate data
for i = 1: length(sigmaRough)
    for j = 1: length(sigmaSpeed)
        % get error standard deviations
        sigma1 = sigmaRough(i);
        sigma2 = sigmaSpeed(j);
        
        parfor k = 1: n
            % TRF simulation
            [nw, leakPipeID(k), lFromStart(k), ...
                uRoughness(k, :), uWavespeed(k, :), multiplier(k)] = ...
                createNetwork(sigma1, sigma2, scale);
            
            % normalization
            normalizedFRF(k, :) = getNormalizedFRF(nw);
        end
        
        % save .mat file
        fileName = strcat(...
            './SpeedAndFric/R', num2str(i), 'S', num2str(j),'.mat');
        save(fileName,...
            'leakPipeID', 'lFromStart', 'normalizedFRF',...
            'uRoughness', 'uWavespeed', 'multiplier');
    end
end

% end timing
elapsedTime = toc(timerVal);

fprintf('Data generated in %4.2f seconds.\n\n', elapsedTime);
clear;


%% TRF simulation
function [nw, leakPipeID, lFromStart, uRoughness, uWavespeed, multiplier] = ...
    createNetwork(sigmaRough, sigmaSpeed, scale)
    % basic settings
    t = 400;  % total duration in time domain
    dt = 1/125;  % time interval

    % perturbation settings
    tauUnknown = [1,1.0001,1.0002,1.0005,1.0008,1.0017,1.0062,1.01154,1.0017];
    tauUnknown = (tauUnknown - 1) * 1 + 1;
    n = length(tauUnknown) - 1;
    tauUnknown = [tauUnknown, ones(1, t/dt - n)];

    while true
        % create a new network
        nw = Network();

        % network settings
        % # steady state solver parameter
        nw.steadyTolerance = 1e-12;
        nw.steadyMaxTrial = 400;

        % # time settings
        nw.t = t;
        nw.dt = dt;

        % # incidence matrix
        A = [
            -1  1  0  0  0  0  0  0  0;
            0  -1  1  0  0  0  0  0  0;
            0  -1  0  1  0  0  0  0  0;
            0  0  0  -1  1  0  0  0  0;
            0  0  0  0  -1  1  0  0  0;
            0  0  0  -1  0  1  0  0  0;
            0  0  -1  1  0  0  0  0  0;
            0  0  -1  0  0  0  1  0  0;
            0  0  0  0  0  -1  1  0  0;
            0  0  0  0  0  0  -1  1  0;
            0  0  0  0  0  -1  0  1  0;
            0  0  0  0  0  0  0  -1  1;
            ];

        % # incidence matrix of fixed nodes
        nw.matrix.A10 = A(:, 1);

        % # incidence matrix of unknown nodes
        nw.matrix.A12 = A(:, 2:8);

        % # incidence matrix of leak nodes
        nw.matrix.A13 = A(:, 9);

        % # friction mode
        nw.frictionMode = 'D-W2';

        % # pipe properties
        nw.pipe.length = 10 * [100 60 80 60.42 25 55 40 55 40 50 70 150]';
        nw.pipe.diameter = 1e-3 * 10 * [40 40 32 25 25 40 40 40 40 40 32 40]';

        % # uncertainty of roughness
        while true
            roughness = [0.04 0.04 0.032 0.05 0.05 0.04 0.04 0.04 0.04 0.04 0.032 0.04]';
            uRoughness = sigmaRough * randn(length(roughness), 1);
            roughness = roughness + uRoughness;
            if min(roughness) > 0
                break;
            end
        end
        nw.pipe.roughness = roughness;

        % # uncertainty of wavespeed
        while true
            wavespeed = [1431.13 1431.13 1451.53 1247.34 1247.34...
                1431.13 1431.13 1431.13 1431.13 1431.13 1451.53 1431.13]';
            uWavespeed = sigmaSpeed * randn(length(wavespeed), 1);
            wavespeed = wavespeed + uWavespeed;
            if min(wavespeed) > 0
                break;
            end
        end
        nw.pipe.wavespeed = wavespeed;

        % # pipe area
        nw.pipe.setArea();

        % # fixed node properties
        nw.nodeFixed.elevation = 0;
        nw.nodeFixed.initHead = 40;

        % # unknown node properties
        nw.nodeUnknown.elevation = [0 0 0 0 0 0 0]';
        nw.nodeUnknown.initDemand = [0 0 0 0 0 0 0]';

        % # leak node properties
        nw.nodeLeak.elevation = 0;
        nw.nodeLeak.pressure = 0;
        nw.nodeLeak.initDischargeCoeff = 1.73065e-2;

        % # set numbers of network components
        nw.setnNode();
        nw.setnPipe();

        % # randomly add one leak
        nPipe = nw.nPipe;
        leakPipeID = randi(nPipe);
        pipeLength = nw.pipe.length(leakPipeID);
        lFromStart = pipeLength * rand();
        multiplier = scale(1) + (scale(2) - scale(1)) * rand();
        coeff = multiplier * 1.73065e-4;
        pressure = 0;
        nw.addOneLeak(leakPipeID, lFromStart, coeff, pressure);

        % perturbataion
        % # basic settings
        nw.nodeFixed.transientTau = ones(nw.nNodeFixed, nw.t/nw.dt+1);
        nw.nodeUnknown.transientTau = ones(nw.nNodeUnknown, nw.t/nw.dt+1);
        nw.nodeLeak.transientTau = ones(nw.nNodeLeak, nw.t/nw.dt+1);

        % # valve perturbation settings
        nw.nodeLeak.transientTau(end, :) = tauUnknown;

        % run steady solver
        convergent = nw.steadySolver();
        if convergent
            break;
        end
    end

    % run TFR solver
    nw.adjustment();
    nw.FRMSolver(640);
end

%% min-max normalization
function normalizedFRF = getNormalizedFRF(nw)
    XFRM = nw.nodeLeak.FRMDemand(end, :);
    YFRM = nw.nodeLeak.FRMHead(end, :);
    FRF = YFRM./XFRM;
    absFRF = abs(FRF);
    normalizedFRF = absFRF/max(absFRF);
end
