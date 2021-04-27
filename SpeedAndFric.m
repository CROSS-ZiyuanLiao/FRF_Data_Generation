fclose('all');
clc, clear;

n = 4;

leakPipeID = zeros(n, 1);
lFromStart = zeros(n, 1);
normalizedFRF = zeros(n, 640);
uRoughness = zeros(n, 12);
uWavespeed = zeros(n, 12);
% nonconvergent1 = cell(n, 1);
% nonconvergent2 = cell(n, 1);
% nonconvergent3 = cell(n, 1);

sigmaRough = [5e-3, 6.75e-3, 8.5e-3, 10.25e-3, 12e-3, 13.75e-3];
sigmaSpeed = [50, 70, 90, 110, 130, 150];

fprintf('Generating data. Please wait...\n');

timerVal = tic;

for i = 1: 6
    for j = 1: 6
        sigma1 = sigmaRough(i);
        sigma2 = sigmaSpeed(j);
        parfor k = 1: n
%             [nw, leakPipeID(k), lFromStart(k),...
%                 uRoughness(k, :), uWavespeed(k, :), nc]...
%                 = createNetwork(sigma1, sigma2);
            [nw, leakPipeID(k), lFromStart(k),...
                uRoughness(k, :), uWavespeed(k, :)]...
                = createNetwork(sigma1, sigma2);
            normalizedFRF(k, :) = getNormalizedFRF(nw);
%             nonconvergent1{k} = [nonconvergent1{k}; nc{1}];
%             nonconvergent2{k} = [nonconvergent2{k}; nc{2}];
%             nonconvergent3{k} = [nonconvergent3{k}; nc{3}];
        end
        fileName = strcat('./SpeedAndFric/R',num2str(i), 'S', num2str(j),'.mat');
        save(fileName,...
            'leakPipeID', 'lFromStart', 'normalizedFRF', 'uRoughness', 'uWavespeed');
    end
end

% save('./SpeedAndFric/NonConvergent.mat',...
%     'nonconvergent1', 'nonconvergent2', 'nonconvergent3');

elapsedTime = toc(timerVal);
fprintf('Data generated in %4.2f seconds.\n\n', elapsedTime);

clear;








% function [nw, leakPipeID, lFromStart, uRoughness, uWavespeed, nonconvergent]...
%     = createNetwork(sigmaRough, sigmaSpeed)
function [nw, leakPipeID, lFromStart, uRoughness, uWavespeed]...
    = createNetwork(sigmaRough, sigmaSpeed)

%------------------------------------basic settings---------------------------------------
% nonconvergent = {[], [], []};
t = 400;
dt = 1/125;

% *** perturbation settings ***
tauUnknown = [1,1.0001,1.0002,1.0005,1.0008,1.0017,1.0062,1.01154,1.0017];
tauUnknown = (tauUnknown - 1) * 1 + 1;
n = length(tauUnknown) - 1;
tauUnknown = [tauUnknown, ones(1, t/dt - n)];

while true
    %-----------------------------------------start---------------------------------------
    % *** create a new network ***
    nw = Network();
    
    %------------------------------------basic settings-----------------------------------
    % *** steady state solver parameter ***
    nw.steadyTolerance = 1e-12;
    nw.steadyMaxTrial = 400;
    
    % *** time settings ***
    nw.t = t;
    nw.dt = dt;
    
    %----------------------------------network parameters---------------------------------
    % *** incidence matrix ***
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
    
    % *** incidence matrix of fixed nodes ***
    nw.matrix.A10 = A(:, 1);
    
    % *** incidence matrix of unknown nodes ***
    nw.matrix.A12 = A(:, 2:8);
    
    % *** incidence matrix of leak nodes ***
    nw.matrix.A13 = A(:, 9);
    
    % *** friction mode ***
    nw.frictionMode = 'D-W2';
    
    % *** pipe properties ***
    nw.pipe.length = 10 * [100 60 80 60.42 25 55 40 55 40 50 70 150]';
    nw.pipe.diameter = 1e-3 * 10 * [40 40 32 25 25 40 40 40 40 40 32 40]';
    % ************ uncertainty of roughness *************
    while true
        roughness = [0.04 0.04 0.032 0.05 0.05 0.04 0.04 0.04 0.04 0.04 0.032 0.04]';
        uRoughness = sigmaRough * randn(length(roughness), 1);
        roughness = roughness + uRoughness;
        if min(roughness) > 0
            break;
        end
    end
    nw.pipe.roughness = roughness;
    % ************ uncertainty of wavespeed *************
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
    % *********************** end ***********************
    nw.pipe.setArea();
    
    % *** fixed node properties ***
    nw.nodeFixed.elevation = 0;
    nw.nodeFixed.initHead = 40;
    
    % *** unknown node properties ***
    nw.nodeUnknown.elevation = [0 0 0 0 0 0 0]';
    nw.nodeUnknown.initDemand = [0 0 0 0 0 0 0]';
    
    % *** leak node properties ***
    nw.nodeLeak.elevation = 0;
    nw.nodeLeak.pressure = 0;
    nw.nodeLeak.initDischargeCoeff = 1.73065e-2;
    
    % *** set numbers of network components ***
    nw.setnNode();
    nw.setnPipe();
    
    % *** randomly add 1 leak ***
    nPipe = nw.nPipe;
    leakPipeID = randi(nPipe);
    pipeLength = nw.pipe.length(leakPipeID);
    lFromStart = pipeLength * rand();
    coeff = 3 * 1.73065e-4;
    pressure = 0;
    
    nw.addOneLeak(leakPipeID, lFromStart, coeff, pressure);
    
    %-------------------------------------perturbataion-----------------------------------
    % *** basic settings ***
    nw.nodeFixed.transientTau = ones(nw.nNodeFixed, nw.t/nw.dt+1);
    nw.nodeUnknown.transientTau = ones(nw.nNodeUnknown, nw.t/nw.dt+1);
    nw.nodeLeak.transientTau = ones(nw.nNodeLeak, nw.t/nw.dt+1);
    
    % *** perturbation settings ***
    nw.nodeLeak.transientTau(end, :) = tauUnknown;
    
    %--------------------------------------run solvers------------------------------------
    % *** run solvers ***
    convergent = nw.steadySolver();
    %     convergent = randi(2) - 1;
    if convergent
        break;
%     else
%         nonconvergent{1} = [nonconvergent{1}; leakPipeID'];
%         nonconvergent{2} = [nonconvergent{2}; lFromStart'];
%         nonconvergent{3} = [nonconvergent{3}; uRoughness'];
    end
end
nw.adjustment();
nw.FRMSolver(640);
%------------------------------------------end--------------------------------------------
end


function normalizedFRF = getNormalizedFRF(nw)
%-----------------------------------------start-------------------------------------------
XFRM = nw.nodeLeak.FRMDemand(end, :);
YFRM = nw.nodeLeak.FRMHead(end, :);
FRF = YFRM./XFRM;
absFRF = abs(FRF);
normalizedFRF = absFRF/max(absFRF);

% *** test line, to be commented ***
% nw.MOCSolver();
% xMOC = nw.nodeLeak.transientDemand(end, 1: end-1)...
%     - nw.nodeLeak.initDemand(end);
% XMOC = fft(xMOC);
% yMOC = nw.nodeLeak.transientHead(end, 1: nw.nTimeGrid)...
%     - nw.nodeLeak.initHead(end);
% YMOC = fft(yMOC);
% FRFMOC = YMOC./XMOC;
% FRFMOC = FRFMOC(1: 640);
% absFRFMOC = abs(FRFMOC);
% nFRFMOC = absFRFMOC/max(absFRFMOC);
% omega = nw.omega(1: 640);
% plot(omega, normalizedFRF, '.', omega, nFRFMOC);
%------------------------------------------end--------------------------------------------
end

