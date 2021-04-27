% *** function FRMSolver ***

function FRMSolver(obj, n, adjust)
%-----------------------------------------------------------------------------------------
if nargin < 3
    adjust = 0;
end

if nargin < 2
    n = [];
end

% gravitational acceleration
g = obj.G;

% incidence matrix
M10 = obj.matrix.A10;
M10u = (M10==-1);
M10d = (M10==1);

M12 = obj.matrix.A12;
M12u = (M12==-1);
M12d = (M12==1);

M13 = obj.matrix.A13;
M13u = (M13==-1);
M13d = (M13==1);

% number of pipes and nodes
np = obj.nPipe;

nu = obj.nNodeUnknown;
uExist = 0;
if nu > 0
    uExist = 1;
end

% nf = obj.nNodeFixed;

nl = obj.nNodeLeak;
lExist = 0;
if nl > 0
    lExist = 1;
end

% properties of pipes
d = obj.pipe.diameter;
Area = obj.pipe.area;
a = obj.pipe.wavespeed;
QSteady = obj.pipe.initFlow;

if adjust == 1
    l = obj.pipe.lengthAdjusted;
    f = obj.pipe.fAdjusted;
else
    l = obj.pipe.length;
    f = obj.pipe.f;
end

% initial conditions
q3Steady = obj.nodeLeak.initDemand;
H3Steady = obj.nodeLeak.initHead;
e3 = obj.nodeLeak.elevation;
p3 = obj.nodeLeak.pressure;
Z3 = e3+p3;
q2Steady = obj.nodeUnknown.initDemand;
H0Steady = obj.nodeFixed.initHead;

% time and freq grid
nTimeGrid = obj.nTimeGrid;
omega = obj.omega;

% tau of unknown nodes
tauUnknown = obj.nodeUnknown.transientTau;
q2Hat = zeros(1, nTimeGrid);
if uExist
    q2Aster = q2Steady.*(tauUnknown(:, 1:nTimeGrid) - 1);
    q2Hat = fft(q2Aster, nTimeGrid, 2);
end

% tau of fixed nodes
tauFixed = obj.nodeFixed.transientTau;
H0Aster = H0Steady.*(tauFixed(:, 1:nTimeGrid) - 1);
H0Hat = fft(H0Aster, nTimeGrid, 2);

% tau of leak nodes
tauLeak = obj.nodeLeak.transientTau;
tauLeakHat = zeros(1, nTimeGrid);
if lExist
    tauLeakAster = tauLeak(:, 1:nTimeGrid) - 1;
    tauLeakHat = fft(tauLeakAster, nTimeGrid, 2);
end

if ~isempty(n)
    nTimeGrid = n;
    q2Hat = q2Hat(:, 1: nTimeGrid);
    H0Hat = H0Hat(:, 1: nTimeGrid);
    tauLeakHat = tauLeakHat(:, 1: nTimeGrid);
    omega = omega(:, 1: nTimeGrid);
end

% fprintf('FRM running....\n');
% tic;
% friction component
Rs = f .* abs(QSteady) ./ (1i*omega.*d.*Area + eps); % steady friciton component
Ru = 0; % unsteady friciton component
Rv = 0; % viscoelastic component

% propagation constant gamma for all pipes and omegas
gamma = 1i * omega ./ a .* sqrt((1 + Rs + Ru) * (1 + Rv));

% ************ step 1 ************
% characteristic impedance Z for all pipes and omegas, where row i represents pipe i
ZAll = a ./ (g * Area) .* sqrt((1 + Rs + Ru) * (1 + Rv)^(-1));
SAll = sinh(gamma .* l);
CAll = cosh(gamma .* l);

% solution vector
tempResult = zeros(np*2+nu+nl, nTimeGrid);

% construct multiple leak formulation
[LHS11,LHS12,LHS13,LHS14] = deal([]);
if uExist
    LHS11 = M12d';
    LHS12 = -M12u';
    LHS13 = zeros(nu);
    LHS14 = zeros(nu,nl);
end
[LHS21,LHS22,LHS23,LHS24] = deal([]);
if lExist
    LHS21 = M13d';
    LHS22 = -M13u';
    LHS23 = zeros(nl,nu);
    LHS24 = -diag(0.5*q3Steady./(H3Steady-Z3));
end
LHS31 = eye(np);
LHS41 = zeros(np);

if uExist && lExist         % unknowns and leaks
    exist = 1;
elseif uExist && (~lExist)  % unknowns but no leaks
    exist = 2;
elseif (~uExist) && lExist  % no unknowns but leaks
    exist = 3;
end

% solve equation for each omega
for j = 1: nTimeGrid
    
    Z = ZAll(:, j);
    invZ = diag(1./Z);
    Z = diag(Z);
    S = diag(SAll(:, j));
    C = diag(CAll(:, j));
    
    if exist == 1
        q2 = q2Hat(:, j);
        tau = tauLeakHat(:, j);
        r = H0Hat(:, j);
        
        LHS32 = -C;
        LHS33 = invZ*S*M12u;
        LHS34 = invZ*S*M13u;
        LHS42 = Z*S;
        LHS43 = M12d-C*M12u;
        LHS44 = M13d-C*M13u;
        
        RHS1 = q2;
        RHS2 = diag(q3Steady)*tau;
        RHS3 = invZ*S*M10u*r;
        RHS4 = (C*M10u-M10d)*r;
        
        % form matrix
        LHS1 = [LHS11, LHS12, LHS13, LHS14];
        LHS2 = [LHS21, LHS22, LHS23, LHS24];
        LHS3 = [LHS31, LHS32, LHS33, LHS34];
        LHS4 = [LHS41, LHS42, LHS43, LHS44];
        
        LHS = [LHS1; LHS2; LHS3; LHS4];
        RHS = [RHS1; RHS2; RHS3; RHS4];
        
    elseif exist == 2
        q2 = q2Hat(:, j);
        r = H0Hat(:, j);
        
        LHS32 = -C;
        LHS33 = invZ*S*M12u;
        LHS42 = Z*S;
        LHS43 = M12d-C*M12u;
        
        RHS1 = q2;
        RHS3 = invZ*S*M10u*r;
        RHS4 = (C*M10u-M10d)*r;
        
        % form matrix
        LHS1 = [LHS11, LHS12, LHS13];
        LHS3 = [LHS31, LHS32, LHS33];
        LHS4 = [LHS41, LHS42, LHS43];
        
        LHS = [LHS1; LHS3; LHS4];
        RHS = [RHS1; RHS3; RHS4];
        
    elseif exist == 3
        tau = tauLeakHat(:, j);
        r = H0Hat(:, j);
        
        LHS32 = -C;
        LHS34 = invZ*S*M13u;
        LHS42 = Z*S;
        LHS44 = M13d-C*M13u;
        
        RHS2 = diag(q3Steady)*tau;
        RHS3 = invZ*S*M10u*r;
        RHS4 = (C*M10u-M10d)*r;
        
        % form matrix
        LHS2 = [LHS21, LHS22, LHS24];
        LHS3 = [LHS31, LHS32, LHS34];
        LHS4 = [LHS41, LHS42, LHS44];
        
        LHS = [LHS2; LHS3; LHS4];
        RHS = [RHS2; RHS3; RHS4];
    end
    
    x = pinv(LHS)*RHS;
    tempResult(:, j) = x;
end

% ************ step 2 ************
obj.nodeFixed.FRMHead = H0Hat;

qD = tempResult(1:np, :);
qU = tempResult((np+1):(2*np), :);

obj.pipe.FRMFlowDown = qD;
obj.pipe.FRMFlowUp = qU;

if exist == 1
    obj.nodeUnknown.FRMDemand = q2Hat;
    obj.nodeUnknown.FRMHead = tempResult((2*np+1):(2*np+nu), :);
    
    obj.nodeLeak.FRMCoeff = tauLeakHat;
    obj.nodeLeak.FRMHead = tempResult((2*np+nu+1):(2*np+nu+nl), :);
    obj.nodeLeak.FRMDemand = M13d'*qD - M13u'*qU;
    
%     obj.pipe.FRMFlowDown = tempResult(1:np, :);
%     obj.pipe.FRMFlowUp = tempResult((np+1):(2*np), :);
    
elseif exist == 2
    obj.nodeUnknown.FRMDemand = q2Hat;
    obj.nodeUnknown.FRMHead = tempResult((2*np+1):(2*np+nu), :);
    
%     obj.pipe.FRMFlowDown = tempResult(1:np, :);
%     obj.pipe.FRMFlowUp = tempResult((np+1):(2*np), :);
    
elseif exist == 3
    obj.nodeLeak.FRMCoeff = tauLeakHat;
    obj.nodeLeak.FRMHead = tempResult((2*np+1):(2*np+nl), :);
    obj.nodeLeak.FRMDemand = M13d'*qD - M13u'*qU;
    
%     obj.pipe.FRMFlowDown = tempResult(1:np, :);
%     obj.pipe.FRMFlowUp = tempResult((np+1):(2*np), :);    
end

% toc;
% fprintf('FRM done.\n\n');
%-----------------------------------------------------------------------------------------
end