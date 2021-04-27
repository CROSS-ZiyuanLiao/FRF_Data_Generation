% *** function MOCSolver ***

function MOCSolver(obj)
%-----------------------------------------------------------------------------------------
% gravitational acceleration
g = obj.G;

% incidence matrix
M10 = obj.matrix.A10;
M10u = (M10==-1);
M10d = (M10==1);

M12 = obj.matrix.A12;
M12u = (M12==-1);
M12d = (M12==1);
% M21 = M12';
M13 = obj.matrix.A13;
M13u = (M13==-1);
M13d = (M13==1);

% A12 = [M12, M13];
% A21 = A12';

% number of pipes and nodes
np = obj.nPipe;
% nf = obj.nNodeFixed;

nu = obj.nNodeUnknown;
uExist = 0;
if nu > 0
    uExist = 1;
end

nl = obj.nNodeLeak;
lExist = 0;
if nl > 0
    lExist = 1;
end

% properties of pipes
% l = obj.pipe.length;
d = obj.pipe.diameter;
Area = obj.pipe.area;
a = obj.pipe.wavespeed;

% friction mode
% HW = 1;
% DW = 2;
% modes = strcmpi(obj.frictionMode, {'H-W', 'D-W'});
% frictionMode = find(modes==1);
% if frictionMode == HW
%     nExp = 1.852*ones(np, 1);  
% else
%     nExp = 2*ones(np, 1);
% end

% initial conditions
QSteady = obj.pipe.initFlow;
% f = obj.pipe.f;
H0Steady = obj.nodeFixed.initHead;
H2Steady = obj.nodeUnknown.initHead;
q2Steady = obj.nodeUnknown.initDemand;
H3Steady = obj.nodeLeak.initHead;
q3Steady = obj.nodeLeak.initDemand;
CvSteady = obj.nodeLeak.initDischargeCoeff;
e3 = obj.nodeLeak.elevation;
p3 = obj.nodeLeak.pressure;
Z3 = e3+p3;

% time configuration
% t = obj.t;
dt = obj.dt;

% perturbation configuration
tauH0 = obj.nodeFixed.transientTau;
tauq2 = obj.nodeUnknown.transientTau;
tauCv = obj.nodeLeak.transientTau;

H0 = H0Steady.*tauH0;

q2 = [];
if uExist
    q2 = q2Steady.*tauq2;
end

Cv = [];
if lExist
    Cv = CvSteady.*tauCv;
end

% time and space grid set-up
nTimeGrid = obj.nTimeGrid;
nPipeGrid = obj.nPipeGrid;
% t = nTimeGrid * dt;

% l adjustment
% l = nPipeGrid.*a*dt;

% f adjustment
% f = -([M12, M13]*[H2Steady; H3Steady] + M10*H0Steady)*2*g.*d.*Area.*Area...
%     ./(QSteady.*abs(QSteady).*l);
f = obj.pipe.fAdjusted;

% initial flow in each pipe
Q = cell(np, 1);
for i = 1: np
    Q{i} = zeros(nPipeGrid(i)+1, nTimeGrid+1);
    Q{i}(:, 1) = QSteady(i);
end

% initial head in each pipe
H = cell(np, 1);
HStart = ([M10, M12, M13]==-1) * [H0Steady; H2Steady; H3Steady];
HEnd = ([M10, M12, M13]==1) * [H0Steady; H2Steady; H3Steady];
for i = 1: np
    H{i} = zeros(nPipeGrid(i)+1, nTimeGrid+1);
    H{i}(:, 1) = linspace(HStart(i), HEnd(i), nPipeGrid(i)+1);
end

% initial unknown node head
H2 = [];
if uExist    
    H2 = zeros(nu, nTimeGrid+1);
    H2(:, 1) = H2Steady;
end

% initial leak node head and demand
H3 = [];
q3 = [];
if lExist
    H3 = zeros(nl, nTimeGrid+1);
    q3 = zeros(nl, nTimeGrid+1);
    H3(:, 1) = H3Steady;
    q3(:, 1) = q3Steady;
end

fprintf('MOC running....\n');
tic;
% method of characteristic
R = f./(2*d.*Area);
Ca = g*Area./a;
Cp = cell(np, 1);
Cn = cell(np, 1);

for j = 2: (nTimeGrid+1)
    % update Cp & Cn
    for i = 1: np
        Cp{i} = Q{i}(:, j-1) + Ca(i)*H{i}(:, j-1)...
            -R(i)*dt*Q{i}(:, j-1).*abs(Q{i}(:, j-1));
        Cn{i} = Q{i}(:, j-1) - Ca(i)*H{i}(:, j-1)...
            -R(i)*dt*Q{i}(:, j-1).*abs(Q{i}(:, j-1));
    end
    
    % interior grids
    for i = 1: np
        tempCp = Cp{i}(1: nPipeGrid(i)-1);
        tempCn = Cn{i}(3: nPipeGrid(i)+1);
        H{i}(2:nPipeGrid(i), j) = (tempCp-tempCn)/2/Ca(i);
        Q{i}(2:nPipeGrid(i), j) = (tempCp+tempCn)/2;
    end
    
    % boundaries, i.e. nodes
    tempCp = zeros(np, 1);
    tempCn = zeros(np, 1);
    for i = 1: np
        tempCp(i) = Cp{i}(end-1);
        tempCn(i) = Cn{i}(2);
    end
    
    % fixed nodes
    tempH0 = H0(:, j);
    tempQ0 = tempCp.*M10d + tempCn.*M10u - Ca.*M10.*tempH0';
    
    % unknown nodes
    if uExist
        tempH2 = (sum(M12.*(tempCp.*M12d+tempCn.*M12u),1) - q2(:, j)')...
            ./sum(Ca.*M12.*M12,1);
        tempQ2 = tempCp.*M12d + tempCn.*M12u - Ca.*M12.*tempH2;
        tempH2 = tempH2';
        H2(:, j) = tempH2;
    end
    
    % leak nodes   
    if lExist
        aa = sum(M13.*(tempCp.*M13d+tempCn.*M13u),1);
        bb = sum(Ca.*M13.*M13,1);
        tempH3 = H3(:, j-1)';
        tempCv = Cv(:, j)';
        for N = 1: 40
            tempH3 = tempH3 - ...
                (aa - bb.*tempH3 - tempCv.*sqrt(tempH3-Z3').*(tempH3>Z3'))...
                ./(-bb - 0.5*tempCv.*(tempH3-Z3'+eps).^(-0.5).*(tempH3>Z3'));
            err = aa - bb.*tempH3 - tempCv.*sqrt(tempH3-Z3').*(tempH3>Z3');
            if err < 1e-18
                break;
            end
        end
        tempQ3 = tempCp.*M13d + tempCn.*M13u - Ca.*M13.*tempH3;        
        tempq3 = sum(M13.*(tempCp.*M13d+tempCn.*M13u-Ca.*M13.*tempH3),1);
        tempH3 = tempH3';
        tempq3 = tempq3';
        H3(:, j) = tempH3;
        q3(:, j) = tempq3;
    end
    
    % give node data to connect pipes
    if uExist && ~lExist
        Qu = sum(tempQ0.*M10u,2) + sum(tempQ2.*M12u,2);
        Qd = sum(tempQ0.*M10d,2) + sum(tempQ2.*M12d,2);
        Hu = sum(tempH0'.*M10u,2) + sum(tempH2'.*M12u,2);
        Hd = sum(tempH0'.*M10d,2) + sum(tempH2'.*M12d,2);
    elseif ~uExist && lExist
        Qu = sum(tempQ0.*M10u,2) + sum(tempQ3.*M13u,2);
        Qd = sum(tempQ0.*M10d,2) + sum(tempQ3.*M13d,2);
        Hu = sum(tempH0'.*M10u,2) + sum(tempH3'.*M13u,2);
        Hd = sum(tempH0'.*M10d,2) + sum(tempH3'.*M13d,2);
    elseif uExist && lExist
        Qu = sum(tempQ0.*M10u,2) + sum(tempQ2.*M12u,2) + sum(tempQ3.*M13u,2);
        Qd = sum(tempQ0.*M10d,2) + sum(tempQ2.*M12d,2) + sum(tempQ3.*M13d,2);
        Hu = sum(tempH0'.*M10u,2) + sum(tempH2'.*M12u,2) + sum(tempH3'.*M13u,2);
        Hd = sum(tempH0'.*M10d,2) + sum(tempH2'.*M12d,2) + sum(tempH3'.*M13d,2);
    end
    
    for i = 1: np
        Q{i}(1, j) = Qu(i);
        Q{i}(nPipeGrid(i)+1, j) = Qd(i);
        H{i}(1, j) = Hu(i);
        H{i}(nPipeGrid(i)+1, j) = Hd(i);
    end
end

% obj.tAdjusted = t;

% obj.pipe.f = f;
% obj.pipe.lengthAdjusted = l;
obj.pipe.transientFlow = Q;
obj.pipe.transientHead = H;

obj.nodeFixed.transientHead = H0;

obj.nodeUnknown.transientDemand = q2;
obj.nodeUnknown.transientHead = H2;

obj.nodeLeak.transientDemand = q3;
obj.nodeLeak.transientHead = H3;

toc;
fprintf('MOC done.\n\n');
%-----------------------------------------------------------------------------------------
end