% *** function addOneLeak ***

function addOneLeak(obj, pipeID, lFromStart, coeff, pressure, varargin)
%-----------------------------------------------------------------------------------------
% check pipe id
np = obj.nPipe;
if pipeID > np
    fprintf('The pipe ID exceeds the number of pipes.\n');
    fprintf('The network has %d pipes.\n', np);
    return;
end

% check the distance from the start node
if lFromStart <= 0
    fprintf('The distance from start node should not be less than zero.\n');
    return;
end

pipeLength = obj.pipe.length(pipeID);
if pipeLength <= lFromStart
    fprintf('The distance from start node exceeds length of the pipe.\n');
    fprintf('Pipe %d has a length of %2.2f meters.\n', np, pipeLength);
    return;
end

% check discharge coeff
if coeff < 0
    fprintf('The lump leak coeff should not be less than zero.\n');
    return;
end

% check external pressure
if pressure < 0
    fprintf('The pressure outside the leak should not be less than zero.\n');
    return;
end
%-----------------------------------------------------------------------------------------
% get network scale
nf = obj.nNodeFixed;
nu = obj.nNodeUnknown;
nl = obj.nNodeLeak;
nn = nf + nu + nl;

% arrange incidence matrix
M10 = obj.matrix.A10;
M12 = obj.matrix.A12;
M13 = obj.matrix.A13;
A = [M10, M12, M13];

% if lFromStart == 0
%     UnknownNodeID = find(A(pipeID, :)==-1);
%     obj.changeIntoLeak(UnknownNodeID, coeff, pressure);
%     return;
% end
% 
% if lFromStart == pipeLength
%     UnknownNodeID = find(A(pipeID, :)==1);
%     obj.changeIntoLeak(UnknownNodeID, coeff, pressure);
%     return;
% end

% change pipe properties
l = obj.pipe.length;
l(pipeID) = lFromStart;
l = [l; pipeLength-lFromStart];

d = obj.pipe.diameter;
d = [d; d(pipeID)];

rough = obj.pipe.roughness;
rough = [rough; rough(pipeID)];

a = obj.pipe.wavespeed;
a = [a; a(pipeID)];

% change leak node properties
Cv = obj.nodeLeak.initDischargeCoeff;
Cv = [coeff; Cv];

Z3 = obj.nodeLeak.pressure;
Z3 = [pressure; Z3];

p1 = A(pipeID, :);
iEnd = find(p1==1);
if nargin == 5
    e1 = obj.nodeFixed.elevation;
    e2 = obj.nodeUnknown.elevation;
    e3 = obj.nodeLeak.elevation;
    e = [e1; e2; e3];
    
    eStart = e(p1==-1);
    eEnd = e(iEnd);
    eLeak = eStart + (eEnd-eStart)*(lFromStart/pipeLength);
else
    eLeak = varargin{1};
end
elev = obj.nodeLeak.elevation;
elev = [eLeak; elev];

% change incidence matrix
node = zeros(np+1, 1);
node(pipeID) = 1;
node(np+1) = -1;

p2 = zeros(1, nn);
p2(iEnd) = 1;
p1(iEnd) = 0;

A(pipeID, :) = p1;
A = [A; p2];
A10 = A(:, nf);
A12 = A(:, (nf+1): (nf+nu));
A13 = A(:, (nf+nu+1): nn);
A13 = [node, A13];
%-----------------------------------------------------------------------------------------
% update the network
obj.matrix.A10 = A10;
obj.matrix.A12 = A12;
obj.matrix.A13 = A13;

obj.pipe.length = l;
obj.pipe.diameter = d;
obj.pipe.roughness = rough;
obj.pipe.wavespeed = a;

obj.nodeLeak.initDischargeCoeff = Cv;
obj.nodeLeak.pressure = Z3;
obj.nodeLeak.elevation = elev;

obj.setnNode();
obj.setnPipe();
obj.pipe.setArea();
%-----------------------------------------------------------------------------------------
end