% *** function changeIntoLeak ***

function changeIntoLeak(obj, UnknownNodeID, coeff, pressure)
%-----------------------------------------------------------------------------------------
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

% check node ID
nu = obj.nNodeUnknown;
if UnknownNodeID < 0 || UnknownNodeID > nu
    fprintf('Unknown node ID not found.\n');
    return;
end
if UnknownNodeID == 0
    return;
end
%-----------------------------------------------------------------------------------------
% arrange incidence matrix
M12 = obj.matrix.A12;
M13 = obj.matrix.A13;

% change leak node properties
Cv = obj.nodeLeak.initDischargeCoeff;
Cv = [coeff; Cv];

Z3 = obj.nodeLeak.pressure;
Z3 = [pressure; Z3];

eLeak = obj.nodeUnknown.elevation(UnknownNodeID);
elev = obj.nodeLeak.elevation;
elev = [eLeak; elev];

% change incidence matrix
node = M12(:, UnknownNodeID);
A13 = [node, M13];
A12 = M12;
A12(:, UnknownNodeID) = [];
%-----------------------------------------------------------------------------------------
% update the network
obj.matrix.A12 = A12;
obj.matrix.A13 = A13;

obj.nodeUnknown.elevation(UnknownNodeID) = [];
obj.nodeUnknown.initDemand(UnknownNodeID) = [];

obj.nodeLeak.initDischargeCoeff = Cv;
obj.nodeLeak.pressure = Z3;
obj.nodeLeak.elevation = elev;

obj.setnNode();
%-----------------------------------------------------------------------------------------
end