% *** function adjustment ***

function adjustment(obj)
%-----------------------------------------------------------------------------------------
% gravitational acceleration
g = obj.G;

% incidence matrix
M10 = obj.matrix.A10;
M12 = obj.matrix.A12;
M13 = obj.matrix.A13;

% properties of pipes
l = obj.pipe.length;
d = obj.pipe.diameter;
Area = obj.pipe.area;
a = obj.pipe.wavespeed;

% initial conditions
QSteady = obj.pipe.initFlow;
H0Steady = obj.nodeFixed.initHead;
H2Steady = obj.nodeUnknown.initHead;
H3Steady = obj.nodeLeak.initHead;

% time configuration
t = obj.t;
dt = obj.dt;

% time and space grid set-up
nTimeGrid = round(t/dt, 0);
nPipeGrid = round(l./(a*dt), 0);
t = nTimeGrid * dt;

% l adjustment
l = nPipeGrid.*a*dt;

% f adjustment
% f = obj.pipe.roughness;
f = -([M12, M13]*[H2Steady; H3Steady] + M10*H0Steady)*2*g.*d.*Area.*Area...
    ./(QSteady.*abs(QSteady).*l);

% Record adjustment
obj.pipe.lengthAdjusted = l;
obj.pipe.fAdjusted = f;
obj.tAdjusted = t;
obj.nTimeGrid = nTimeGrid;
obj.nPipeGrid = nPipeGrid;

obj.T = (0: nTimeGrid)*dt;
obj.freq = (0: nTimeGrid-1)/t;
obj.omega = 2*pi() * ((0: nTimeGrid-1)/t);
%-----------------------------------------------------------------------------------------
end