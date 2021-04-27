% *** function DWCoeff ***

function f = DWCoeff(obj,Q,K,D)
%-----------------------------------------------------------------------------------------
% brief: Darcy-Weisbach friction factor using SI units
%
% output:
% f ---- column vector of Darcy-Weisbach friction factors
%
% input:
% Q ---- pipe flow
% K ---- pipe equivalent roughness
% D ---- pipe diameter
%-----------------------------------------------------------------------------------------
% column vector of Reynolds numbers
Re = 4* Q./D/pi/obj.VISCOS;

% Cubic Interpolation From Moody Diagram for 2000<Re<4000 (Dunlop, 1991):
e  = K./D; % column vector of relative roughness
AB = 5.74 / 4000^0.9;
Y2 = e/3.7 + AB;
Y3 = - 2 * log10(Y2);
FA = 1./Y3./Y3;
FB = FA .* (2 + 1.8*AB./Y2./log(Y2));
r = Re/2000;
X1 = 7*FA - FB;
X2 = 0.128 - 17*FA + 2.5*FB;
X3 = -0.128 + 13*FA - 2*FB;
X4 = r.*(0.032 - 3*FA + 0.5*FB);

% Calculating f using different equations depending on Re

% Hagen-Poiseuille formula for Re<2000 (Bhave, 1991)
% Swamee and Jain approximation to the Colebrook-White equation for Re>4000 (Bhave, 1991)
% Cubic Interpolation From Moody Diagram for 2000<Re<4000 (Dunlop, 1991)
f = 64./Re .* (Re <= 2000)...
    + 0.25 ./ (log10(e/3.7 + 5.74./Re.^0.9)).^2 .* (Re > 4000)...
    + (X1 + r.*(X2 + r.*(X3 + X4))) .* (Re>2000) .* (Re <= 4000);
%-----------------------------------------------------------------------------------------
end