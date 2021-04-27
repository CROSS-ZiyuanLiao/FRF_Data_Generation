% *** function steadySolver ***

function convergent = steadySolver(obj)
%-----------------------------------------------------------------------------------------
convergent = 1;

% friction mode
HW = 1;
DW1 = 2;
DW2 = 3;

% friction mode
% 'H-W'  ---- Hazen-Williams
% 'D-W1' ---- Darcy-Weisbach roughness
% 'D-W2' ---- Darcy-Weisbach f
modes = strcmpi(obj.frictionMode, {'H-W', 'D-W1', 'D-W2'});
frictionMode = find(modes==1);

% gravitational acceleration
g = obj.G;

% incidence matrix
A10 = obj.matrix.A10;
M12 = obj.matrix.A12;
M13 = obj.matrix.A13;
A12 = [M12, M13];
A21 = A12';

% number of pipes and nodes
np = obj.nPipe;
nu = obj.nNodeUnknown;
% nl = obj.nNodeLeak;

% properties of pipes
l = obj.pipe.length;
d = obj.pipe.diameter;
rough = obj.pipe.roughness;
Area = obj.pipe.area;

% properties of fixed node
H0 = obj.nodeFixed.initHead;

% properties of unknown nodes
Z2 = obj.nodeUnknown.elevation;
q2 = obj.nodeUnknown.initDemand;

% properties if leak nodes
e3 = obj.nodeLeak.elevation;
p3 = obj.nodeLeak.pressure;
Z3 = e3+p3;
Cv = obj.nodeLeak.initDischargeCoeff;

% calculation options
tol = obj.steadyTolerance;
N = obj.steadyMaxTrial;

% initial value before iteration
H2 = Z2;
H3 = 2*Z3;
q3 = Cv.*(Z3.^0.5);
Q = Area;
H = [H2; H3];
q = [q2; q3];

% pipe resistances and loss exponents
if frictionMode == HW
    R = 10.67*l./(rough.^1.852.*d.^4.87);
    nExp = 1.852*ones(np, 1);
elseif frictionMode == DW1
    f = obj.DWCoeff(Q, rough, d);
    R1 = (l./d)./(2*g.*Area.*Area);
    R = f.*R1;
    nExp = 2*ones(np, 1);
elseif frictionMode == DW2
    R = rough.*(l./d)./(2*g.*Area.*Area);
    nExp = 2*ones(np, 1);
end

% initial error
norm_err_energy = norm(R.*Q.*abs(Q).^(nExp - 1) + A12*H + A10*H0);
norm_err_mass = norm(A21*Q - q);
norm_err = norm([norm_err_energy, norm_err_mass]);

% error and H trace in iteration
errtrace = zeros(1, N+1);
Htrace = zeros(length(H), N+1);
errtrace(1) = norm_err;
Htrace(:, 1) = H;

% immutable matrix components
tol_A11 = 1e-5;                  % tolerence for matrix A11 to be invertible
invN = diag(1./nExp);
D22Diag2 = zeros(nu, 1);

% fprintf('Steady Solver running....\n');
% begin iterations
for i = 1: N
    % ************************************************************************************
    % matrix components
    A11Diag = R.*abs(Q).^(nExp-1);
    A11Diag(A11Diag < tol_A11) = tol_A11;      % small values in A11 make convergence unsteady therefore we need to define a lower bound, see Todini (1988), page 7.
    invA11 = diag(1./A11Diag);
    invD11 = diag(diag(invN).*diag(invA11));
    
    D22Diag3 = -0.5 * Cv .* (H3-Z3+eps).^(-0.5) .* (H3>Z3);
    D22Diag = [D22Diag2; D22Diag3];
    D22 = diag(D22Diag);
    % ************************************************************************************
    % relaxtion coeff = 1
    LHS = -(D22 - A21*invD11*A12);
    RHS = (A21*invN*(Q + invA11*A10*H0) + (q - A21*Q) + D22*H);
    dH = LHS\RHS + H;
    dQ = invN*Q + invD11*(A12*(H - dH) + A10*H0);
    Htemp = H - dH;
    Qtemp = Q - dQ;
    
    % update pipe resistances if Darcy-Weisbach formula used
    if frictionMode == DW1
        f = obj.DWCoeff(Qtemp+eps,rough,d);
        R = f.*R1;
    end
    
    % update q
    H3 = Htemp(nu+1: end);
    q3 = Cv .* sqrt(H3-Z3) .* (H3>Z3);
    q = [q2; q3];
    
    % update error
    norm_err_energy_1 = norm(R.*Qtemp.*abs(Qtemp).^(nExp - 1) + A12*Htemp + A10*H0);
    norm_err_mass_1 = norm(A21*Qtemp - q);
    norm_err_1 = norm([norm_err_energy_1, norm_err_mass_1]);
    
    % if error decrease, update and begin next iteration
    if norm_err_1 < norm_err
        H = Htemp;
        Q = Qtemp;
        norm_err = norm_err_1;
        %         norm_err_energy = norm_err_energy_1;
        
        errtrace(i+1) = norm_err;
        Htrace(:, i+1) = H;
        
        % convergence check
        if norm_err < tol
%             fprintf('Norm error of the energy conservation equations is %2.10f mH2O.\n', norm_err_energy_1);
%             fprintf('Norm error of the mass conservation equations is %2.10f m^3/s.\n', norm_err_mass_1);
%             fprintf('Final iteration: %d.\n', i);
            break;
        end
        % next iteration
        continue;
    end
    % ************************************************************************************
    % relaxtion coeff = 0.5
    Htemp = H - 0.5*dH;
    Qtemp = Q - 0.5*dQ;
    
    % update pipe resistances if Darcy-Weisbach formula used
    if frictionMode == DW1
        f = obj.DWCoeff(Qtemp,rough,d);
        R = f.*R1;
    end
    
    % update q
    H3 = Htemp(nu+1: end);
    q3 = Cv .* sqrt(H3-Z3) .* (H3>Z3);
    q = [q2; q3];
    
    % update error
    norm_err_energy_05 = norm(R.*Qtemp.*abs(Qtemp).^(nExp - 1) + A12*Htemp + A10*H0);
    norm_err_mass_05 = norm(A21*Qtemp - q);
    norm_err_05 = norm([norm_err_energy_05, norm_err_mass_05]);
    
    % if error decrease, update and begin next iteration
    if norm_err_05 < norm_err
        H = Htemp;
        Q = Qtemp;
        norm_err = norm_err_05;
        %         norm_err_energy = norm_err_energy_05;
        
        errtrace(i+1) = norm_err;
        Htrace(:, i+1) = H;
        
        % convergence check
        if norm_err < tol
%             fprintf('Norm error of the energy conservation equations is %2.10f mH2O.\n', norm_err_energy_05);
%             fprintf('Norm error of the mass conservation equations is %2.10f m^3/s.\n', norm_err_mass_05);
%             fprintf('Final iteration: %d.\n', i);
            break;
        end
        % next iteration
        continue;
    end
    % ************************************************************************************
    % polyfit to decide relaxation coeff
    x = [0, 0.5, 1];
    y = [norm_err, norm_err_05, norm_err_1];
    %     y = [norm_err_energy, norm_err_energy_05, norm_err_energy_1];
    p = polyfit(x, y, 2);
    rc = -0.5 * p(2)/p(1);
    
    % relaxtion coeff = rc
    Htemp = H - rc*dH;
    Qtemp = Q - rc*dQ;
    
    % update pipe resistances if Darcy-Weisbach formula used
    if frictionMode == DW1
        f = obj.DWCoeff(Qtemp,rough,d);
        R = f.*R1;
    end
    
    % update q
    H3 = Htemp(nu+1: end);
    q3 = Cv .* sqrt(H3-Z3) .* (H3>Z3);
    q = [q2; q3];
    
    % update error
    norm_err_energy_rc = norm(R.*Qtemp.*abs(Qtemp).^(nExp - 1) + A12*Htemp + A10*H0);
    norm_err_mass_rc = norm(A21*Qtemp - q);
    norm_err_rc = norm([norm_err_energy_rc, norm_err_mass_rc]);
    
    if rc > 0          % update and begin next iteration
        H = Htemp;
        Q = Qtemp;
        norm_err = norm_err_rc;
        %         norm_err_energy = norm_err_energy_rc;
        
        errtrace(i+1) = norm_err;
        Htrace(:, i+1) = H;
        
        % convergence check
        if norm_err < tol
%             fprintf('Norm error of the energy conservation equations is %2.10f mH2O.\n', norm_err_energy_rc);
%             fprintf('Norm error of the mass conservation equations is %2.10f m^3/s.\n', norm_err_mass_rc);
%             fprintf('Final iteration: %d.\n', i);
            break;
        end
        % next iteration
        continue;
    end
    % ************************************************************************************
    % relaxation coeff LEQ 0
    while rc <= 0
        x = [rc, 0, 0.5];
        y = [norm_err_rc, norm_err, norm_err_05];
        %         y = [norm_err_energy, norm_err_energy_05, norm_err_energy_1];
        p = polyfit(x, y, 2);
        rc = -0.5 * p(2)/p(1);
        
        %         fprintf('relaxation coeff is: %2.4f.\n', rc);
        
        Htemp = H - rc*dH;
        Qtemp = Q - rc*dQ;
        
        % update pipe resistances if Darcy-Weisbach formula used
        if frictionMode == DW1
            f = obj.DWCoeff(Qtemp,rough,d);
            R = f.*R1;
        end
        
        % update q
        H3 = Htemp(nu+1: end);
        q3 = Cv .* sqrt(H3-Z3) .* (H3>Z3);
        q = [q2; q3];
        
        % update error
        norm_err_energy_rc = norm(R.*Qtemp.*abs(Qtemp).^(nExp - 1) + A12*Htemp + A10*H0);
        norm_err_mass_rc = norm(A21*Qtemp - q);
        norm_err_rc = norm([norm_err_energy_rc, norm_err_mass_rc]);
        
        if norm_err_rc < norm_err
            break;
        end
    end
    
    %     fprintf('relaxation coeff is: %2.4f.\n', rc);
    
    % update and begin next iteration
    H = Htemp;
    Q = Qtemp;
    norm_err = norm_err_rc;
    %     norm_err_energy = norm_err_energy_rc;
    
    errtrace(i+1) = norm_err;
    Htrace(:, i+1) = H;
    
    % convergence check
    if norm_err < tol
%         fprintf('Norm error of the energy conservation equations is %2.10f mH2O.\n', norm_err_energy_rc);
%         fprintf('Norm error of the mass conservation equations is %2.10f m^3/s.\n', norm_err_mass_rc);
%         fprintf('Final iteration: %d.\n', i);
        break;
    end
    % ************************************************************************************
    
    %     LHS = D22 - A21*invD11*A12;
    %     RHS = A21*invN*(Q + invA11*A10*H0) + (q - A21*Q) + D22*H;
    %     H = LHS\RHS;
    %     Q = (speye(np,np) - invN)*Q - invD11*(A12*H + A10*H0);
    %     Q = (eye(np) - invN)*Q - invD11*(A12*H + A10*H0);
    
    % update pipe resistances if Darcy-Weisbach formula used
    %     if frictionMode == DW1
    %         f = obj.DWCoeff(Q,rough,d);
    %         R = f.*R1;
    %     end
    
    % update q
    %     H3 = H(nu+1: end);
    %     q3 = Cv .* sqrt(H3-Z3) .* (H3>Z3);
    %     q = [q2; q3];
    %
    % calculate 2-norm error of the energy equations, used to check convergence and print progress
    %     norm_err_energy = norm(R.*Q.*abs(Q).^(nExp - 1) + A12*H + A10*H0);
    
    % calculate 2-norm error of the mass equations, used to check convergence and print progress
    %     norm_err_mass = norm(A21*Q - q);
    
    % print progress
    %     fprintf('Iteration %d: Norm error of the energy conservation equations is %2.10f mH2O \n',i,norm_err_energy);
    %     fprintf('Iteration %d: Norm error of the mass conservation equations is %2.10f m^3/s \n',i,norm_err_mass);
    % ************************************************************************************
end

% adjustment for energy and mass convergence
% Q = round(full(Q), 8);
% Q = full(Q);

% H = round(H, 8);

% H3 = H(nu+1: end);

if norm_err < tol
    
    H2 = H(1: nu);
    q = A21 * Q;
    q2 = q(1: nu);
    q3 = q(nu+1: end);
    f = -(A12*H + A10*H0) .*d*2*g.*Area.*Area ./ (Q.*abs(Q).*l);
    
    obj.pipe.f = f;
    obj.pipe.initFlow = Q;
    
    obj.nodeUnknown.initHead = H2;
    obj.nodeLeak.initHead = H3;
    
    obj.nodeUnknown.initDemand = q2;
    obj.nodeLeak.initDemand = q3;
    
%     fprintf('Initial condition calculated.\n\n');
else
    convergent = 0;
    %     if frictionMode == DW1
    %         norm_err_energy = norm(f.*R1.*Q.*abs(Q).^(nExp - 1) + A12*H + A10*H0);
    %     else
    norm_err_energy = norm(R.*Q.*abs(Q).^(nExp - 1) + A12*H + A10*H0);
    %     end
    norm_err_mass = norm(A21*Q - q);
    
    fprintf('Convergence not reached after %d iterations.\n', N);
    fprintf('Norm error of the energy conservation equations is %2.10f mH2O.\n', norm_err_energy);
    fprintf('Norm error of the mass conservation equations is %2.10f m^3/s.\n', norm_err_mass);
end

% ******************** to be commented ********************
% yyaxis right;
% plot(1:i+1, errtrace(1:i+1));
% yyaxis left;
% plot(1:i+1, Htrace(end, 1:i+1));
% grid on;
% ******************** to be commented ********************
%-----------------------------------------------------------------------------------------
end