% -------------------------------------------------------------------------
% File: D4_SmallGain_condition.m
%
% Description:
%   Setting up gain matrices and evaluating spectral radii.
%
% Reference:
%   Christian Gatke, Julian D. Schiller, and Matthias A. Müller,
%   "Small-gain analysis of incremental input/output-to-state stability
%   for large-scale distributed systems", 2026 (not yet published).
%
% -------------------------------------------------------------------------


clear; clc;

% choose the file you want to check
LyapOpti = true;
if LyapOpti == true
    outer = load('LyapOpti_outer.mat','P_e','G_e','eta_e');
    inner  = load('LyapOpti_inner.mat', 'P_i','G_i','eta_i');
else
    outer = load('TrajectOpti_outer.mat','P_e','G_e','eta_e');
    inner  = load('TrajectOpti_inner.mat', 'P_i','G_i','eta_i');
end

Pe = (outer.P_e+outer.P_e')/2;
Ge = (outer.G_e+outer.G_e')/2;
Pi = (inner.P_i+inner.P_i')/2;
Gi = (inner.G_i+inner.G_i')/2;
eta_e = outer.eta_e;
eta_i = inner.eta_i;

M_grid = [3 5 10]; % number of subsystems
rho_L = zeros(numel(M_grid),1);
rho_T = zeros(numel(M_grid),1);

a = true; rhoMax_L = inf; rhoMax_T = inf;
for im = 1:numel(M_grid)
    M = M_grid(im);
    Gamma = zeros(M); G = zeros(M);
    etaVec = eta_i*ones(M,1);
    etaVec([1 M]) = eta_e;
    for i=1:M
        if i == 1 || i == M
            % both outer carriages are identic
            if i == 1, j = 2; else, j = M-1; end
            Gamma(i,j) = max(eig(Ge,Pi));
            g = computeGain(Ge,eye(size(Pe,2)));
            G(i,j) =sqrt(g/min(eig(Pe)));
        else
            % inner carriages: neighboring carriage can be an inner or
            % outer carriage
            jL = i-1; jR = i+1;
            if jL == 1, PL = Pe; else, PL = Pi; end
            if jR == M, PR = Pe; else, PR = Pi; end
            [Gamma(i,jL),Gamma(i,jR)] = computeGains(Gi,PL,PR);
            [gL,gR] = computeGains(Gi,eye(size(PL,2)),eye(size(PR,2)));
            G(i,jL) = sqrt(gL/min(eig(Pi)));
            G(i,jR) = sqrt(gR/min(eig(Pi)));
        end
    end
    % scale with corresponding eta
    Gamma = diag(1./(1-etaVec))*Gamma;
    G = diag(1./(1-sqrt(etaVec)))*G;

    % compute spectral radius
    rho_L(im) = max(abs(eig(Gamma)));
    rho_T(im) = max(abs(eig(G)));
    
    % compute maximum row sum
    if M >= 5 && a
        a = false;
        rhoMax_L = max(sum(Gamma,2));
        rhoMax_T = max(sum(G,2));
    end

    fprintf('\n========== M = %d ==========\n',M);
    fprintf('Lyapunov:  rho = %.1f\n', rho_L(im));
    fprintf('Trajectory:  rho = %.1f\n', rho_T(im));
end

fprintf('-----------------------------------\n')
if a
    fprintf('No spectral radius upper bound for M --> inf is computed. Choose M >= 5.\n')
else
    fprintf('Spectral radius upper bound for M --> inf:\n')
    fprintf('Lyapunov:  rho_max = %.2f\n', rhoMax_L);
    fprintf('Trajectory:  rho_max = %.2f\n', rhoMax_T);
end

% functions to compute gains
function [gL,gR] = computeGains(G,PL,PR)
    gL_var = sdpvar(1,1);
    gR_var = sdpvar(1,1);
    constr = [gL_var >= 0, gR_var >= 0];
    constr = [constr, ...
        blkdiag(gL_var*PL, gR_var*PR) >= G];
    optimize(constr, gL_var + gR_var, ...
        sdpsettings('solver','mosek','verbose',0));
    gL = value(gL_var);
    gR = value(gR_var);
    if gL < 0 || gR < 0, fprintf('Error: g < 0\n'); end
end

function g = computeGain(G,P)
    g_var = sdpvar(1,1);
    constr = g_var >= 0;
    constr = [constr, g_var*P >= G];
    optimize(constr, g_var, ...
        sdpsettings('solver','mosek','verbose',0));
    g = value(g_var);
    if g < 0, fprintf('Error: g < 0\n'); end
end