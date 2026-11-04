% -------------------------------------------------------------------------
% File: D2_2_LMI_SDP_inner.m
%
% Description:
%   LMI verification of the intermediate carriages.
%
% Reference:
%   Christian Gatke, Julian D. Schiller, and Matthias A. Müller,
%   "Small-gain analysis of incremental input/output-to-state stability
%   for large-scale distributed systems", 2026 (not yet published).
%
% -------------------------------------------------------------------------


clear; clc; yalmip('clear');
% select optimization strategy:
% p.LyapOpti = true; --> minimizes Lyapunov gains
% p.LyapOpti = false; --> minimizes trajectory-based gains
p.LyapOpti = true;

%% parameters and setting
p.dt = 0.01; p.m = 90e3; p.k = 1e7; p.d = 5e5;
p.bound = 0.8; % bound for relative velocity between carriages
etaGrid = [0 0.1 0.15 0.2 0.22 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.85 0.90 0.95 0.98 0.99];
betaStart = 0.5; betaMax = 1e4;
nBisect = 16; betaTol = 1e-3;

%% Bisection: finding the smallest beta for each eta
best.beta = inf;
for ie = 1:numel(etaGrid)
    eta = etaGrid(ie);
    fprintf('\neta = %.3f\n',eta);
    high = betaStart;
    cert = [];

    % find a feasible solution
    while high <= betaMax
        [status,sol] = solveSDP(eta,high,p);
        if status == 1
            fprintf('  Solution found for eta = %.3f and beta = %.3f.\n', eta, high);
            cert = sol;
            break;
        end
        high = 2*high;
    end

    if isempty(cert)
        fprintf('  no beta found.\n');
        continue
    end

    % find smallest upper bound
    fprintf('  Start bisection.\n')
    low = 0;
    for k = 1:nBisect
        if high-low <= betaTol, break; end
        test = (low+high)/2;
        [status,sol] = solveSDP(eta,test,p);
        if status == 1
            high = test; cert = sol;
        elseif status == 0
            low = test;
        else
            % abort bisection; we keep the last feasible upper bound
            break
        end
    end
    fprintf('  feasible beta = %.6g\n',high);
    if high < best.beta
        best = cert;
        best.beta = high;
        best.eta = eta;
    end
end

if ~isfinite(best.beta)
    error('no local LMI solved: nothing is saved.');
end

%% save data
P_i = best.P;
G_i = best.G;
Q_i = best.Q;
R_i = best.R;
eta_i = best.eta;

if p.LyapOpti == true
    save('LyapOpti_inner.mat','P_i','G_i','eta_i', 'Q_i', 'R_i');
else
    save('TrajectOpti_inner.mat','P_i','G_i','eta_i', 'Q_i', 'R_i');
end
fprintf('eta_i = %.6f, beta = %.6g\n',eta_i,best.beta);
fprintf('P_i =\n'); disp(P_i);
fprintf('G_i =\n'); disp(G_i);


%% function solveSDP
function [status,sol] = solveSDP(eta,beta,p)
    yalmip('clear');
    status = -1; sol = [];

    % decision variables
    P = sdpvar(2,2,'symmetric');
    Q = sdpvar(3,3,'diagonal');
    R = sdpvar(1,1);
    G = sdpvar(4,4,'symmetric');
    c = sdpvar(1,1);
    
    % Jacobian matrices
    B = [1 0 0; 0 1 0]; C = [1 0]; D = [0 0 1];

    constr = [trace(P)==1, P>=1e-6*eye(2), G>=0, Q>=0, R>=0];
    
    % optimization strategy
    if p.LyapOpti == true
         constr = [constr, G<=beta*(1-eta)/2*blkdiag(P,P)];
    else
        constr = [constr, c >= 1e-10, P >= c*eye(size(P,1))];
        constr = [constr, G <= (beta^2*(1-sqrt(eta))^2/4)*c*eye(size(G,1))];
    end

    % cubic damping: h = (max. relative velocity)^2
    % coupled dynamics: same h-value in A and E
    h = p.bound^2;
    vertices = [0 0; 0 h; h 0; h h];

    a = p.dt*p.k/p.m;
    b = 3*p.dt*p.d/p.m;
    for k = 1:size(vertices,1)
        A = [1 p.dt; -2*a 1-b*(vertices(k,1)+vertices(k,2))];
        E = [0 0 0 0; a b*vertices(k,1) a b*vertices(k,2)];
        
        % LMI conditions
        J = [A B E];
        Y = [C D zeros(1,4)];
        L = J'*P*J - blkdiag(eta*P,Q,G) - Y'*R*Y;
        constr = [constr, L + 1e-10 <= 0];
    end
    
    % minimize trace of G to account for small coupling gains
    try
        result = optimize(constr,trace(G),sdpsettings('solver','mosek','verbose',0));
    catch ME
        fprintf('  solver error: %s\n',ME.message);
        return
    end
    if result.problem == 0
        status = 1;
        sol.P = value(P); sol.G = value(G); sol.Q = value(Q); sol.R = value(R);
    elseif result.problem == 1
        status = 0;
    else
        fprintf('  solver unclear: %s\n',result.info);
    end
end