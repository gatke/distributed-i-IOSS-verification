% -------------------------------------------------------------------------
% File: D3_1_LMI_check_outer.m
%
% Description:
%   Check satisfaction of LMI conditions for the boundary carriages.
%
% Reference:
%   Christian Gatke, Julian D. Schiller, and Matthias A. Müller,
%   "Small-gain analysis of incremental input/output-to-state stability
%   for large-scale distributed systems", 2026 (not yet published).
%
% -------------------------------------------------------------------------


clear; clc;

% choose the file you want to check
load('LyapOpti_outer.mat');
% load('TrajectOpti_outer.mat');

eta = eta_e; G = G_e; P = P_e; Q = Q_e; R = R_e;

% check symmetry
if ~(issymmetric(P) && issymmetric(Q) && issymmetric(G))
    fprintf('at least one matrix is not symmetric\n')
end

% check positive (semi)definiteness
if min(eig(P)) <= 0 || min(eig(Q)) < 0 || min(eig(G)) < 0
    fprintf('at least one matrix is not positive (semi)definite.\n')
end

% parameters and setting
p.dt = 0.01; p.m = 90e3; p.k = 1e7; p.d = 5e5;
a = p.dt*p.k/p.m; b = 3*p.dt*p.d/p.m;
p.bound = 0.8; % bound for relative velocity between carriages
vertices = [0 p.bound^2];

B = [1 0 0; 0 1 0]; C = [1 0]; D = [0 0 1];

maxEig = -inf; lmiOK = true;
for k = 1:size(vertices,2)
    A = [1 p.dt; -a 1-b*vertices(k)];
    E = [0 0; a b*vertices(k)];

    % LMI conditions
    J = [A B E];
    Y = [C D zeros(1,2)];
    L = J'*P*J - blkdiag(eta*P,Q,G) - Y'*R*Y;
    maxEig = max(maxEig,max(eig(L)));
    passed = maxEig <= 0;
    if ~passed
        lmiOK = false;
        fprintf('eigenvalue: %.15f\n', maxEig)
        fprintf('LMI condition is not satisfied.\n')
        break;
    end
end

if lmiOK
    fprintf('Largest eigenvalue: %.15f\n', maxEig)
    fprintf('LMI condtions are satisfied.\n')
end

