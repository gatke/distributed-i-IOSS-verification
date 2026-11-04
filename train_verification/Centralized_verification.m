% -------------------------------------------------------------------------
% File: Centralized_verification.m
%
% Description:
%   Centralized verification of the train model for different numbers
%   of carriages.
%
% Reference:
%   Christian Gatke, Julian D. Schiller, and Matthias A. Müller,
%   "Small-gain analysis of incremental input/output-to-state stability
%   for large-scale distributed systems", 2026 (not yet published).
%
% -------------------------------------------------------------------------


clear; clc; yalmip('clear');

M_grid = [3 4 5 6 7 8 9 10 11 12 13 14]; % number of subsystems
eta = 0.99;
p.dt = 0.01; p.m = 90e3; p.k = 1e7; p.d = 5e5;
p.bound = 0.8; % bound for relative velocity between carriages

diagnostics = cell(numel(M_grid),1);
for ic = 1:numel(M_grid)
    M = M_grid(ic);
    fprintf('\n--- M = %d ---\n',M);
    [Avertices,B,C,D] = overall_system(M,p);
    n = 2*M; nw = 3*M; ny = M;
    nv = numel(Avertices);

    P = sdpvar(n,n,'symmetric');
    Q = sdpvar(nw,nw,'diagonal');
    R = sdpvar(ny,ny,'diagonal');

    constr = [trace(P)==n, P>=1e-6*eye(n), Q>=0, R>=0];
    Y = [C D];
    for iv = 1:nv
        J = [Avertices{iv} B];
        L = J'*P*J - blkdiag(eta*P,Q) - Y'*R*Y;
        constr = [constr, L <= 0];
    end

    tOpt = tic;
    try
        diagnostics{ic} = optimize(constr,[],sdpsettings('solver','mosek','verbose',0));
        optimizeTime = toc(tOpt);
        status = diagnostics{ic}.problem;
    catch ME
        status = -1;
        fprintf('MATLAB error: %s\n',ME.message);
    end

    if status == 0
        fprintf('SDP successfully solved.\n');
        fprintf('Optimization time: %.3f\n', optimizeTime);
    else
        fprintf('No solution found.\n');
    end
end



function [Avertices,B,C,D] = overall_system(M,p)
    nEdges = M-1;
    nVertices = 2^nEdges;

    a = p.dt*p.k/p.m;
    b = 3*p.dt*p.d/p.m;
    A0 = kron(eye(M),[1 p.dt; 0 1]);
    B = kron(eye(M),[1 0 0; 0 1 0]);
    C = kron(eye(M),[1 0]);
    D = kron(eye(M),[0 0 1]);

    Avertices = cell(nVertices,1);
    for v = 0:nVertices-1
        A = A0;
        for j = 1:nEdges
            % one bit per coupling: 0 or p.bound^2
            theta = double(bitget(uint32(v),j))*p.bound^2;

            posL = 2*j-1; velL = 2*j;
            posR = 2*(j+1)-1; velR = 2*(j+1);

            A(velL,posL) = A(velL,posL)-a;
            A(velL,posR) = A(velL,posR)+a;
            A(velR,posR) = A(velR,posR)-a;
            A(velR,posL) = A(velR,posL)+a;

            A(velL,velL) = A(velL,velL)-b*theta;
            A(velL,velR) = A(velL,velR)+b*theta;
            A(velR,velR) = A(velR,velR)-b*theta;
            A(velR,velL) = A(velR,velL)+b*theta;
        end
        Avertices{v+1} = A;
    end
end
