% -------------------------------------------------------------------------
% File: D1_train_simulation.m
%
% Description:
%   Simulation of the acceleration process of a train.
%
% Reference:
%   Christian Gatke, Julian D. Schiller, and Matthias A. Müller,
%   "Small-gain analysis of incremental input/output-to-state stability
%   for large-scale distributed systems", 2026 (not yet published).
%
% -------------------------------------------------------------------------

clear; close all; clc;

dt = 0.01; % discretization step
N = 3e4; % number of time steps
Fmax = 15e4; % maximum force
m = 90e3; k = 1e7; d = 5e5;

M = 5; % number of subsystems
n = 2; tt = N*dt; disp("simulation time [min]:"); disp(tt/60);
x = zeros(2*M,N+1);

% disturbance
w = zeros(2*M,N);
w(1:2:end,:) = dt*0.01*randn(M,N);
w(2:2:end,:) = dt*0.001*randn(M,N);

% linear increase of the Force until Fmax
F = zeros(1,N); tF = round(0.75*tt); tN = tF/dt;
for t = 1:N
    if t <= tN
        F(t) = Fmax*t/tN;
    else
        F(t) = Fmax;
    end
    x(1:2,t+1) = f1(x(1:2,t),w(1:2,t),x(3:4,t),F(t),dt,m,k,d);
    x(M*n-1:M*n,t+1) = fM(x(M*n-1:M*n,t),w(M*n-1:M*n,t),x(M*n-3:M*n-2,t),dt,m,k,d);
    for j = 2:M-1
        z = [x(j*n-3:j*n-2,t); x(j*n+1:j*n+2,t)];
        x(j*n-1:j*n,t+1) = fi(x(j*n-1:j*n,t),w(j*n-1:j*n,t),z(:),dt,m,k,d);
    end
end

disp("distance [km]:"); disp(x(1,end)/1000);
disp("velocity [km/h]:"); disp(x(2,end)*3.6);

% determine maximum bound for relative velocity between carriages
bound = 0;
for j = 1:M-1
    dv = x(2*j,:) - x(2*j+2,:);
    bound = max(bound,max(abs(dv)));
end
disp('max. relative velocity[m/s]:'); disp(bound);



function x_next = f1(x,w,z,F,dt,m,k,d)
    x_next = [x(1) + dt*x(2) + w(1); ...
        x(2) + dt/m*(F+k*(z(1)-x(1))+d*(z(2)-x(2))^3) + w(2)];
end

function x_next = fi(x,w,z,dt,m,k,d)
    x_next = [x(1) + dt*x(2) + w(1); ...
        x(2) + dt/m*(k*(z(1)-x(1))+d*(z(2)-x(2))^3+k*(z(3)-x(1))+d*(z(4)-x(2))^3) + w(2)];
end

function x_next = fM(x,w,z,dt,m,k,d)
    x_next = [x(1) + dt*x(2) + w(1); ...
        x(2) + dt/m*(k*(z(1)-x(1))+d*(z(2)-x(2))^3) + w(2)];
end