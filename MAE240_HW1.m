%% Homework 1 Problem 4
%% Setup — Constants
clc;clear;close all
G = 6.67430e-20;   % km^3 kg^-1 s^-2
AU = 1.496e8;       % km
yr = 365.25 * 24 * 3600;  % s

%% Mass, Distance and Velocity
% Sun + Earth + Mars
n = 3;
m_bodies = [1.989e30, 5.972e24, 6.42e23];  % kg
M_sun=m_bodies(1);
M_earth=m_bodies(2);
M_mars=m_bodies(3);
r_AU=[0, 1, 1.524];

% Preallocate state-vector 
r0_all = zeros(1, 3*n);
v0_all = zeros(1, 3*n);

% Sun at the origin, initially at rest
r0_all(1:3) = [0, 0, 0];
v0_all(1:3) = [0, 0, 0];

% Planets start on circular, prograde, coplanar orbits at random phase
for i = 2:n
    r_km  = r_AU(i) * AU;
    theta = 2*pi * rand;
    r0_all((3*i-2):(3*i)) = [r_km*cos(theta), r_km*sin(theta), 0];
    v_c = sqrt(G * M_sun / r_km);
    v0_all((3*i-2):(3*i)) = [-v_c*sin(theta), v_c*cos(theta), 0];
end

% Assemble the initial condition as a column vector for ode45
x0 = [r0_all, v0_all].';

%% Integrate
T_end  = 10 * yr;
tspan  = linspace(0, T_end, 5000);
opts   = odeset('AbsTol', 1e-13, 'RelTol', 1e-11);

[T, X] = ode45(@(t,x) eom_NBP(t, x, m_bodies, G), tspan, x0, opts);

% Split state history into positions and velocities
r_all = X(:, 1:3*n);
v_all = X(:, 3*n+1:end);

%% Part A - Trajectory
figure; hold on; grid on; axis equal;
colors = {'y','b','r'};
names  = {'Sun','Earth','Mars'};
for i = 1:n
    ri = r_all(:, (3*i-2):(3*i));
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, ...
          'Color', colors{i}, 'DisplayName', names{i});
end
% Plot start and end markers separately so they appear in legend
h_start = gobjects(n,1);
h_end   = gobjects(n,1);
for i = 1:n
    ri = r_all(:, (3*i-2):(3*i));
    h_start(i) = plot3(ri(1,1)/AU, ri(1,2)/AU, ri(1,3)/AU, ...
        'o', 'Color', colors{i}, 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s start', names{i}));
    h_end(i) = plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, ...
        's', 'Color', colors{i}, 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s end', names{i}));
end

xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
legend('Location', 'best');
title('Trajectory of 3-Bodies');
view(3);
%% Part B- Center of Mass
m_total = sum(m_bodies);
cm_all = centerOfMassHistory(r_all,m_bodies);

% Linear Momentum Prediction
p_all= linearMomentumHistory(v_all,m_bodies);
Vcm0 = p_all(1,:)/m_total;
cm_ref= cm_all(1,:) + T*Vcm0;

% Max Deviation
cm_err= sqrt(sum((cm_all-cm_ref).^2,2)); %Normalized Error
fprintf('Center of mass at t=0   (AU): [%+.4e  %+.4e  %+.4e]\n', ...
    cm_all(1,1)/AU, cm_all(1,2)/AU, cm_all(1,3)/AU);
fprintf('Center of mass at t=T   (AU): [%+.4e  %+.4e  %+.4e]\n', ...
    cm_all(end,1)/AU, cm_all(end,2)/AU, cm_all(end,3)/AU);
fprintf('Max deviation from R_c(t_0)+V_c(t_0)(t-t_0): %.4e km  (%.4e AU)\n', ...
    max(cm_err), max(cm_err)/AU);

% Plot
figure;
labels = {'x','y','z'};
for k = 1:3
    subplot(3,1,k);
    plot(T/yr, cm_all(:,k)/AU, 'k-'); hold on;
    plot(T/yr, cm_ref(:,k)/AU, 'r--');
    ylabel(sprintf('R_{c,%s} [AU]', labels{k})); grid on;
    if k==1, legend('numerical','linear prediction'); end
end
xlabel('Time [yr]');
sgtitle("COM Numerical History vs Linear Prediction")

%% Part C — Angular Momentum
H_all = angularMomentumHistory(r_all, v_all, m_bodies);
dH    = H_all - H_all(1,:);
max_dH=max(abs(H_all),[],1);
rel_drift = sqrt(sum(dH.^2,2)) / norm(H_all(1,:));
fprintf('Max relative angular momentum drift: %.4e\n', max(rel_drift));
fprintf(['Max absolute angular momentum drfit:' ...
    ' Δp_x = %.4e, Δp_y = %.4e, Δp_z = %.4e\n'], max_dH(1), max_dH(2), max_dH(3));

figure;
labels = {'x','y','z'};
for k = 1:3
    subplot(3,1,k);
    plot(T/yr, H_all(:,k), 'k-'); grid on;
    ylabel(sprintf('H_%s [kg km^2/s]', labels{k}));
end
xlabel('Time [yr]');
sgtitle('Angular Momentum')

%% Part D — Total Energy
[KE_all, PE_all, E_all] = totalEnergyHistory(r_all, v_all, m_bodies, G);
rel_err_E = abs(E_all - E_all(1)) / abs(E_all(1));

figure;
subplot(2,1,1);
plot(T/yr, KE_all/abs(E_all(1)), 'b-'); hold on;
plot(T/yr, PE_all/abs(E_all(1)), 'r-');
plot(T/yr, E_all/abs(E_all(1)),  'k-');
legend('T','U','E'); ylabel('Energy / |E(t_0)|'); grid on;

subplot(2,1,2);
semilogy(T/yr, rel_err_E + eps, 'k-'); grid on;
ylabel('|E(t) - E(t_0)| / |E(t_0)|');
xlabel('Time [yr]'); title('Relative Energy Error');
sgtitle('Energy vs Time');

%% Helper Functions
function dxdt = eom_NBP(~,x,m_bodies,G)
%EOM_NBP Equations of motion for the Newtonian N-body problem.
n = numel(m_bodies);
r = x(1:3*n);
v = x(3*n+1:end);
dvdt = zeros(3*n,1);
for i = 1:n
    ri = r(3*(i-1)+1:3*i);
    ai = zeros(3,1);
    for j = 1:n
        if j ~= i
            mj = m_bodies(j);
            rj = r(3*(j-1)+1:3*j);
            r_ij = ri - rj;
            ai = ai - G * mj / norm(r_ij)^3 * r_ij;
        end
    end
    dvdt(3*(i-1)+1:3*i) = ai;
end
dxdt = [v; dvdt];
end

function cm_all = centerOfMassHistory(r_all, m_bodies)
%CENTEROFMASSHISTORY Center-of-mass history from stacked position history.
Nt = size(r_all,1);
n  = numel(m_bodies);
M  = sum(m_bodies);
cm_all = zeros(Nt,3);
for i = 1:n
    ri = r_all(:, (3*i-2):(3*i));
    cm_all = cm_all + m_bodies(i) * ri;
end
cm_all = cm_all / M;
end

function p_all = linearMomentumHistory(v_all, m_bodies)
%LINEARMOMENTUMHISTORY Total linear momentum history.
Nt = size(v_all,1);
n  = numel(m_bodies);
p_all = zeros(Nt,3);
for i = 1:n
    vi = v_all(:, (3*i-2):(3*i));
    p_all = p_all + m_bodies(i) * vi;
end
end

function H_all = angularMomentumHistory(r_all, v_all, m_bodies)
%ANGULARMOMENTUMHISTORY Total angular momentum history.
Nt = size(r_all,1);
n  = numel(m_bodies);
H_all = zeros(Nt,3);
for i = 1:n
    ri = r_all(:, (3*i-2):(3*i));
    vi = v_all(:, (3*i-2):(3*i));
    H_all = H_all + m_bodies(i) * cross(ri, vi, 2);
end
end

function [KE_all, PE_all, E_all] = totalEnergyHistory(r_all, v_all, m_bodies, G)
%TOTALENERGYHISTORY Kinetic, potential, and total energy histories.
Nt = size(r_all,1);
n  = numel(m_bodies);
KE_all = zeros(Nt,1);
for i = 1:n
    vi = v_all(:, (3*i-2):(3*i));
    KE_all = KE_all + 0.5 * m_bodies(i) * sum(vi.^2, 2);
end
PE_all = zeros(Nt,1);
for i = 1:n
    ri = r_all(:, (3*i-2):(3*i));
    for j = i+1:n
        rj = r_all(:, (3*j-2):(3*j));
        rij = rj - ri;
        d_ij = sqrt(sum(rij.^2, 2));
        PE_all = PE_all - G * m_bodies(i) * m_bodies(j) ./ d_ij;
    end
end
E_all = KE_all + PE_all;
end