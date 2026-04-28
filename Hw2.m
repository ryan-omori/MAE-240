% MAE 240 -- Part (a): Barycentric N-body integrator
% Equations of motion (eq. 2):
%   r_i_ddot = -sum_{j~=i} G*m_j / |r_i - r_j|^3 * (r_i - r_j)
% Verified by checking conservation of energy (eq. 3) and angular momentum (eq. 4).

close all; clearvars; clc;

%% Constants
G      = 6.67430e-20;           % [km^3 / (kg s^2)]
M_sun  = 1.989e30;              % [kg]
AU     = 1.496e8;               % [km]
yr     = 365.25 * 24 * 3600;   % [s]
mu_sun = G * M_sun;

%% Body setup: Sun, Venus, Earth, Mars, Jupiter
names = {'Sun','Venus','Earth','Mars','Jupiter'};
r_AU  = [0.00; 0.72; 1.00; 1.52; 5.20];
mass  = [M_sun; 4.867e24; 5.972e24; 6.417e23; 1.898e27];
N     = numel(mass);

%% Initial conditions: circular coplanar orbits, then shift to barycenter
rng(240);
r0 = zeros(3*N, 1);
v0 = zeros(3*N, 1);

for i = 2:N
    rmag  = r_AU(i) * AU;
    theta = 2*pi*(i-2)/(N-1);
    vc    = sqrt(mu_sun / rmag);
    r0(3*i-2:3*i) = [rmag*cos(theta); rmag*sin(theta); 0];
    v0(3*i-2:3*i) = [-vc*sin(theta);  vc*cos(theta);  0];
end

% Shift to true barycenter (makes frame inertial)
M_tot = sum(mass);
r_cm0 = zeros(3,1);
v_cm0 = zeros(3,1);
for i = 1:N
    r_cm0 = r_cm0 + mass(i) * r0(3*i-2:3*i);
    v_cm0 = v_cm0 + mass(i) * v0(3*i-2:3*i);
end
r_cm0 = r_cm0 / M_tot;
v_cm0 = v_cm0 / M_tot;

for i = 1:N
    r0(3*i-2:3*i) = r0(3*i-2:3*i) - r_cm0;
    v0(3*i-2:3*i) = v0(3*i-2:3*i) - v_cm0;
end

x0 = [r0; v0];

%% Integrate
tEnd  = 12 * yr;
tspan = linspace(0, tEnd, 5000);
opts  = odeset('AbsTol',1e-13,'RelTol',1e-11);

fprintf('Integrating %d bodies for 12 years...\n', N);
tic
[T, X] = ode45(@(t,x) eom_NBP_barycentric(t, x, mass, G), tspan, x0, opts);
fprintf('Done in %.2f s\n', toc);

r_all = X(:, 1:3*N);
v_all = X(:, 3*N+1:end);

%% Compute conservation integrals
[E, H] = compute_invariants(r_all, v_all, mass, G);

rel_dE = abs(E - E(1)) / abs(E(1));
rel_dH = vecnorm(H - H(1,:), 2, 2) / norm(H(1,:));

fprintf('Max relative energy drift    : %.4e\n', max(rel_dE));
fprintf('Max relative ang-mom drift   : %.4e\n', max(rel_dH));

%% Plots
colors = [0.95 0.85 0.10;
          0.60 0.80 0.95;
          0.20 0.70 0.30;
          0.85 0.40 0.20;
          0.70 0.50 0.85];

% Figure 1: Trajectories
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all(:, 3*i-2:3*i);
    plot(ri(:,1)/AU, ri(:,2)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', names{i});
    plot(ri(1,1)/AU, ri(1,2)/AU, 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'HandleVisibility','off');
end
xlabel('x [AU]'); ylabel('y [AU]');
title('Barycentric trajectories'); legend('Location','best');

nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all(:, 3*i-2:3*i);
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', names{i});
    % Initial position: filled circle
    plot3(ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   'o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
    % Final position: filled square
    plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, 's', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
end
% Dummy entries for marker legend
plot3(nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','Start');
plot3(nan,nan,nan,'ks','MarkerFaceColor','k','MarkerSize',6,'DisplayName','End');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('3D view'); legend('Location','best'); view(3);

sgtitle('Part (a) -- N-body barycentric integration');

% Figure 2: Conservation
figure('Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(T/yr, rel_dE, 'b-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|dE/E_0|');
title('Relative energy drift');

nexttile;
semilogy(T/yr, rel_dH, 'r-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|dH/H_0|');
title('Relative angular momentum drift');

nexttile;
plot(T/yr, E, 'b-', 'LineWidth', 1.5); hold on;
yline(E(1), 'k--', 'LineWidth', 1.0);
grid on; box on; xlabel('Time [yr]'); ylabel('E [km^2 kg s^{-2}]');
title('Total energy E(t)');

nexttile;
plot(T/yr, vecnorm(H,2,2), 'r-', 'LineWidth', 1.5); hold on;
yline(norm(H(1,:)), 'k--', 'LineWidth', 1.0);
grid on; box on; xlabel('Time [yr]'); ylabel('|H| [km^2 kg s^{-1}]');
title('Angular momentum |H(t)|');

sgtitle('Part (a) -- Conservation integrals (eqs. 3-4)');


%% Problem 2
% PROBLEM 2(a) -- Body-centric N-body integrator (eq. 6)
% =========================================================================
%
% Body-centric EOM for the N massive bodies (Sun = body 1 is the origin):
%   r_i_ddot = -G(m_N + m_i)/|r_i|^3 * r_i
%              - sum_{j~=i} G*m_j * [ (r_i-r_j)/|r_i-r_j|^3 + r_j/|r_j|^3 ]
%
% The frame is NOT inertial, so energy and angular momentum are NOT
% conserved -- this is compared directly against the barycentric results.
 
% Reuse the same planet initial conditions but Sun-centered (no barycenter shift)
r0_bc = zeros(3*N, 1);
v0_bc = zeros(3*N, 1);
for i = 2:N
    rmag  = r_AU(i) * AU;
    theta = 2*pi*(i-2)/(N-1);
    vc    = sqrt(mu_sun / rmag);
    r0_bc(3*i-2:3*i) = [rmag*cos(theta); rmag*sin(theta); 0];
    v0_bc(3*i-2:3*i) = [-vc*sin(theta);  vc*cos(theta);  0];
end
% Sun fixed at origin
r0_bc(1:3) = [0;0;0];
v0_bc(1:3) = [0;0;0];
 
x0_bc = [r0_bc; v0_bc];
 
fprintf('\nIntegrating body-centric N-body (eq. 6) for 12 years...\n');
tic
[T_bc, X_bc] = ode45(@(t,x) eom_NBP_bodycentric(t, x, mass, G), tspan, x0_bc, opts);
fprintf('Done in %.2f s\n', toc);
 
r_all_bc = X_bc(:, 1:3*N);
v_all_bc = X_bc(:, 3*N+1:end);
 
% Energy and angular momentum in the body-centric frame (will NOT be conserved)
[E_bc, H_bc] = compute_invariants(r_all_bc, v_all_bc, mass, G);
 
rel_dE_bc = abs(E_bc - E_bc(1)) / abs(E_bc(1));
rel_dH_bc = vecnorm(H_bc - H_bc(1,:), 2, 2) / norm(H_bc(1,:));
 
fprintf('Body-centric max relative energy drift    : %.4e\n', max(rel_dE_bc));
fprintf('Body-centric max relative ang-mom drift   : %.4e\n', max(rel_dH_bc));
fprintf('(Expected to be large -- frame is non-inertial)\n');
 
% Figure 3: Body-centric trajectories
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 
nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bc(:, 3*i-2:3*i);
    plot(ri(:,1)/AU, ri(:,2)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', names{i});
    plot(ri(1,1)/AU, ri(1,2)/AU, 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'HandleVisibility','off');
end
xlabel('x [AU]'); ylabel('y [AU]');
title('Body-centric trajectories (Sun-fixed)'); legend('Location','best');
 
nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bc(:, 3*i-2:3*i);
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', names{i});
    plot3(ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
    plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, 's', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
end
plot3(nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','Start');
plot3(nan,nan,nan,'ks','MarkerFaceColor','k','MarkerSize',6,'DisplayName','End');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('3D view'); legend('Location','best'); view(3);
 
sgtitle('Part 2(a) -- Body-centric N-body integration');
 
% Figure 4: Compare conservation between barycentric and body-centric
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 
nexttile; hold on; grid on; box on;
semilogy(T/yr,    rel_dE,    'b-',  'LineWidth', 1.5, 'DisplayName', 'Barycentric');
semilogy(T_bc/yr, rel_dE_bc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Body-centric');
xlabel('Time [yr]'); ylabel('|dE/E_0|');
title('Relative energy drift comparison'); legend('Location','best');
 
nexttile; hold on; grid on; box on;
semilogy(T/yr,    rel_dH,    'b-',  'LineWidth', 1.5, 'DisplayName', 'Barycentric');
semilogy(T_bc/yr, rel_dH_bc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Body-centric');
xlabel('Time [yr]'); ylabel('|dH/H_0|');
title('Relative angular momentum drift comparison'); legend('Location','best');
 
sgtitle('Part 2(a) -- Conservation: barycentric vs body-centric');
 
% =========================================================================
% PROBLEM 2(b) -- Restricted (N+1)st-body spacecraft in the body-centric frame
% =========================================================================
%
% Spacecraft EOM (eq. 5, restricted form -- m_sc << m_N):
%   r_sc_ddot = -G*m_N/|r_sc|^3 * r_sc
%               - sum_j G*m_j * [ (r_sc-r_j)/|r_sc-r_j|^3 + r_j/|r_j|^3 ]
%
% The massive-body ephemeris from 2(a) is interpolated to drive the spacecraft.
 
% Build interpolants from the body-centric integration (same as demo)
r_interp = makeInterpolants(T_bc, r_all_bc, N);
 
% Spacecraft: depart from Earth on a Hohmann-like transfer toward Jupiter
iEarth  = 3;
iJup    = 5;
r0_earth = r0_bc(3*iEarth-2:3*iEarth);
v0_earth = v0_bc(3*iEarth-2:3*iEarth);
r0_jup   = r0_bc(3*iJup-2:3*iJup);
 
r1_sc   = norm(r0_earth);
r2_goal = norm(r0_jup);
a_trans = 0.5*(r1_sc + r2_goal);
v_dep   = sqrt(mu_sun*(2/r1_sc - 1/a_trans));
u_t     = v0_earth / norm(v0_earth);
 
x0_sc = [r0_earth; v_dep * u_t];
 
fprintf('\nIntegrating restricted spacecraft (eq. 5)...\n');
tic
[T_sc, X_sc] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp, mass, G), tspan, x0_sc, opts);
fprintf('Done in %.2f s\n', toc);
 
r_sc = X_sc(:,1:3);
 
% Distance to Jupiter over time
r_jup_sc = reshape(r_interp{iJup}(T_sc), [], 3);
d_sc_jup = vecnorm(r_sc - r_jup_sc, 2, 2);
[~, idxClose] = min(d_sc_jup);
fprintf('Closest approach to Jupiter: %.4e AU at t = %.3f yr\n', ...
    d_sc_jup(idxClose)/AU, T_sc(idxClose)/yr);
 
% Figure 5: Spacecraft trajectory overlaid on body-centric planets
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 
nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bc(:, 3*i-2:3*i);
    plot(ri(:,1)/AU, ri(:,2)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', names{i});
end
plot(r_sc(:,1)/AU, r_sc(:,2)/AU, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Spacecraft');
plot(r_sc(1,1)/AU, r_sc(1,2)/AU, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 7, 'HandleVisibility','off');
xlabel('x [AU]'); ylabel('y [AU]');
title('Body-centric frame: planets + spacecraft'); legend('Location','best');
 
nexttile; hold on; grid on; box on;
plot(T_sc/yr, vecnorm(r_sc,2,2)/AU, 'k-', 'LineWidth', 1.5, 'DisplayName', 'SC radius');
yline(r2_goal/AU, 'r--', 'Jupiter orbit', 'LineWidth', 1.1, 'LabelHorizontalAlignment','left');
plot(T_sc(idxClose)/yr, d_sc_jup(idxClose)/AU, 'ro', 'MarkerFaceColor','r', 'DisplayName','Closest approach');
xlabel('Time [yr]'); ylabel('Distance [AU]');
title('Spacecraft heliocentric radius'); legend('Location','best');
 
sgtitle('Part 2(b) -- Restricted (N+1)^{st}-body spacecraft trajectory');
 
% -------------------------------------------------------------------------
% Helper functions
% -------------------------------------------------------------------------
 
function dx = eom_NBP_barycentric(~, x, mass, G)
% Barycentric N-body EOM (eq. 2). No indirect terms -- frame is inertial.
%   r_i_ddot = -sum_{j~=i} G*m_j / |r_i-r_j|^3 * (r_i-r_j)
N = numel(mass);
r = x(1:3*N);
v = x(3*N+1:end);
a = zeros(3*N, 1);
for i = 1:N
    ri = r(3*i-2:3*i);
    for j = 1:N
        if j == i, continue; end
        rj  = r(3*j-2:3*j);
        rij = ri - rj;
        a(3*i-2:3*i) = a(3*i-2:3*i) - G*mass(j)*rij / norm(rij)^3;
    end
end
dx = [v; a];
end
 
function [E, H] = compute_invariants(r_all, v_all, mass, G)
% Total energy (eq. 3) and angular momentum (eq. 4) at each time step.
Nt = size(r_all, 1);
N  = numel(mass);
E  = zeros(Nt, 1);
H  = zeros(Nt, 3);
for k = 1:Nt
    rk = reshape(r_all(k,:), 3, N);
    vk = reshape(v_all(k,:), 3, N);
    T  = 0.5 * sum(mass.' .* sum(vk.^2, 1));
    U  = 0;
    for i = 1:N
        for j = i+1:N
            U = U + G*mass(i)*mass(j) / norm(rk(:,i)-rk(:,j));
        end
    end
    E(k) = T - U;
    Hk = zeros(3,1);
    for i = 1:N
        Hk = Hk + mass(i)*cross(rk(:,i), vk(:,i));
    end
    H(k,:) = Hk.';
end
end
 
function dx = eom_NBP_bodycentric(~, x, mass, G)
% Body-centric N-body EOM (eq. 6). Body 1 (Sun) is fixed at the origin.
%   r_i_ddot = -G(m_1+m_i)/|r_i|^3 * r_i
%              - sum_{j~=i} G*m_j*[ (r_i-r_j)/|r_i-r_j|^3 + r_j/|r_j|^3 ]
N = numel(mass);
r = x(1:3*N);
v = x(3*N+1:end);
a = zeros(3*N, 1);
% Body 1 (Sun) stays fixed
a(1:3) = [0;0;0];
for i = 2:N
    ri   = r(3*i-2:3*i);
    % Direct term
    a_i = -G*(mass(1)+mass(i)) * ri / norm(ri)^3;
    % Indirect terms from other non-central bodies
    for j = 2:N
        if j == i, continue; end
        rj  = r(3*j-2:3*j);
        rij = ri - rj;
        a_i = a_i - G*mass(j) * ( rij/norm(rij)^3 + rj/norm(rj)^3 );
    end
    a(3*i-2:3*i) = a_i;
end
dx = [v; a];
end
 
function dx = eom_sc_bodycentric(t, x, r_interp, mass, G)
% Restricted (N+1)st-body EOM in the body-centric frame (eq. 5).
%   r_sc_ddot = -G*m_1/|r_sc|^3 * r_sc
%               - sum_j G*m_j*[ (r_sc-r_j)/|r_sc-r_j|^3 + r_j/|r_j|^3 ]
N    = numel(mass);
r_sc = x(1:3);
v_sc = x(4:6);
% Direct solar term
a_sc = -G*mass(1) * r_sc / norm(r_sc)^3;
% Indirect planetary terms
for j = 2:N
    rj   = reshape(r_interp{j}(t), [3,1]);
    rscj = r_sc - rj;
    a_sc = a_sc - G*mass(j) * ( rscj/norm(rscj)^3 + rj/norm(rj)^3 );
end
dx = [v_sc; a_sc];
end
 
function r_interp = makeInterpolants(T, r_all, N)
% Build spline interpolants for each body's position history.
r_interp = cell(N,1);
for i = 1:N
    ri = r_all(:, 3*i-2:3*i);
    r_interp{i} = griddedInterpolant(T, ri, 'spline', 'nearest');
end
end