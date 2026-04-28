% MAE 240 -- Problems 1 and 2: Barycentric and Body-centric N-body integrators
close all; clearvars; clc;

%% Constants
G      = 6.67430e-20;
M_sun  = 1.989e30;
AU     = 1.496e8;
yr     = 365.25 * 24 * 3600;
mu_sun = G * M_sun;

%% Body setup: Sun, Venus, Earth, Mars, Jupiter
names = {'Sun','Venus','Earth','Mars','Jupiter'};
r_AU  = [0.00; 0.72; 1.00; 1.52; 5.20];
mass  = [M_sun; 4.867e24; 5.972e24; 6.417e23; 1.898e27];
N     = numel(mass);

colors = [0.95 0.85 0.10;
          0.60 0.80 0.95;
          0.20 0.70 0.30;
          0.85 0.40 0.20;
          0.70 0.50 0.85];


%% PROBLEM 1(a) -- Barycentric N-body integrator
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

% Shift to barycenter
M_tot = sum(mass);
r_cm0 = zeros(3,1); v_cm0 = zeros(3,1);
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

tEnd  = 12 * yr;
tspan = linspace(0, tEnd, 5000);
opts  = odeset('AbsTol',1e-13,'RelTol',1e-11);

fprintf('=== Problem 1(a): Barycentric N-body ===\n');
tic
[T_bary, X_bary] = ode45(@(t,x) eom_NBP_barycentric(t, x, mass, G), tspan, x0, opts);
fprintf('Done in %.2f s\n', toc);

r_all_bary = X_bary(:, 1:3*N);
v_all_bary = X_bary(:, 3*N+1:end);

% Conservation -- barycentric (should be tightly conserved)
[E_bary, H_bary] = compute_invariants(r_all_bary, v_all_bary, mass, G);
rel_dE_bary = abs(E_bary - E_bary(1)) / abs(E_bary(1));
rel_dH_bary = vecnorm(H_bary - H_bary(1,:), 2, 2) / norm(H_bary(1,:));
fprintf('Max relative energy drift    : %.4e\n', max(rel_dE_bary));
fprintf('Max relative ang-mom drift   : %.4e\n', max(rel_dH_bary));

%% PROBLEM 1(b) -- Restricted barycentric spacecraft
iEarth_1b  = 3;
iJup_1b    = 5;
r0_earth_1b = r0(3*iEarth_1b-2:3*iEarth_1b);
v0_earth_1b = v0(3*iEarth_1b-2:3*iEarth_1b);
r0_jup_1b   = r0(3*iJup_1b-2:3*iJup_1b);

u_r_1b  = r0_earth_1b / norm(r0_earth_1b);
u_t_1b  = v0_earth_1b / norm(v0_earth_1b);
r0_sc_1b = r0_earth_1b + 1e-3 * AU * u_r_1b;

r1_sc_1b   = norm(r0_sc_1b);
r2_goal_1b = norm(r0_jup_1b);
a_trans_1b = 0.5*(r1_sc_1b + r2_goal_1b);
T_trans_1b = pi * sqrt(a_trans_1b^3 / mu_sun);
v_dep_1b   = sqrt(mu_sun*(2/r1_sc_1b - 1/a_trans_1b));

r_interp_bary = makeInterpolants(T_bary, r_all_bary, N);

% Bisection tuning -- 1(b)
sc_lo_1b = 0.98; sc_hi_1b = 1.05;
bestScale_1b = sc_lo_1b; bestGap_1b = inf;

fprintf('\n=== Problem 1(b): Barycentric spacecraft tuning ===\n');
for k = 1:20
    sc = 0.5*(sc_lo_1b + sc_hi_1b);
    x0_sc_1b = [r0_sc_1b; sc * v_dep_1b * u_t_1b];
    [T_tmp, X_tmp] = ode45(@(t,x) eom_sc_barycentric(t, x, r_interp_bary, mass, G), tspan, x0_sc_1b, opts);
    [apoR, ~, ~]   = firstApoapsis(T_tmp, X_tmp(:,1:3), X_tmp(:,4:6), T_trans_1b, 0.35, 1.75);
    if apoR <= r2_goal_1b
        gap = r2_goal_1b - apoR;
        if gap < bestGap_1b, bestGap_1b = gap; bestScale_1b = sc; end
        sc_lo_1b = sc;
    else
        sc_hi_1b = sc;
    end
end
fprintf('Best speed scale 1(b): %.6f\n', bestScale_1b);

x0_sc_1b = [r0_sc_1b; bestScale_1b * v_dep_1b * u_t_1b];
[T_sc_1b, X_sc_1b] = ode45(@(t,x) eom_sc_barycentric(t, x, r_interp_bary, mass, G), tspan, x0_sc_1b, opts);
r_sc_1b = X_sc_1b(:,1:3);
v_sc_1b = X_sc_1b(:,4:6);

% Spacecraft specific energy and angular momentum -- barycentric frame
% Note: uses Sun-only two-body proxy since spacecraft mass is negligible
Nt_1b    = size(r_sc_1b, 1);
E_sc_1b  = zeros(Nt_1b, 1);
H_sc_1b  = zeros(Nt_1b, 3);
for k = 1:Nt_1b
    rk = r_sc_1b(k,:);
    vk = v_sc_1b(k,:);
    % Shift r to be relative to Sun position in barycentric frame
    r_sun_k = reshape(r_interp_bary{1}(T_sc_1b(k)), [1,3]);
    r_rel   = rk - r_sun_k;
    E_sc_1b(k)   = 0.5*dot(vk,vk) - G*mass(1)/norm(r_rel);
    H_sc_1b(k,:) = cross(r_rel, vk);
end
rel_dE_sc_1b = abs(E_sc_1b - E_sc_1b(1)) / abs(E_sc_1b(1));
rel_dH_sc_1b = vecnorm(H_sc_1b - H_sc_1b(1,:), 2, 2) / norm(H_sc_1b(1,:));
fprintf('Max relative SC energy drift 1(b)  : %.4e\n', max(rel_dE_sc_1b));
fprintf('Max relative SC ang-mom drift 1(b) : %.4e\n', max(rel_dH_sc_1b));

%% PROBLEM 2(a) -- Body-centric N-body integrator
% Fresh ICs: Sun at origin, no barycenter shift
r0_bc = zeros(3*N, 1);
v0_bc = zeros(3*N, 1);
for i = 2:N
    rmag  = r_AU(i) * AU;
    theta = 2*pi*(i-2)/(N-1);
    vc    = sqrt(mu_sun / rmag);
    r0_bc(3*i-2:3*i) = [rmag*cos(theta); rmag*sin(theta); 0];
    v0_bc(3*i-2:3*i) = [-vc*sin(theta);  vc*cos(theta);  0];
end
r0_bc(1:3) = [0;0;0];
v0_bc(1:3) = [0;0;0];
x0_bc = [r0_bc; v0_bc];

fprintf('\n=== Problem 2(a): Body-centric N-body ===\n');
tic
[T_bc, X_bc] = ode45(@(t,x) eom_NBP_bodycentric(t, x, mass, G), tspan, x0_bc, opts);
fprintf('Done in %.2f s\n', toc);

r_all_bc = X_bc(:, 1:3*N);
v_all_bc = X_bc(:, 3*N+1:end);

% Conservation -- body-centric (expected to drift: non-inertial frame)
[E_bc, H_bc] = compute_invariants(r_all_bc, v_all_bc, mass, G);
rel_dE_bc = abs(E_bc - E_bc(1)) / abs(E_bc(1));
rel_dH_bc = vecnorm(H_bc - H_bc(1,:), 2, 2) / norm(H_bc(1,:));
fprintf('Body-centric max relative energy drift    : %.4e\n', max(rel_dE_bc));
fprintf('Body-centric max relative ang-mom drift   : %.4e\n', max(rel_dH_bc));

%% PROBLEM 2(b) -- Restricted body-centric spacecraft

r_interp_bc = makeInterpolants(T_bc, r_all_bc, N);

iEarth_2b   = 3;
iJup_2b     = 5;
r0_earth_2b = r0_bc(3*iEarth_2b-2:3*iEarth_2b);
v0_earth_2b = v0_bc(3*iEarth_2b-2:3*iEarth_2b);
r0_jup_2b   = r0_bc(3*iJup_2b-2:3*iJup_2b);

u_r_2b   = r0_earth_2b / norm(r0_earth_2b);
u_t_2b   = v0_earth_2b / norm(v0_earth_2b);
r0_sc_2b = r0_earth_2b + 1e-3 * AU * u_r_2b;

r1_sc_2b   = norm(r0_sc_2b);
r2_goal_2b = norm(r0_jup_2b);
a_trans_2b = 0.5*(r1_sc_2b + r2_goal_2b);
T_trans_2b = pi * sqrt(a_trans_2b^3 / mu_sun);
v_dep_2b   = sqrt(mu_sun*(2/r1_sc_2b - 1/a_trans_2b));

% Bisection tuning -- 2(b), completely fresh bracket
sc_lo_2b = 0.98; sc_hi_2b = 1.05;
bestScale_2b = sc_lo_2b; bestGap_2b = inf;

fprintf('\n=== Problem 2(b): Body-centric spacecraft tuning ===\n');
for k = 1:20
    sc = 0.5*(sc_lo_2b + sc_hi_2b);
    x0_sc_2b = [r0_sc_2b; sc * v_dep_2b * u_t_2b];
    [T_tmp, X_tmp] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp_bc, mass, G), tspan, x0_sc_2b, opts);
    [apoR_bc, ~, ~] = firstApoapsis(T_tmp, X_tmp(:,1:3), X_tmp(:,4:6), T_trans_2b, 0.35, 1.75);
    if apoR_bc <= r2_goal_2b
        gap = r2_goal_2b - apoR_bc;
        if gap < bestGap_2b, bestGap_2b = gap; bestScale_2b = sc; end
        sc_lo_2b = sc;
    else
        sc_hi_2b = sc;
    end
end
fprintf('Best speed scale 2(b): %.6f\n', bestScale_2b);

x0_sc_2b = [r0_sc_2b; bestScale_2b * v_dep_2b * u_t_2b];
[T_sc_2b, X_sc_2b] = ode45(@(t,x) eom_sc_bodycentric(t, x, r_interp_bc, mass, G), tspan, x0_sc_2b, opts);
r_sc_2b = X_sc_2b(:,1:3);
v_sc_2b = X_sc_2b(:,4:6);

% Spacecraft specific energy and angular momentum -- body-centric frame
% r_sc is already relative to Sun (origin), so no shift needed
Nt_2b    = size(r_sc_2b, 1);
E_sc_2b  = zeros(Nt_2b, 1);
H_sc_2b  = zeros(Nt_2b, 3);
for k = 1:Nt_2b
    rk = r_sc_2b(k,:);
    vk = v_sc_2b(k,:);
    E_sc_2b(k)   = 0.5*dot(vk,vk) - G*mass(1)/norm(rk);
    H_sc_2b(k,:) = cross(rk, vk);
end
rel_dE_sc_2b = abs(E_sc_2b - E_sc_2b(1)) / abs(E_sc_2b(1));
rel_dH_sc_2b = vecnorm(H_sc_2b - H_sc_2b(1,:), 2, 2) / norm(H_sc_2b(1,:));
fprintf('Max relative SC energy drift 2(b)  : %.4e\n', max(rel_dE_sc_2b));
fprintf('Max relative SC ang-mom drift 2(b) : %.4e\n', max(rel_dH_sc_2b));

%% FIGURES
% Figure 1: Problem 1 trajectories 
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bary(:, 3*i-2:3*i);
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', names{i});
    plot3(ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
    plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, 's', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
end
plot3(nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','Start');
plot3(nan,nan,nan,'ks','MarkerFaceColor','k','MarkerSize',6,'DisplayName','End');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Problem 1(a) -- Barycentric N-body trajectories');
legend('Location','best'); view(3);

nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bary(:, 3*i-2:3*i);
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', names{i});
    plot3(ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
    plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, 's', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
end
plot3(r_sc_1b(:,1)/AU, r_sc_1b(:,2)/AU, r_sc_1b(:,3)/AU, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Spacecraft');
plot3(r_sc_1b(1,1)/AU,   r_sc_1b(1,2)/AU,   r_sc_1b(1,3)/AU,   'ko', 'MarkerFaceColor','k', 'MarkerSize', 7, 'HandleVisibility','off');
plot3(r_sc_1b(end,1)/AU, r_sc_1b(end,2)/AU, r_sc_1b(end,3)/AU, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 7, 'HandleVisibility','off');
plot3(nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','Start');
plot3(nan,nan,nan,'ks','MarkerFaceColor','k','MarkerSize',6,'DisplayName','End');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Problem 1(b) -- Restricted (N+1)^{th}-body spacecraft, barycentric');
legend('Location','best'); view(3);

sgtitle('Problem 1 -- Barycentric formulation trajectories');

% Figure 2: Problem 1 conservation
figure('Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(T_bary/yr, rel_dE_bary, 'b-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaE / E_0|');
title('1(a) N-body: relative energy drift');

nexttile;
semilogy(T_bary/yr, rel_dH_bary, 'r-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaH / H_0|');
title('1(a) N-body: relative ang-mom drift');

nexttile;
semilogy(T_sc_1b/yr, rel_dE_sc_1b, 'b-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaE / E_0|');
title('1(b) Spacecraft: specific energy drift (not conserved -- planetary perturbations)');

nexttile;
semilogy(T_sc_1b/yr, rel_dH_sc_1b, 'r-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaH / H_0|');
title('1(b) Spacecraft: specific ang-mom drift (not conserved -- planetary perturbations)');

sgtitle('Problem 1 -- Conservation: 1(a) integration accuracy and 1(b) non-conservation');

% Figure 3: Problem 2 trajectories
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

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
title('Problem 2(a) -- Body-centric N-body trajectories');
legend('Location','best'); view(3);

nexttile; hold on; grid on; box on; axis equal;
for i = 1:N
    ri = r_all_bc(:, 3*i-2:3*i);
    plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, '-', 'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', names{i});
    plot3(ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
    plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, 's', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 6, 'HandleVisibility','off');
end
plot3(r_sc_2b(:,1)/AU, r_sc_2b(:,2)/AU, r_sc_2b(:,3)/AU, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Spacecraft');
plot3(r_sc_2b(1,1)/AU,   r_sc_2b(1,2)/AU,   r_sc_2b(1,3)/AU,   'ko', 'MarkerFaceColor','k', 'MarkerSize', 7, 'HandleVisibility','off');
plot3(r_sc_2b(end,1)/AU, r_sc_2b(end,2)/AU, r_sc_2b(end,3)/AU, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 7, 'HandleVisibility','off');
plot3(nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6,'DisplayName','Start');
plot3(nan,nan,nan,'ks','MarkerFaceColor','k','MarkerSize',6,'DisplayName','End');
xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
title('Problem 2(b) -- Restricted (N+1)^{th}-body spacecraft, body-centric');
legend('Location','best'); view(3);

sgtitle('Problem 2 -- Body-centric formulation trajectories');

% Figure 4: Problem 2 conservation
figure('Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(T_bc/yr, rel_dE_bc, 'b-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaE / E_0|');
title('2(a) N-body: relative energy drift (not conserved -- non-inertial frame)');

nexttile;
semilogy(T_bc/yr, rel_dH_bc, 'r-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaH / H_0|');
title('2(a) N-body: relative ang-mom drift (not conserved -- non-inertial frame)');

nexttile;
semilogy(T_sc_2b/yr, rel_dE_sc_2b, 'b-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaE / E_0|');
title('2(b) Spacecraft: specific energy drift');

nexttile;
semilogy(T_sc_2b/yr, rel_dH_sc_2b, 'r-', 'LineWidth', 1.5);
grid on; box on; xlabel('Time [yr]'); ylabel('|\DeltaH / H_0|');
title('2(b) Spacecraft: specific ang-mom drift');

sgtitle('Problem 2 -- Conservation: 2(a) non-conservation expected, 2(b) spacecraft drift');

%% HELPER FUNCTIONS
function dx = eom_NBP_barycentric(~, x, mass, G)
N = numel(mass);
r = x(1:3*N); v = x(3*N+1:end);
a = zeros(3*N, 1);
for i = 1:N
    ri = r(3*i-2:3*i);
    for j = 1:N
        if j == i, continue; end
        rij = ri - r(3*j-2:3*j);
        a(3*i-2:3*i) = a(3*i-2:3*i) - G*mass(j)*rij / norm(rij)^3;
    end
end
dx = [v; a];
end

function dx = eom_NBP_bodycentric(~, x, mass, G)
% Body 1 (Sun) fixed at origin. Direct + indirect terms (eq. 6).
N = numel(mass);
r = x(1:3*N); v = x(3*N+1:end);
a = zeros(3*N, 1);
for i = 2:N
    ri  = r(3*i-2:3*i);
    a_i = -G*(mass(1)+mass(i)) * ri / max(norm(ri)^3, eps);
    for j = 2:N
        if j == i, continue; end
        rj  = r(3*j-2:3*j);
        rij = ri - rj;
        a_i = a_i - G*mass(j) * ( rij/max(norm(rij)^3,eps) + rj/max(norm(rj)^3,eps) );
    end
    a(3*i-2:3*i) = a_i;
end
dx = [v; a];
end

function dx = eom_sc_barycentric(t, x, r_interp, mass, G)
% No indirect terms -- inertial frame (eq. 1).
N    = numel(mass);
r_sc = x(1:3); v_sc = x(4:6);
a_sc = zeros(3,1);
for i = 1:N
    ri   = reshape(r_interp{i}(t), [3,1]);
    r_sci = r_sc - ri;
    a_sc = a_sc - G*mass(i) * r_sci / norm(r_sci)^3;
end
dx = [v_sc; a_sc];
end

function dx = eom_sc_bodycentric(t, x, r_interp, mass, G)
% Direct + indirect terms -- non-inertial frame (eq. 5).
N    = numel(mass);
r_sc = x(1:3); v_sc = x(4:6);
a_sc = -G*mass(1) * r_sc / max(norm(r_sc)^3, eps);
for j = 2:N
    rj   = reshape(r_interp{j}(t), [3,1]);
    rscj = r_sc - rj;
    a_sc = a_sc - G*mass(j) * ( rscj/max(norm(rscj)^3,eps) + rj/max(norm(rj)^3,eps) );
end
dx = [v_sc; a_sc];
end

function [E, H] = compute_invariants(r_all, v_all, mass, G)
Nt = size(r_all, 1);
N  = numel(mass);
E  = zeros(Nt, 1);
H  = zeros(Nt, 3);
for k = 1:Nt
    rk = reshape(r_all(k,:), 3, N);
    vk = reshape(v_all(k,:), 3, N);
    KE = 0.5 * sum(mass.' .* sum(vk.^2, 1));
    U  = 0;
    for i = 1:N
        for j = i+1:N
            U = U + G*mass(i)*mass(j) / norm(rk(:,i)-rk(:,j));
        end
    end
    E(k) = KE - U;
    Hk = zeros(3,1);
    for i = 1:N
        Hk = Hk + mass(i)*cross(rk(:,i), vk(:,i));
    end
    H(k,:) = Hk.';
end
end

function r_interp = makeInterpolants(T, r_all, N)
r_interp = cell(N,1);
for i = 1:N
    ri = r_all(:, 3*i-2:3*i);
    r_interp{i} = griddedInterpolant(T, ri, 'spline', 'nearest');
end
end

function [rApo, idxApo, tApo] = firstApoapsis(T, r_hist, v_hist, T_guess, fracLo, fracHi)
r_norm = vecnorm(r_hist, 2, 2);
rdot   = sum(r_hist .* v_hist, 2) ./ max(r_norm, eps);
mask   = (T >= fracLo*T_guess) & (T <= fracHi*T_guess);
idxWindow = find(mask);
idxApo = [];
if numel(idxWindow) >= 3
    for kk = idxWindow(1)+1 : idxWindow(end)-1
        if rdot(kk-1) >= 0 && rdot(kk+1) <= 0
            i1 = max(1, kk-2); i2 = min(numel(T), kk+2);
            [~, loc] = max(r_norm(i1:i2));
            idxApo = i1 + loc - 1;
            break
        end
    end
end
if isempty(idxApo)
    [~, loc] = max(r_norm(idxWindow));
    idxApo   = idxWindow(loc);
end
rApo = r_norm(idxApo);
tApo = T(idxApo);
end