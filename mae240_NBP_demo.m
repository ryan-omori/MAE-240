%{
FILENAME:  mae240_NBP_demo.m

COURSE:    MAE 240 — Space Flight Mechanics (Spring 2026)
REFERENCE: mae240_s26_newtonian_annotated.pdf  (lecture notes only)

DESCRIPTION:
    A sectioned (%%) MATLAB script for the MAE 240 lecture on the Newtonian
    N-body problem. The script follows the lecture sequence:
      • Newtonian equations of motion and state-space form
      • numerical integration of the initial-value problem
      • barycenter / center-of-mass behavior
      • conservation of linear momentum, angular momentum, and energy
      • comparison with a separate randomized many-body realization

    The main run (Sections 3–9) is intentionally Solar-System-like:
      one dominant central body plus four pseudo-planets.

    Section 10 then switches to a randomized compact cluster so the final
    figure is more dynamically interesting than simple nearly straight-line
    departures.

INSTRUCTOR NOTES:
    - Run section-by-section (Ctrl+Enter) in lecture, or run the whole file.
    - All figures are guarded by doPlot.
    - eom_NBP is included below as a local helper function.
    - The lecture notes are the reference; there is no textbook for MAE 240.
%}

%% Section 0 — Global controls (read first)
close all
clearvars
clc

HK.enabled         = true;
HK.doClc           = false;
HK.doClear         = false;
HK.doClose         = false;
HK.pauseAtSections = true;
HK.useKeyboard     = true;
HK.clcAfterSection = true;
HK.clearVarsKeep   = {'HK','fmtTitle','doPlot'};

doPlot = true;

format compact;

fmtTitle = defaultFmtTitle();
fmtTitle('MAE 240 DEMO — Newtonian N-Body Problem');

disp('Reference: mae240_s26_newtonian_annotated.pdf');
disp('Lecture concepts: EOM/state vector, barycenter, linear momentum, angular momentum, energy.');
disp('Sections 3–9 keep the Solar-System-like run. Section 10 uses separate randomized ICs.');

hkPause(HK);

%% Section 1 — The N-body problem: equations of motion and the state vector
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 1 — The N-body problem: equations of motion and the state vector');

%{
  LECTURE CONNECTION
  ==================
  In the lecture notes, the Newtonian N-body equations are written as

      m_i r¨_i = Σ_{j≠i} -G m_i m_j (r_i - r_j) / |r_i - r_j|^3

  which gives

      r¨_i = Σ_{j≠i} -G m_j (r_i - r_j) / |r_i - r_j|^3 .

  To integrate numerically, we rewrite the problem in first-order form with

      x = [r_1; ...; r_n; v_1; ...; v_n] ∈ R^{6n}

  and then dx/dt = [v_1; ...; v_n; a_1; ...; a_n].
%}

disp('Local helper at end of file: dxdt = eom_NBP(t, x, m_bodies, G)');
disp('  x = [r_all; v_all]  with size 6n x 1');
disp('  r_all = [r1; r2; ...; rn],  v_all = [v1; v2; ...; vn]');

hkPause(HK);

%% Section 2 — Physical constants and unit bookkeeping
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 2 — Physical constants and unit system (km, kg, s)');

G     = 6.67430e-20;          % km^3 kg^-1 s^-2
M_sun = 1.989e30;             % kg
AU    = 1.496e8;              % km
yr    = 365.25 * 24 * 3600;   % s

v_earth_check = sqrt(G * M_sun / AU);
fprintf('Sanity check — Earth circular speed: %.2f km/s  (expect ~29.8)\n', v_earth_check);

fprintf('\nReference orbital radii and circular speeds:\n');
fprintf('  %-10s  %8s  %10s\n', 'Body', 'r [AU]', 'v_c [km/s]');
fprintf('  %s\n', repmat('-',1,32));
bodies_ref = {'Venus','Earth','Mars','Jupiter'};
radii_ref  = [0.723, 1.000, 1.524, 5.203];
for k = 1:numel(bodies_ref)
    vc = sqrt(G * M_sun / (radii_ref(k) * AU));
    fprintf('  %-10s  %8.3f  %10.2f\n', bodies_ref{k}, radii_ref(k), vc);
end

hkPause(HK);

%% Section 3 — Solar-System-like initial conditions
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 3 — Solar-System-like initial conditions (Sun + 4 pseudo-planets)');

%{
  LECTURE CONNECTION
  ==================
  This section sets up the initial-value problem discussed in the notes.
  We keep the main lecture run intentionally simple:
    • one dominant central mass (the Sun)
    • four lighter pseudo-planets
    • random orbital phases
    • circular-orbit speeds about the Sun
    • all motion initially in the x-y plane

  This produces a recognizable Solar-System-like simulation and gives a clean
  setting for the conserved-quantity checks in Sections 6–9.
%}

rng(42);    % reproducible lecture run

n = 5;      % Sun + 4 pseudo-planets
body_names = {'Sun','P2','P3','P4','P5'};

% Orbital radii (AU): Venus, Earth, Mars, Jupiter style spacing
r_AU = [0, 0.72, 1.00, 1.52, 5.20];

% Random planet masses in a realistic planetary range
mass_min  = 1e24;
mass_max  = 1e27;
m_planets = mass_min + (mass_max - mass_min) * rand(1, n-1);
m_bodies  = [M_sun, m_planets];

% Preallocate state-vector pieces
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

fprintf('\nInitial conditions summary (n = %d bodies):\n', n);
fprintf('  %-6s  %-10s  %-30s  %-30s\n', 'Body', 'Mass [kg]', 'Position [km]', 'Velocity [km/s]');
fprintf('  %s\n', repmat('-',1,82));
for i = 1:n
    ri = r0_all((3*i-2):(3*i));
    vi = v0_all((3*i-2):(3*i));
    fprintf('  %-6s  %10.3e  [%+.3e %+.3e %+.3e]  [%+.3e %+.3e %+.3e]\n', ...
        body_names{i}, m_bodies(i), ri(1),ri(2),ri(3), vi(1),vi(2),vi(3));
end

hkPause(HK);

%% Section 4 — Integrate with ode45
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 4 — Integrating the N-body equations with ode45');

T_end = 5 * yr;
tspan = linspace(0, T_end, 5000);
options = odeset('AbsTol', 1e-13, 'RelTol', 1e-11);

fprintf('Integrating %d bodies over %.1f years (%d output times)...\n', ...
    n, T_end/yr, numel(tspan));

tic
[T, X] = ode45(@(t,x) eom_NBP(t, x, m_bodies, G), tspan, x0, options);
elapsed = toc;

fprintf('Integration complete in %.2f s\n', elapsed);
fprintf('Output size: %d × %d  (time steps × state dimension)\n', size(X,1), size(X,2));

% Split the state history into positions and velocities
r_all = X(:, 1:3*n);
v_all = X(:, 3*n+1:end);

hkPause(HK);

%% Section 5 — 3-D trajectory plot
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 5 — 3-D trajectory plot: Solar-System-like motion and the barycenter');

%{
  WHAT TO LOOK FOR
  ================
  • The Sun should move very little compared with the planets.
  • The planets should remain in roughly Keplerian-looking trajectories.
  • The dashed black curve is the system center of mass (barycenter).
  • The legend identifies trajectories only; start/end markers are omitted.
%}

M_total = sum(m_bodies);
cm_all  = centerOfMassHistory(r_all, m_bodies);

if doPlot
    body_colors = [0.95 0.85 0.10;   % Sun: gold
                   0.60 0.80 0.95;   % P2 : light blue
                   0.20 0.70 0.30;   % P3 : green
                   0.85 0.40 0.20;   % P4 : rust-red
                   0.70 0.50 0.85];  % P5 : purple

    figure('Color','w','Name','MAE 240 — Solar-System-like trajectories');
    hold on; grid on; box on; axis equal;

    h_traj = gobjects(n+1,1);
    for i = 1:n
        ri = r_all(:, (3*i-2):(3*i));

        h_traj(i) = plot3(ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, ...
            '-', 'LineWidth', 1.5, 'Color', body_colors(i,:), 'DisplayName', body_names{i});

        plot3(ri(1,1)/AU, ri(1,2)/AU, ri(1,3)/AU, ...
            'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', body_colors(i,:), ...
            'HandleVisibility', 'off');

        plot3(ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, ...
            's', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', body_colors(i,:), ...
            'HandleVisibility', 'off');
    end

    h_traj(end) = plot3(cm_all(:,1)/AU, cm_all(:,2)/AU, cm_all(:,3)/AU, ...
        'k--', 'LineWidth', 1.75, 'DisplayName', 'Center of Mass');

    xlabel('x  [AU]', 'FontSize', 14);
    ylabel('y  [AU]', 'FontSize', 14);
    zlabel('z  [AU]', 'FontSize', 14);
    title(sprintf('N-body trajectories over %.0f years  (n = %d)', T_end/yr, n), ...
        'FontSize', 14);
    legend(h_traj, 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 13);
    view(3);

    annotation('textbox', [0.13 0.78 0.18 0.08], 'String', ...
        {'o = start', 's = finish'}, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w');
end

hkPause(HK);

%% Section 6 — First integral: center of mass / barycenter
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 6 — Barycenter from linear momentum: R_c(t) = R_c(t_0) + V_c(t_0)(t-t_0)');

%{
  LECTURE CONNECTION
  ==================
  The lecture notes define the barycenter by

      R_c = (1/M) Σ m_i r_i,   V_c = (1/M) Σ m_i v_i

  and show that conservation of linear momentum implies

      R_c(t) = R_c(t_0) + V_c(t_0) (t - t_0).
%}

p_all  = linearMomentumHistory(v_all, m_bodies);
Vcm0   = p_all(1,:) / M_total;
cm_ref = cm_all(1,:) + T(:) * Vcm0;
cm_err = cm_all - cm_ref;
cm_err_norm = sqrt(sum(cm_err.^2, 2));

fprintf('Center of mass at t=0   (AU): [%+.4e  %+.4e  %+.4e]\n', ...
    cm_all(1,1)/AU, cm_all(1,2)/AU, cm_all(1,3)/AU);
fprintf('Center of mass at t=T   (AU): [%+.4e  %+.4e  %+.4e]\n', ...
    cm_all(end,1)/AU, cm_all(end,2)/AU, cm_all(end,3)/AU);
fprintf('Max deviation from R_c(t_0)+V_c(t_0)(t-t_0): %.4e km  (%.4e AU)\n', ...
    max(cm_err_norm), max(cm_err_norm)/AU);

if doPlot
    figure('Color','w','Name','MAE 240 — Barycenter');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    labels = {'x','y','z'};
    for k = 1:3
        nexttile;
        plot(T/yr, cm_all(:,k)/AU, 'k-', 'LineWidth', 1.75); hold on;
        plot(T/yr, cm_ref(:,k)/AU, 'r--', 'LineWidth', 1.20);
        grid on; box on;
        ylabel(sprintf('R_{c,%s}  [AU]', labels{k}), 'FontSize', 12);
        if k == 1
            title('Center of mass: numerical history vs. linear prediction from the notes', 'FontSize', 13);
            legend('numerical', 'R_c(t_0)+V_c(t_0)(t-t_0)', 'Location', 'best');
        end
        if k == 3
            xlabel('Time  [yr]', 'FontSize', 12);
        end
    end
    title(tl, 'Barycenter / center-of-mass motion', 'FontSize', 13);
end

hkPause(HK);

%% Section 7 — First integral: linear momentum
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 7 — First integral: total linear momentum is constant');

%{
  LECTURE CONNECTION
  ==================
      P = Σ m_i v_i = const.

  Because the Solar-System-like ICs are not initialized in the barycentric
  frame, P(0) is generally not zero. That is fine: the conserved quantity is
  P itself, not necessarily the zero vector.
%}

dp_all = p_all - p_all(1,:);
max_dp = max(abs(dp_all), [], 1);
relP   = zeros(size(p_all));

for k = 1:3
    p0k = p_all(1,k);
    if abs(p0k) > eps
        relP(:,k) = p_all(:,k) / p0k;
    else
        relP(:,k) = NaN(size(T));
    end
end

fprintf('Initial total linear momentum P(0) [kg km/s]:\n');
fprintf('  p_x = %+.4e,  p_y = %+.4e,  p_z = %+.4e\n', p_all(1,1), p_all(1,2), p_all(1,3));
fprintf('Max absolute drift in the components [kg km/s]:\n');
fprintf('  Δp_x = %.4e,  Δp_y = %.4e,  Δp_z = %.4e\n', max_dp(1), max_dp(2), max_dp(3));

if doPlot
    figure('Color','w','Name','MAE 240 — Linear momentum');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    labels = {'x','y','z'};
    for k = 1:3
        nexttile;
        p0k = p_all(1,k);
        if abs(p0k) > eps
            plot(T/yr, relP(:,k), 'k-', 'LineWidth', 1.75);
            ylabel(sprintf('p_%s / p_%s(t_0)', labels{k}, labels{k}), 'FontSize', 12);
        else
            plot(T/yr, p_all(:,k), 'k-', 'LineWidth', 1.75);
            ylabel(sprintf('p_%s  [kg km/s]', labels{k}), 'FontSize', 12);
        end
        grid on; box on;
        if k == 1
            title('Linear momentum components should remain constant', 'FontSize', 13);
        end
        if k == 3
            xlabel('Time  [yr]', 'FontSize', 12);
        end
    end
end

hkPause(HK);

%% Section 8 — First integral: angular momentum
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 8 — First integral: total angular momentum is constant');

%{
  LECTURE CONNECTION
  ==================
      H = Σ r_i × (m_i v_i) = const.

  For this planar Solar-System-like run, H_x and H_y should stay near zero,
  while H_z carries the dominant angular momentum. Reporting only an absolute
  ΔH_z can look misleadingly large because H_z itself is enormous, so we also
  report relative drift measures.
%}

H_all = angularMomentumHistory(r_all, v_all, m_bodies);
dH_all = H_all - H_all(1,:);
dH_abs = max(abs(dH_all), [], 1);
H0_norm = sqrt(sum(H_all(1,:).^2));
rel_H_norm = sqrt(sum(dH_all.^2, 2)) / max(H0_norm, eps);
rel_Hz = abs(H_all(:,3) - H_all(1,3)) / max(abs(H_all(1,3)), eps);

fprintf('Initial angular momentum H(0) [kg km^2/s]:\n');
fprintf('  H_x = %+.4e,  H_y = %+.4e,  H_z = %+.4e\n', H_all(1,1), H_all(1,2), H_all(1,3));
fprintf('Max absolute drift in the components [kg km^2/s]:\n');
fprintf('  ΔH_x = %.4e,  ΔH_y = %.4e,  ΔH_z = %.4e\n', dH_abs(1), dH_abs(2), dH_abs(3));
fprintf('Max relative drift in total angular-momentum norm: %.4e\n', max(rel_H_norm));
fprintf('Max relative drift in H_z: %.4e\n', max(rel_Hz));

if doPlot
    figure('Color','w','Name','MAE 240 — Angular momentum');
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(T/yr, H_all(:,1), 'k-', 'LineWidth', 1.75);
    grid on; box on;
    ylabel('H_x  [kg km^2/s]', 'FontSize', 12);
    title('Angular momentum: H_x and H_y stay near zero; H_z stays constant', 'FontSize', 13);

    nexttile;
    plot(T/yr, H_all(:,2), 'k-', 'LineWidth', 1.75);
    grid on; box on;
    ylabel('H_y  [kg km^2/s]', 'FontSize', 12);

    nexttile;
    plot(T/yr, H_all(:,3) / H_all(1,3), 'k-', 'LineWidth', 1.75);
    grid on; box on;
    ylabel('H_z / H_z(t_0)', 'FontSize', 12);
    xlabel('Time  [yr]', 'FontSize', 12);
end

hkPause(HK);

%% Section 9 — First integral: total energy
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 9 — First integral: total mechanical energy is constant');

%{
  LECTURE CONNECTION
  ==================
  In the lecture notes, the positive force potential U is introduced so that

      E = T - U.

  Here we compute the equivalent conventional gravitational potential energy

      PE = -U

  and therefore use

      E = KE + PE.
%}

[KE_all, PE_all, E_all] = totalEnergyHistory(r_all, v_all, m_bodies, G);
rel_err_E = abs(E_all - E_all(1)) / abs(E_all(1));

fprintf('E(0)   = %.6e  [kg km^2/s^2]\n', E_all(1));
fprintf('KE(0)  = %.6e  [kg km^2/s^2]\n', KE_all(1));
fprintf('PE(0)  = %.6e  [kg km^2/s^2]\n', PE_all(1));
fprintf('E is bound (negative): %s\n', string(E_all(1) < 0));
fprintf('Max relative energy error: %.4e\n', max(rel_err_E));

if doPlot
    figure('Color','w','Name','MAE 240 — Energy conservation');
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(T/yr, KE_all/abs(E_all(1)), 'b-', 'LineWidth', 1.75); hold on;
    plot(T/yr, PE_all/abs(E_all(1)), 'r-', 'LineWidth', 1.75);
    plot(T/yr, E_all/abs(E_all(1)),  'k-', 'LineWidth', 2.00);
    grid on; box on;
    ylabel('Energy / |E(t_0)|', 'FontSize', 12);
    legend('KE','PE','E = KE+PE', 'Location', 'best', 'FontSize', 11);
    title('Kinetic, potential, and total energy', 'FontSize', 13);

    nexttile;
    semilogy(T/yr, rel_err_E + eps, 'k-', 'LineWidth', 1.75);
    grid on; box on;
    ylabel('|ΔE/E_0|', 'FontSize', 12);
    xlabel('Time  [yr]', 'FontSize', 12);
    title('Relative energy error', 'FontSize', 13);
end

hkPause(HK);

%% Section 10 — Randomized compact ICs: sensitivity and chaos preview
if ~exist('HK','var') || ~isstruct(HK), HK = defaultHK(); end
if ~exist('fmtTitle','var') || ~isa(fmtTitle,'function_handle'), fmtTitle = defaultFmtTitle(); end
if HK.enabled, hkSection(HK); end
fmtTitle('Section 10 — Randomized compact ICs: sensitivity to initial conditions');

%{
  LECTURE CONNECTION
  ==================
  To make the final figure more dynamically interesting, we do NOT use the
  original large random cube with large random velocities. That choice often
  sends the bodies off on nearly straight lines.

  Instead, we generate a randomized compact cluster:
    • comparable random masses
    • random 3-D positions in a small sphere
    • random 3-D velocities with the net momentum removed
    • velocity rescaling until the total energy is negative

  The result is still randomized, but much more likely to show strong mutual
  interaction over the integration window.
%}

n_rand = 5;
T_rand     = 2 * yr;
tspan_rand = linspace(0, T_rand, 3000);
opt_fast   = odeset('AbsTol',1e-8,'RelTol',1e-6);   % loose tol for seed search only

% Search for a seed that ejects a body other than Body 1, for variety in outcomes.
fprintf('Randomized compact IC run: n = %d bodies — searching for diverse-outcome seed...\n', n_rand);
seed_use   = 7;     % fallback
far_body   = 1;
for s_try = [7 13 17 23 31 37 41 47 53 59 61 67 71 79 83 89 97]
    rng(s_try);
    [m_s, r0_s, v0_s] = makeRandomCompactICs(n_rand, G, AU, mass_min, mass_max);
    [~, X_s] = ode45(@(t,x) eom_NBP(t, x, m_s, G), tspan_rand, [r0_s, v0_s].', opt_fast);
    r_end = X_s(end, 1:3*n_rand);
    cm_s  = zeros(1,3);
    for ii = 1:n_rand
        cm_s = cm_s + m_s(ii) * r_end((3*ii-2):(3*ii));
    end
    cm_s = cm_s / sum(m_s);
    d_end = zeros(1,n_rand);
    for ii = 1:n_rand
        d_end(ii) = norm(r_end((3*ii-2):(3*ii)) - cm_s);
    end
    [~, far_body_try] = max(d_end);
    if far_body_try ~= 1
        seed_use = s_try;
        far_body = far_body_try;
        break;
    end
    seed_use = s_try;   % keep last tried if all eject Body 1
    far_body = far_body_try;
end
fprintf('Seed %d selected — most distant body at t_end: Body %d\n', seed_use, far_body);

rng(seed_use);
[m_rand, r0_rand, v0_rand] = makeRandomCompactICs(n_rand, G, AU, mass_min, mass_max);
x0_rand = [r0_rand, v0_rand].';

tic
[T_r, X_r] = ode45(@(t,x) eom_NBP(t, x, m_rand, G), tspan_rand, x0_rand, options);
elapsed_rand = toc;
fprintf('Integration complete in %.2f s\n', elapsed_rand);

r_rand = X_r(:, 1:3*n_rand);
v_rand = X_r(:, 3*n_rand+1:end);
cm_rand = centerOfMassHistory(r_rand, m_rand);
[KE_r, PE_r, E_r] = totalEnergyHistory(r_rand, v_rand, m_rand, G);
rel_err_r = abs(E_r - E_r(1)) / abs(E_r(1));

fprintf('Random IC energy: E(0) = %.4e  kg km^2/s^2  (bound if < 0: %s)\n', ...
    E_r(1), string(E_r(1) < 0));
fprintf('Random IC max relative energy error: %.4e\n', max(rel_err_r));

if doPlot
    rand_colors = lines(n_rand);

    figure('Color','w','Name','MAE 240 — Randomized compact trajectories');
    tl_rand = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(tl_rand, ...
        sprintf('Randomized compact ICs: %d bodies over %.0f years  (o = start,  s = finish)', ...
        n_rand, T_rand/yr), 'FontSize', 13);

    ax3d = nexttile; hold(ax3d,'on'); grid(ax3d,'on'); box(ax3d,'on');
    axXY = nexttile; hold(axXY,'on'); grid(axXY,'on'); box(axXY,'on'); axis(axXY,'equal');
    axXZ = nexttile; hold(axXZ,'on'); grid(axXZ,'on'); box(axXZ,'on'); axis(axXZ,'equal');
    axYZ = nexttile; hold(axYZ,'on'); grid(axYZ,'on'); box(axYZ,'on'); axis(axYZ,'equal');

    h_rand = gobjects(n_rand,1);
    for i = 1:n_rand
        ri = r_rand(:, (3*i-2):(3*i));
        c  = rand_colors(i,:);
        nm = sprintf('Body %d', i);

        % 3-D view
        h_rand(i) = plot3(ax3d, ri(:,1)/AU, ri(:,2)/AU, ri(:,3)/AU, ...
            '-', 'LineWidth', 1.5, 'Color', c, 'DisplayName', nm);
        plot3(ax3d, ri(1,1)/AU,   ri(1,2)/AU,   ri(1,3)/AU,   ...
            'o','MarkerSize',7,'LineWidth',1.5,'Color',c,'HandleVisibility','off');
        plot3(ax3d, ri(end,1)/AU, ri(end,2)/AU, ri(end,3)/AU, ...
            's','MarkerSize',7,'LineWidth',1.5,'Color',c,'HandleVisibility','off');

        % x-y projection
        plot(axXY, ri(:,1)/AU, ri(:,2)/AU, '-','LineWidth',1.5,'Color',c,'HandleVisibility','off');
        plot(axXY, ri(1,1)/AU,   ri(1,2)/AU,   'o','MarkerSize',6,'Color',c,'HandleVisibility','off');
        plot(axXY, ri(end,1)/AU, ri(end,2)/AU, 's','MarkerSize',6,'Color',c,'HandleVisibility','off');

        % x-z projection
        plot(axXZ, ri(:,1)/AU, ri(:,3)/AU, '-','LineWidth',1.5,'Color',c,'HandleVisibility','off');
        plot(axXZ, ri(1,1)/AU,   ri(1,3)/AU,   'o','MarkerSize',6,'Color',c,'HandleVisibility','off');
        plot(axXZ, ri(end,1)/AU, ri(end,3)/AU, 's','MarkerSize',6,'Color',c,'HandleVisibility','off');

        % y-z projection
        plot(axYZ, ri(:,2)/AU, ri(:,3)/AU, '-','LineWidth',1.5,'Color',c,'HandleVisibility','off');
        plot(axYZ, ri(1,2)/AU,   ri(1,3)/AU,   'o','MarkerSize',6,'Color',c,'HandleVisibility','off');
        plot(axYZ, ri(end,2)/AU, ri(end,3)/AU, 's','MarkerSize',6,'Color',c,'HandleVisibility','off');
    end

    % Center of mass on 3-D panel only (constant position; shown for reference)
    plot3(ax3d, cm_rand(:,1)/AU, cm_rand(:,2)/AU, cm_rand(:,3)/AU, ...
        'k--', 'LineWidth', 1.75, 'HandleVisibility','off');

    xlabel(ax3d,'x  [AU]','FontSize',12); ylabel(ax3d,'y  [AU]','FontSize',12);
    zlabel(ax3d,'z  [AU]','FontSize',12); title(ax3d,'3-D view','FontSize',12);
    view(ax3d,3); axis(ax3d,'equal'); set(ax3d,'FontSize',11);
    legend(ax3d, h_rand, 'Location','best','FontSize',11);

    xlabel(axXY,'x  [AU]','FontSize',12); ylabel(axXY,'y  [AU]','FontSize',12);
    title(axXY,'x-y  projection','FontSize',12); set(axXY,'FontSize',11);

    xlabel(axXZ,'x  [AU]','FontSize',12); ylabel(axXZ,'z  [AU]','FontSize',12);
    title(axXZ,'x-z  projection','FontSize',12); set(axXZ,'FontSize',11);

    xlabel(axYZ,'y  [AU]','FontSize',12); ylabel(axYZ,'z  [AU]','FontSize',12);
    title(axYZ,'y-z  projection','FontSize',12); set(axYZ,'FontSize',11);

    figure('Color','w','Name','MAE 240 — Energy: Solar-System-like vs Randomized');
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    semilogy(T/yr, rel_err_E + eps, 'k-', 'LineWidth', 1.75);
    grid on; box on;
    xlabel('Time  [yr]', 'FontSize', 12);
    ylabel('|ΔE/E_0|', 'FontSize', 12);
    title('Solar-System-like ICs', 'FontSize', 13);

    nexttile;
    semilogy(T_r/yr, rel_err_r + eps, 'r-', 'LineWidth', 1.75);
    grid on; box on;
    xlabel('Time  [yr]', 'FontSize', 12);
    ylabel('|ΔE/E_0|', 'FontSize', 12);
    title('Randomized compact ICs', 'FontSize', 13);

    title(tl, 'Relative energy error: integration quality comparison', 'FontSize', 13);
end

disp(' ');
disp('Wrap-up takeaways:');
disp('  • x = [r; v] is the 6n × 1 state vector for the Newtonian N-body IVP.');
disp('  • The Solar-System-like run is used to illustrate the conserved quantities.');
disp('  • The barycenter follows the linear relation implied by linear momentum conservation.');
disp('  • Linear momentum, angular momentum, and total energy remain nearly constant numerically.');
disp('  • The randomized compact run illustrates sensitivity to initial conditions without simple straight-line escape.');

hkPause(HK);

%% Local functions

function HK = defaultHK()
%DEFAULTHK Default housekeeping configuration.
HK.enabled         = true;
HK.doClc           = false;
HK.doClear         = false;
HK.doClose         = false;
HK.pauseAtSections = false;
HK.useKeyboard     = false;
HK.clcAfterSection = false;
HK.clearVarsKeep   = {'HK','fmtTitle','doPlot'};
end

function hkSection(HK)
%HKSECTION Apply housekeeping at the top of a section.
if HK.doClose
    close all
end
if HK.doClear
    keep = HK.clearVarsKeep;
    clearvars('-except', keep{:});
end
if HK.doClc
    clc
end
end

function hkPause(HK)
%HKPAUSE Optional pacing and workspace inspection at the end of a section.
if ~(isfield(HK,'pauseAtSections') && HK.pauseAtSections)
    return
end
fprintf('\n');
if isfield(HK,'useKeyboard') && HK.useKeyboard
    disp('--- INSPECT at K>> (dbcont to continue; dbquit to stop) ---');
    try
        w = evalin('caller','whos');
        for kk = 1:numel(w)
            name = w(kk).name;
            if isvarname(name)
                val = evalin('caller', name);
                eval([name ' = val;']); %#ok<EVLDIR>
            end
        end
        clear w kk name val
    catch
    end
    keyboard
else
    disp('--- press any key to continue ---');
    pause
end
if isfield(HK,'clcAfterSection') && HK.clcAfterSection
    clc
end
end

function fmtTitle = defaultFmtTitle()
%DEFAULTFMTTITLE Returns a function handle for consistent section title banners.
fmtTitle = @(s) localFmtTitle(s);
end

function localFmtTitle(s)
%LOCALFMTTITLE Print a consistent 78-char banner.
line = repmat('-', 1, 78);
fprintf('\n%s\n%s\n%s\n', line, s, line);
end

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

function [m_rand, r0_rand, v0_rand] = makeRandomCompactICs(n_rand, G, AU, mass_min, mass_max)
%MAKERANDOMCOMPACTICS Randomized compact cluster with negative total energy.

% Comparable masses; no dominant central body in the randomized case.
m_rand = mass_min + (mass_max - mass_min) * rand(1, n_rand);

% Random 3-D positions inside a compact sphere.
cluster_radius = 0.08 * AU;
r0_mat = zeros(n_rand,3);
for i = 1:n_rand
    placed = false;
    while ~placed
        trial = (2*rand(1,3) - 1) * cluster_radius;
        if norm(trial) <= cluster_radius
            r0_mat(i,:) = trial;
            placed = true;
        end
    end
end

% Random 3-D velocities with modest scale.
vel_scale = 3.5;   % km/s
v0_mat = (2*rand(n_rand,3) - 1) * vel_scale;

% Remove net linear momentum so the cluster is initialized in its CoM frame.
P0 = [0, 0, 0];
M_total = sum(m_rand);
for i = 1:n_rand
    P0 = P0 + m_rand(i) * v0_mat(i,:);
end
Vcm0 = P0 / M_total;
for i = 1:n_rand
    v0_mat(i,:) = v0_mat(i,:) - Vcm0;
end

% Rescale velocities until the system is bound: E = KE + PE < 0.
E0 = randomClusterEnergy(r0_mat, v0_mat, m_rand, G);
iter = 0;
while E0 >= 0 && iter < 25
    v0_mat = 0.75 * v0_mat;
    E0 = randomClusterEnergy(r0_mat, v0_mat, m_rand, G);
    iter = iter + 1;
end

% Pack into the same 1 × 3n layout used in the script.
r0_rand = reshape(r0_mat.', 1, []);
v0_rand = reshape(v0_mat.', 1, []);
end

function E0 = randomClusterEnergy(r0_mat, v0_mat, m_rand, G)
%RANDOMCLUSTERENERGY Total energy for a random compact cluster.

n_rand = numel(m_rand);
KE0 = 0;
for i = 1:n_rand
    KE0 = KE0 + 0.5 * m_rand(i) * dot(v0_mat(i,:), v0_mat(i,:));
end

PE0 = 0;
for i = 1:n_rand
    for j = i+1:n_rand
        d_ij = norm(r0_mat(j,:) - r0_mat(i,:));
        PE0 = PE0 - G * m_rand(i) * m_rand(j) / d_ij;
    end
end

E0 = KE0 + PE0;
end
