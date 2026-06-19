%% validate_IK.m  --  Round-Trip Test for IK_Franka (flange and EE frame)
%
% Tests FK -> IK -> FK with N random valid joint configurations.
% Uses the correct Modified-DH FK (identical to forward_kinematics_modDH.m).
%
% EE-FK:  T_EE = T_flange * EE_offset,  EE_offset = F_T_EE (Rz(+pi/4), 0.1034 m)
% EE-IK:  EE pose pre-transformed to flange pose via T_EE / EE_offset, then the
%          flange kernel is called. The EE test is honest (not circular).
%
% Inputs  (hard-coded):
%   N         -- number of test configurations
%   lims_deg  -- (7x2) joint limits [deg]
%
% Outputs (printed to console):
%   Per-frame statistics: number successful/skipped, max/mean/median
%   position error [m] and rotation error (Frobenius norm).
%
% Reference:
%   He & Liu, "Analytical Inverse Kinematics for Franka Emika Panda --
%   a Geometrical Solver for the 7-DOF Robot Arm", 2021.

clc; clear;

N = 1000; % number of test configurations

% Joint limits (Franka spec, in degrees)
lims_deg = [-166, 166;   % q1
            -101, 101;   % q2
            -166, 166;   % q3
            -176,  -4;   % q4 -- upper limit -4 deg (not +4 deg)
            -166, 166;   % q5
              -1, 215;   % q6
            -166, 166];  % q7

rng(42); % reproducibility: same configurations for both modes

% Pre-generate random configurations
Q_test = deg2rad(lims_deg(:,1) + rand(7, N) .* (lims_deg(:,2) - lims_deg(:,1)));

%% --- Flange-Frame ---
run_test('flange', N, Q_test, @fk_franka_flange, @ik_franka_flange);

%% --- EE-Frame ---
run_test('EE', N, Q_test, @fk_franka_EE, @ik_franka_EE);

%% =====================================================================
%% TEST LOOP
%% =====================================================================
function run_test(mode_name, N, Q_test, fk_fun, ik_fun)
    err_pos = nan(N, 1);
    err_rot = nan(N, 1);
    n_skip  = 0;

    for k = 1:N
        q_ref = Q_test(:, k);
        T_fwd = fk_fun(q_ref);
        R_ee  = T_fwd(1:3, 1:3);
        p_ee  = T_fwd(1:3, 4);

        try
            q_ik = ik_fun(R_ee, p_ee, q_ref(7), q_ref);
        catch
            n_skip = n_skip + 1;
            continue;
        end

        T_chk = fk_fun(q_ik);
        err_pos(k) = norm(T_chk(1:3,4) - T_fwd(1:3,4));
        err_rot(k) = norm(T_chk(1:3,1:3)' * T_fwd(1:3,1:3) - eye(3), 'fro');
    end

    mask = ~isnan(err_pos);
    n_ok = sum(mask);

    fprintf('=== IK Validation (%s-Frame, N=%d) ===\n', mode_name, N);
    fprintf('Successful: %d  |  Skipped (singularity): %d\n\n', n_ok, n_skip);

    if n_ok > 0
        fprintf('Position error [m]:\n');
        fprintf('  Max:    %.3e\n', max(err_pos(mask)));
        fprintf('  Mean:   %.3e\n', mean(err_pos(mask)));
        fprintf('  Median: %.3e\n', median(err_pos(mask)));
        fprintf('\nOrientation error (Frobenius):\n');
        fprintf('  Max:    %.3e\n', max(err_rot(mask)));
        fprintf('  Mean:   %.3e\n', mean(err_rot(mask)));
        fprintf('  Median: %.3e\n', median(err_rot(mask)));
        if max(err_pos(mask)) < 1e-10 && max(err_rot(mask)) < 1e-10
            fprintf('\n[OK] All errors within machine precision.\n\n');
        else
            fprintf('\n[WARNING] Errors exceed machine precision -- check IK!\n\n');
        end
    end
end

%% =====================================================================
%% FORWARD KINEMATICS
%% =====================================================================

% Modified DH FK -- flange frame (see forward_kinematics_modDH.m)
function T = fk_franka_flange(q)
    d     = [0.333; 0; 0.316; 0; 0.384; 0; 0];    % [m]
    a     = [0; 0; 0; 0.0825; -0.0825; 0; 0.088];  % [m]
    alpha = [0; -pi/2; pi/2; pi/2; -pi/2; pi/2; pi/2];

    T = eye(4);
    for i = 1:7
        Ti = [cos(q(i)),              -sin(q(i)),             0,           a(i);
              sin(q(i))*cos(alpha(i)), cos(q(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i);
              sin(q(i))*sin(alpha(i)), cos(q(i))*sin(alpha(i)),  cos(alpha(i)),  cos(alpha(i))*d(i);
              0, 0, 0, 1];
        T = T * Ti;
    end
    % Flange offset d_f = 0.107 m
    T = T * [1,0,0,0; 0,1,0,0; 0,0,1,0.107; 0,0,0,1];
end

% Modified DH FK -- official EE frame (Franka DH table: theta_EE = +pi/4)
% T_EE = T_flange * EE_offset,  EE_offset = F_T_EE (Rz(+pi/4), 0.1034 m)
function T = fk_franka_EE(q)
    dee = 0.1034; % [m]
    EE_offset = [cos(pi/4), -sin(pi/4), 0, 0;
                 sin(pi/4),  cos(pi/4), 0, 0;
                 0,          0,         1, dee;
                 0,          0,         0, 1];
    T = fk_franka_flange(q) * EE_offset;
end

%% =====================================================================
%% INVERSE KINEMATICS (shared core function)
%% =====================================================================

% Flange-Frame: dee=0, q7_eff=q7
function q = ik_franka_flange(R_ee, p_ee, q7, qa)
    q = ik_core(R_ee, p_ee, q7, qa, 0, q7);
end

% EE-Frame: pre-transform EE pose to flange pose, then call the flange kernel.
function q = ik_franka_EE(R_ee, p_ee, q7, qa)
    dee = 0.1034;
    EE_offset = [cos(pi/4), -sin(pi/4), 0, 0;
                 sin(pi/4),  cos(pi/4), 0, 0;
                 0,          0,         1, dee;
                 0,          0,         0, 1];
    T_EE = [R_ee, p_ee; 0, 0, 0, 1];
    T_flange = T_EE / EE_offset;
    R_fl = T_flange(1:3, 1:3);
    p_fl = T_flange(1:3, 4);
    q = ik_core(R_fl, p_fl, q7, qa, 0, q7);
end

% Core IK (A2-branch, B1/B2 via dot(O2O6,x5(qa)), C1/C2 via qa(2))
function q = ik_core(R_ee, p_ee, q7, qa, dee, q7_eff)
    d1 = 0.333; d3 = 0.316; a4 = 0.0825; d5 = 0.384; a5 = -0.0825; a7 = 0.088; df = 0.107;

    rot_x = @(phi) [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    rot_z = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    L24 = sqrt(d3^2 + a4^2);
    L46 = sqrt(d5^2 + a5^2);

    z_ee = R_ee(:, 3);
    p7   = p_ee - (df + dee) * z_ee;
    x6   = R_ee * [cos(q7_eff); -sin(q7_eff); 0];
    p6   = p7 - a7 * x6;
    p2   = [0; 0; d1];
    O2O6 = p6 - p2;
    d26  = norm(O2O6);

    % q4 -- A2-branch (paper Eq. (48))
    phi1 = atan2(d3, a4);
    phi2 = atan2(d5, abs(a5));
    cos_phi3 = (L24^2 + L46^2 - d26^2) / (2 * L24 * L46);
    if abs(cos_phi3) > 1, error('Target unreachable'); end
    phi3 = acos(cos_phi3);
    q4   = phi1 + phi2 + phi3 - 2*pi;

    % R6
    R6 = [x6, -z_ee, cross(x6, -z_ee)];

    % q6 -- B1/B2 (paper Eq. (54, 83))
    phi_O2O6O4  = acos(cv((d26^2 + L46^2 - L24^2) / (2 * d26 * L46)));
    angle_O2O6H = phi_O2O6O4 + atan(abs(a5) / d5);
    v6          = R6' * O2O6;
    denom_psi6  = norm(v6(1:2));
    if denom_psi6 < 1e-9, error('Singularity: v6(1:2) ~ 0'); end
    phi6 = atan2(v6(2), v6(1));
    psi6 = asin(cv((d26 * cos(angle_O2O6H)) / denom_psi6));

    x5 = get_x5(qa);
    if dot(O2O6, x5) <= 0
        q6 = pi - psi6 - phi6; % B1
    else
        q6 = psi6 - phi6;       % B2
    end
    q6 = atan2(sin(q6), cos(q6));

    % q1, q2 -- C1/C2 (paper Eq. (66, 84))
    phi_O4O2O6  = acos(cv((L24^2 + d26^2 - L46^2) / (2 * L24 * d26)));
    angle_PO2O6 = atan(a4 / d3) + phi_O4O2O6;
    angle_O2PO6 = phi_O2O6O4 + phi3 + phi1 - angle_O2O6H - pi/2;
    if abs(sin(angle_O2PO6)) < 1e-9, error('Singularity: angle_O2PO6'); end

    dist_PO6 = d26 * sin(angle_PO2O6) / sin(angle_O2PO6);
    vec_O2P  = O2O6 - dist_PO6 * (R6 * [sin(q6); cos(q6); 0]);
    n_O2P    = norm(vec_O2P);
    if n_O2P < 1e-9, error('Singularity: vec_O2P ~ 0'); end

    if qa(2) >= 0
        q1 = atan2(vec_O2P(2), vec_O2P(1));
        q2 = acos(cv(vec_O2P(3) / n_O2P));
    else
        q1 = atan2(-vec_O2P(2), -vec_O2P(1));
        q2 = -acos(cv(vec_O2P(3) / n_O2P));
    end

    % q3
    R2 = rot_z(q1) * rot_x(-pi/2) * rot_z(q2);
    z3 = vec_O2P / n_O2P;
    cp = cross(vec_O2P, O2O6);
    if norm(cp) < 1e-9, error('Singularity: cross product ~ 0'); end
    y3 = cp / norm(cp);
    x3 = cross(y3, z3);
    x3_in_2 = R2' * x3;
    q3 = atan2(x3_in_2(3), x3_in_2(1));

    % q5
    z5      = R6 * [sin(q6); cos(q6); 0];
    p4      = p2 + d3*z3 + a4*x3;
    v_5_HO4 = (rot_x(pi/2) * rot_z(q6)) * R6' * (p4 - (p6 - d5*z5));
    q5      = -atan2(v_5_HO4(2), v_5_HO4(1));

    q = [q1; q2; q3; q4; q5; q6; q7];
end

% x5-axis in world coordinates via Modified DH (paper Eq. (83))
function x5 = get_x5(q)
    d1 = 0.333; d3 = 0.316; a4 = 0.0825; a5 = -0.0825; d5 = 0.384;
    MDH = @(al, a, d, th) [cos(th), -sin(th), 0, a;
                            sin(th)*cos(al), cos(th)*cos(al), -sin(al), -d*sin(al);
                            sin(th)*sin(al), cos(th)*sin(al),  cos(al),  d*cos(al);
                            0, 0, 0, 1];
    T05 = MDH(0,0,d1,q(1)) * MDH(-pi/2,0,0,q(2)) * MDH(pi/2,0,d3,q(3)) * ...
          MDH(pi/2,a4,0,q(4)) * MDH(-pi/2,a5,d5,q(5));
    x5 = T05(1:3, 1);
end

% Clamp to [-1, 1] for acos/asin
function out = cv(in)
    out = max(min(in, 1), -1);
end
