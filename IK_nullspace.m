%% Ik_nullspace.m  --  Franka Emika Panda IK Null-Space Motion
%
% Sweeps the redundancy parameter q7 over a range and computes the resulting
% joint angles via the analytical IK (A2-branch, flange frame).
% Identifies which configurations satisfy the hardware joint limits.
%
% Frame convention : Flange frame  (dee=0, no pi/4 offset on q7).
%                    Corresponds to forward_kinematics_modDH.m.
%                    A2-branch only: q4 = phi1+phi2+phi3-2*pi.
%
% Inputs  (set directly in script):
%   p_ee       -- (3x1) TCP position [m]
%   roll/pitch/yaw -- ZYX Euler angles [rad]
%   qa         -- (7x1) reference configuration [rad]
%
% Outputs:
%   nullspace_traj -- (Nx7) joint angles [rad] for each q7 step
%   is_valid       -- (Nx1) logical: true if configuration is within joint limits
%   Figure with two subplots (joint angles + validity over q7 sweep)
%
% Reference:
%   He & Liu, "Analytical Inverse Kinematics for Franka Emika Panda --
%   a Geometrical Solver for the 7-DOF Robot Arm", 2021.
%clear; clc; close all;

%% 1. INPUT VALUES
p_ee = [-0.33938; 0.58864; 0.23309]; % TCP position
roll = deg2rad(-36.52); pitch = deg2rad(83.16); yaw = deg2rad(100.05); % TCP orientation
qa = deg2rad([118.68; 35.20; 0.00; -103.83; 158.32; 123.24; 10.00]); % Reference

%% 2. NULL-SPACE INITIALISATION
q7_range = deg2rad(linspace(-120, 120, 100)); % 100 steps for smooth trajectory
nullspace_traj = zeros(length(q7_range), 7);
is_valid = false(length(q7_range), 1);

fprintf('Computing null-space motion for %d steps...\n', length(q7_range));

for i = 1:length(q7_range)
    current_q7 = q7_range(i);

    % Compute all 4 solutions for this q7
    sols = computeAllSolutions(p_ee, roll, pitch, yaw, current_q7);

    if ~isempty(sols)
        % Select the solution closest to reference qa
        % (prevents "jumping" between solutions during the sweep)
        diffs = sum((sols - qa').^2, 2);
        [~, best_idx] = min(diffs);

        nullspace_traj(i, :) = sols(best_idx, :);

        % Check joint limits for this solution
        is_valid(i) = check_franka_limits(nullspace_traj(i, :));
    end
end

%% 3. VISUALISATION
% Commented-out first draft kept for reference:
% figure('Color', 'w', 'Position', [100, 100, 800, 600]);
%
% subplot(2,1,1);
% plot(rad2deg(q7_range), rad2deg(nullspace_traj(:, 1:6)), 'LineWidth', 1.5);
% grid on; hold on;
% title('Null-space motion: joint angles q1-q6 over q7');
% xlabel('Redundancy angle q7 [deg]'); ylabel('Angle [deg]');
% legend('q1','q2','q3','q4','q5','q6', 'Location', 'eastoutside');
%
% subplot(2,1,2);
% area(rad2deg(q7_range), is_valid, 'FaceColor', [0.8 1 0.8], 'EdgeColor', 'g');
% ylim([-0.1 1.1]); grid on;
% title('Configuration validity (within hardware limits)');
% xlabel('q7 [deg]'); ylabel('Valid (1=Yes)');
% yticklabels({'No', 'Yes'});

figure('Color', 'w', 'Position', [100, 100, 900, 600]);

subplot(2,1,1);
plot(rad2deg(q7_range), rad2deg(nullspace_traj(:, 1:6)), 'LineWidth', 1.5);
grid on; hold on;
title('Null-space motion: joint angles q1-q6 over q7');
xlabel('Redundancy angle q7 [deg]'); ylabel('Angle [deg]');
legend('q1','q2','q3','q4','q5','q6', 'Location', 'eastoutside');

subplot(2,1,2);
h_area = area(rad2deg(q7_range), is_valid, 'FaceColor', [0.8 1 0.8], 'EdgeColor', [0 0.5 0]);
ylim([-0.1 1.1]); grid on;
title('Configuration validity (within hardware limits)');
xlabel('q7 [deg]'); ylabel('Valid (1=Yes)');
set(gca, 'YTick', [0 1], 'YTickLabel', {'No', 'Yes'});
legend(h_area, 'Valid', 'Location', 'eastoutside');


linkaxes(findall(gcf,'type','axes'), 'x');

%% --- HELPER FUNCTIONS ---

% computeAllSolutions -- compute all 4 IK solutions (A2 x B1/B2 x C1/C2) for a given q7.
%   p_ee  -- (3x1) TCP position [m]
%   roll, pitch, yaw -- ZYX Euler angles [rad]
%   q7    -- redundancy parameter [rad]
%   Returns ALL_SOLS: (up to 4x7) matrix of joint angles [rad], or [] if unreachable.
function ALL_SOLS = computeAllSolutions(p_ee, roll, pitch, yaw, q7)
    % Constants & DH parameters
    d1 = 0.333; d3 = 0.316; a4 = 0.0825; d5 = 0.384; a5 = -0.0825; a7 = 0.088; df = 0.107;
    L24 = sqrt(d3^2 + a4^2); L46 = sqrt(d5^2 + a5^2);

    rot_x = @(phi) [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    rot_y = @(phi) [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
    rot_z = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    R_ee = rot_z(yaw) * rot_y(pitch) * rot_x(roll);
    z_ee = R_ee(:, 3);
    p6 = (p_ee - df * z_ee) - a7 * (R_ee * [cos(-q7); sin(-q7); 0]);
    p2 = [0; 0; d1];
    dist_O2O6 = norm(p6 - p2);

    % Reachability check
    cos_phi3 = (L24^2 + L46^2 - dist_O2O6^2) / (2 * L24 * L46);
    if abs(cos_phi3) > 1, ALL_SOLS = []; return; end

    phi1 = atan2(d3, a4); phi2 = atan2(d5, abs(a5)); phi3 = acos(cos_phi3);
    y6 = -R_ee(:,3); x6 = R_ee * [cos(-q7); sin(-q7); 0]; z6 = cross(x6, y6); R6 = [x6, y6, z6];
    phi_O2O6O4 = acos(clamp((dist_O2O6^2 + L46^2 - L24^2) / (2 * dist_O2O6 * L46), -1, 1));
    angle_O2O6H = phi_O2O6O4 + atan(abs(a5) / d5);
    vec_6_O2O6 = R6' * (p6 - p2);
    phi6 = atan2(vec_6_O2O6(2), vec_6_O2O6(1));
    psi6 = asin(clamp((dist_O2O6 * cos(angle_O2O6H)) / norm(vec_6_O2O6(1:2)), -1, 1));

    ALL_SOLS = [];
    % A2-branch only: phi1+phi2+phi3-2*pi  (paper Eq. (48))
    cq4 = phi1 + phi2 + phi3 - 2*pi;
    q6_opts = [(pi - psi6 - phi6), (psi6 - phi6)]; % B1, B2
    for s6 = 1:2
            cq6 = atan2(sin(q6_opts(s6)), cos(q6_opts(s6)));
            phi_O3O2O4 = atan(a4 / d3);
            phi_O4O2O6 = acos(clamp((L24^2 + dist_O2O6^2 - L46^2) / (2 * L24 * dist_O2O6), -1, 1));

            % Denominator check for dist_PO6 (singularity protection)
            denominator = sin(phi_O2O6O4 + phi3 + phi1 - angle_O2O6H - pi/2);
            if abs(denominator) < 1e-6, denominator = 1e-6; end
            dist_PO6 = dist_O2O6 * sin(phi_O3O2O4 + phi_O4O2O6) / denominator;

            vec_O2P = (p6 - p2) - dist_PO6 * (R6 * [sin(cq6); cos(cq6); 0]);
            q12 = [atan2(vec_O2P(2), vec_O2P(1)), acos(clamp(vec_O2P(3)/norm(vec_O2P), -1, 1));
                   atan2(-vec_O2P(2), -vec_O2P(1)), -acos(clamp(vec_O2P(3)/norm(vec_O2P), -1, 1))];

            for s12 = 1:2
                cq1 = q12(s12,1); cq2 = q12(s12,2);
                R2 = rot_z(cq1) * rot_x(-pi/2) * rot_z(cq2);
                z3 = vec_O2P / norm(vec_O2P);
                % Cross product for y3 must be numerically stable
                cp = cross(vec_O2P, (p6 - p2));
                if norm(cp) < 1e-9, y3 = [0;0;1]; % fallback
                else, y3 = cp / norm(cp); end
                x3 = cross(y3, z3);

                x3_in_2 = R2' * x3;
                cq3 = atan2(x3_in_2(3), x3_in_2(1));

                % q5 computation
                z5 = R6 * [sin(cq6); cos(cq6); 0];
                vec_5_HO4 = (rot_x(pi/2)*rot_z(cq6)) * R6' * (p2 + d3*z3 + a4*x3 - (p6 - d5*z5));
                cq5 = -atan2(vec_5_HO4(2), vec_5_HO4(1));

                ALL_SOLS = [ALL_SOLS; cq1, cq2, cq3, cq4, cq5, cq6, q7];
            end
        end
end

function out = clamp(in, low, high)
    out = max(min(in, high), low);
end

function ok = check_franka_limits(q_rad)
    q_deg = rad2deg(q_rad);
    lims = [-166,166; -101,101; -166,166; -176,-4; -166,166; -1,215; -166,166];
    ok = all(q_deg >= lims(:,1)' & q_deg <= lims(:,2)');
end
