%% IK_full.m  --  Franka Emika Panda IK: All 4 Solutions (A2 x B1/B2 x C1/C2)
%
% Computes all four closed-form IK solutions for the Franka Research 3 / Panda
% (7-DOF) using the geometric method of He & Liu (2021), A2-branch only.
%
% Frame convention : Flange frame  (dee=0, no pi/4 offset on q7).
%                    Corresponds to forward_kinematics_modDH.m.
%
% Inputs  (set directly in script):
%   p_ee       -- (3x1) TCP position [m]
%   roll/pitch/yaw -- ZYX Euler angles [rad]
%   qa         -- (7x1) reference configuration [rad]  (fixes q7; selects B branch)
%
% Outputs (printed to console):
%   Table of all 4 solutions [deg] with joint-limit validity flag.
%   Best valid solution (minimum joint-space distance to qa).
%
% Reference:
%   He & Liu, "Analytical Inverse Kinematics for Franka Emika Panda --
%   a Geometrical Solver for the 7-DOF Robot Arm", 2021.
%clear; clc;

%% 0. EXPLANATION
% q7: redundancy parameter (fixed to qa(7))
% q6: 2 solutions (Case B1 / B2) -- both stored in SOLS_RAD
% q1 and q2: 2 solutions (Case C1 / C2) -- both stored in SOLS_RAD
% A2-branch only: q4 = phi1+phi2+phi3-2*pi

%% 1. INPUT VALUES
% p_ee = [-0.17478; 0.42552; 0.23519]; % TCP position
% roll = deg2rad(-162.55); pitch = deg2rad(-14.53); yaw = deg2rad(-59.07); % TCP orientation
% qa = deg2rad([118.68; 35.20; 0.00; -103.83; 158.32; -123.24; 10.00]); % Reference

p_ee = [-0.33938; 0.58864; 0.23309]; % TCP position
roll = deg2rad(-36.52); pitch = deg2rad(83.16); yaw = deg2rad(100.05); % TCP orientation
qa = deg2rad([118.68;   35.20;    0.00; -103.83;  158.32; 123.24;   10.00]); % Reference for redundancy q7 and branch selection q6

%% 2. CONSTANTS & DH GEOMETRY
d1 = 0.333; d3 = 0.316; a4 = 0.0825; d5 = 0.384; a5 = -0.0825; a7 = 0.088;
df = 0.107; dee = 0;

L24 = sqrt(d3^2 + a4^2);
L46 = sqrt(d5^2 + a5^2);

% Helper functions
rot_x = @(phi) [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
rot_y = @(phi) [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
rot_z = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

%% 3. BASE CALCULATIONS
R_ee = rot_z(yaw) * rot_y(pitch) * rot_x(roll);
z_ee = R_ee(:, 3);
p7 = p_ee - (df + dee) * z_ee;

q7 = qa(7); % Redundancy parameter fixed
x6_ee = [cos(-q7); sin(-q7); 0];
x6 = R_ee * x6_ee;
p6 = p7 - a7 * x6;
p2 = [0; 0; d1];
dist_O2O6 = norm(p6 - p2);

% Geometry for q4
phi1 = atan2(d3, a4);
phi2 = atan2(d5, abs(a5));
phi3 = acos((L24^2 + L46^2 - dist_O2O6^2) / (2 * L24 * L46));

% Geometry for q6
y6 = -z_ee;
z6 = cross(x6, y6);
R6 = [x6, y6, z6];
phi_O2O6O4 = acos((dist_O2O6^2 + L46^2 - L24^2) / (2 * dist_O2O6 * L46));
phi_HO6O4 = atan(abs(a5) / d5);
angle_O2O6H = phi_O2O6O4 + phi_HO6O4;
vec_6_O2O6 = R6' * (p6 - p2);
phi6 = atan2(vec_6_O2O6(2), vec_6_O2O6(1));
psi6 = asin((dist_O2O6 * cos(angle_O2O6H)) / norm(vec_6_O2O6(1:2)));

%% 4. LOOP OVER ALL 4 SOLUTIONS (A2 x B1/B2 x C1/C2)
SOLS_RAD = [];

% q4 -- A2-branch only: phi1+phi2+phi3-2*pi  (paper Eq. (48))
% A1 (phi1+phi2-phi3-2*pi) yields q4 outside joint range -> no valid solution.
cq4 = phi1 + phi2 + phi3 - 2*pi;

% Branch 2: q6 (wrist, Case B1 / B2)
q6_options = [(pi - psi6 - phi6), (psi6 - phi6)];

for s6 = 1:2
        cq6 = atan2(sin(q6_options(s6)), cos(q6_options(s6)));

        % Shoulder preparation
        phi_O3O2O4 = atan(a4 / d3);
        phi_O4O2O6 = acos((L24^2 + dist_O2O6^2 - L46^2) / (2 * L24 * dist_O2O6));
        angle_PO2O6 = phi_O3O2O4 + phi_O4O2O6;
        angle_O2PO6 = phi_O2O6O4 + phi3 + phi1 - angle_O2O6H - pi/2;
        dist_PO6 = dist_O2O6 * sin(angle_PO2O6) / sin(angle_O2PO6);

        z5_in_6 = [sin(cq6); cos(cq6); 0];
        vec_O2P = (p6 - p2) - dist_PO6 * (R6 * z5_in_6);

        % Branch 3: q1/q2 (shoulder, Case C1 / C2)
        q1_2_options = [
            atan2(vec_O2P(2), vec_O2P(1)), acos(vec_O2P(3) / norm(vec_O2P)); % C1
            atan2(-vec_O2P(2), -vec_O2P(1)), -acos(vec_O2P(3) / norm(vec_O2P)) % C2
        ];

        for s12 = 1:2
            cq1 = q1_2_options(s12, 1);
            cq2 = q1_2_options(s12, 2);

            % q3 (upper arm)
            R2_t = rot_z(cq1) * rot_x(-pi/2) * rot_z(cq2);
            z3_t = vec_O2P / norm(vec_O2P);
            y3_t = cross(vec_O2P, (p6 - p2)) / norm(cross(vec_O2P, (p6 - p2)));
            x3_t = cross(y3_t, z3_t);
            x3_in_2 = R2_t' * x3_t;
            cq3 = atan2(x3_in_2(3), x3_in_2(1));

            % q5 (forearm)
            p4_t = p2 + d3*z3_t + a4*x3_t;
            z5_t = R6 * [sin(cq6); cos(cq6); 0];
            vec_HO4 = p4_t - (p6 - d5*z5_t);
            R5_6_t = rot_x(pi/2) * rot_z(cq6);
            vec_5_HO4 = R5_6_t * R6' * vec_HO4;
            cq5 = -atan2(vec_5_HO4(2), vec_5_HO4(1));

            SOLS_RAD = [SOLS_RAD; cq1, cq2, cq3, cq4, cq5, cq6, q7];
        end
    end

%% 5. OUTPUT & FILTERING
% Franka Emika joint limits (in degrees)
limits = [
    -166, 166;   % q1
    -101, 101;   % q2
    -166, 166;   % q3
    -176, -4;    % q4 (Franka spec: upper limit -4 deg, not +4 deg)
    -166, 166;   % q5
      -1, 215;   % q6
    -166, 166    % q7
];

SOLS_DEG = rad2deg(SOLS_RAD);
valid = all(SOLS_DEG >= limits(:,1)' & SOLS_DEG <= limits(:,2)', 2);

T = array2table(SOLS_DEG, 'VariableNames', {'q1','q2','q3','q4','q5','q6','q7'});
T.Valid = valid;

fprintf('\n--- All 4 mathematical solutions A2xB1/B2xC1/C2 (in degrees) ---\n');
disp(T);

if any(valid)
    fprintf('Best solution:\n');
    valid_idx = find(valid);
    [~, best_sub] = min(sum((SOLS_RAD(valid_idx,:) - qa').^2, 2));
    disp(T(valid_idx(best_sub), :));
else
    warning('No solution is within joint limits!');
end
