%% IK_Franka.m  --  Franka Emika Panda Analytical Inverse Kinematics (single solution)
%
% Computes one IK solution for the Franka Research 3 / Panda (7-DOF) using the
% closed-form geometric method of He & Liu (2021), A2-branch only.
%
% Frame convention (selectable via frame_mode):
%   'flange' (default): flange pose used directly. Corresponds to forward_kinematics_modDH.m.
%   'EE':               EE pose pre-transformed to flange pose via T_EE * inv(EE_offset),
%                       then the flange kernel is called. EE_offset = F_T_EE (official
%                       Franka default: Rz(+pi/4), 0.1034 m along z).
%
% Inputs  (set directly in script):
%   p_ee       -- (3x1) TCP position [m]
%   roll/pitch/yaw -- ZYX Euler angles of TCP orientation [rad]
%   qa         -- (7x1) reference configuration [rad]  (fixes q7; selects B/C branch)
%   frame_mode -- 'flange' or 'EE'
%
% Outputs (printed to console):
%   Computed joint angles q1..q7 [deg] alongside reference values.
%
% Reference:
%   He & Liu, "Analytical Inverse Kinematics for Franka Emika Panda --
%   a Geometrical Solver for the 7-DOF Robot Arm", 2021.

%% 0. EXPLANATION
% q7: redundancy parameter (fixed to qa(7))
% q6: Case B1/B2 based on dot(O2O6, x5(qa)) -- paper Eq. (83)
% q1 and q2: Case C1/C2 based on sign of qa(2) -- paper Eq. (84)
% A2-branch only: q4 = phi1 + phi2 + phi3 - 2*pi  (paper Eq. (48))

%% 1. INPUT VALUES
%p_ee = [0.56689; 0.21491; 0.35557];
%roll = deg2rad(-167.15); pitch = deg2rad(19.93); yaw = deg2rad(18.80);
%qa = deg2rad([10; 20; 10; -100; 10; 100; 0]);

p_ee = [-0.33938; 0.58864; 0.23309]; % TCP position [m]
roll = deg2rad(-36.52); pitch = deg2rad(83.16); yaw = deg2rad(100.05);
qa = deg2rad([118.68; 35.20; 0.00; -103.83; 158.32; 123.24; 10.00]); % Reference (q7 + branch selection)

frame_mode = 'flange'; % 'flange' (default) or 'EE'

%% 2. CONSTANTS & DH GEOMETRY
d1 = 0.333; d3 = 0.316; a4 = 0.0825; d5 = 0.384; a5 = -0.0825; a7 = 0.088; df = 0.107;

% F_T_EE: flange-to-EE transform (official Franka default: theta_EE=+pi/4, 0.1034 m).
% Replace with the robot's actual F_T_EE if it differs from the factory default.
theta_ee = pi/4;  dee_ee = 0.1034;
EE_offset = [cos(theta_ee), -sin(theta_ee), 0, 0;
             sin(theta_ee),  cos(theta_ee), 0, 0;
             0,              0,             1, dee_ee;
             0,              0,             0, 1];

L24 = sqrt(d3^2 + a4^2); % Length O2-O4 (upper arm)
L46 = sqrt(d5^2 + a5^2); % Length O4-O6 (forearm)

rot_x = @(phi) [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
rot_z = @(phi) [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

%% 3. BASE CALCULATIONS (R_ee, p6, q4)
R_ee = rot_z(yaw) * rot_y(pitch) * rot_x(roll); % ZYX Euler -> rotation matrix
if strcmp(frame_mode, 'EE')
    % Pre-transform: convert EE pose to flange pose, then run the flange kernel.
    T_EE = [R_ee, p_ee; 0, 0, 0, 1];
    T_flange = T_EE / EE_offset; % T_flange = T_EE * inv(EE_offset)
    R_ee = T_flange(1:3, 1:3);
    p_ee = T_flange(1:3, 4);
end
z_ee = R_ee(:, 3);
p7 = p_ee - df * z_ee; % Frame-7 centre (flange frame, dee=0)

x6 = R_ee * [cos(qa(7)); -sin(qa(7)); 0]; % x6-axis (flange frame, q7_eff=qa(7))
p6 = p7 - a7 * x6; % Wrist position O6
p2 = [0; 0; d1];   % Shoulder position O2
O2O6 = p6 - p2;
dist_O2O6 = norm(O2O6);

% q4 -- A2-branch: phi1+phi2+phi3-2*pi  (paper Eq. (48))
phi1 = atan2(d3, a4);
phi2 = atan2(d5, abs(a5));
phi3 = acos(clamp((L24^2 + L46^2 - dist_O2O6^2) / (2 * L24 * L46), -1, 1));
q4 = phi1 + phi2 + phi3 - 2*pi;

%% 4. q6 -- Case B1/B2 (paper Eq. (54, 83))
y6 = -z_ee;
z6 = cross(x6, y6);
R6 = [x6, y6, z6];

phi_O2O6O4 = acos(clamp((dist_O2O6^2 + L46^2 - L24^2) / (2 * dist_O2O6 * L46), -1, 1));
phi_HO6O4 = atan(abs(a5) / d5);
angle_O2O6H = phi_O2O6O4 + phi_HO6O4;

vec_6_O2O6 = R6' * O2O6;
phi6 = atan2(vec_6_O2O6(2), vec_6_O2O6(1));
psi6 = asin(clamp((dist_O2O6 * cos(angle_O2O6H)) / norm(vec_6_O2O6(1:2)), -1, 1));

% B1/B2: compare current x5-axis with O2O6 -- paper Eq. (83)
x5_curr = get_x5_world(qa);
if dot(O2O6, x5_curr) <= 0
    q6 = pi - psi6 - phi6; % Case B1
else
    q6 = psi6 - phi6;       % Case B2
end
q6 = atan2(sin(q6), cos(q6)); % normalise to (-pi, pi]

%% 5. q1, q2 -- Case C1/C2 (paper Eq. (66, 84))
phi_O3O2O4 = atan(a4 / d3);
phi_O4O2O6 = acos(clamp((L24^2 + dist_O2O6^2 - L46^2) / (2 * L24 * dist_O2O6), -1, 1));
angle_PO2O6 = phi_O3O2O4 + phi_O4O2O6;
angle_O2PO6 = phi_O2O6O4 + phi3 + phi1 - angle_O2O6H - pi/2;

dist_PO6 = dist_O2O6 * sin(angle_PO2O6) / sin(angle_O2PO6);
z5_in_6 = [sin(q6); cos(q6); 0];
vec_O2P = O2O6 - dist_PO6 * (R6 * z5_in_6); % Vector shoulder O2 -> auxiliary point P

% C1/C2: sign of qa(2) determines shoulder configuration -- paper Eq. (84)
if qa(2) >= 0
    q1 = atan2(vec_O2P(2), vec_O2P(1));
    q2 = acos(clamp(vec_O2P(3) / norm(vec_O2P), -1, 1));
else
    q1 = atan2(-vec_O2P(2), -vec_O2P(1));
    q2 = -acos(clamp(vec_O2P(3) / norm(vec_O2P), -1, 1));
end

%% 6. q3 (upper-arm rotation) -- Modified DH frame 2
R2 = rot_z(q1) * rot_x(-pi/2) * rot_z(q2);

z3 = vec_O2P / norm(vec_O2P);
y3 = cross(vec_O2P, O2O6) / norm(cross(vec_O2P, O2O6));
x3 = cross(y3, z3);

x3_in_2 = R2' * x3;
q3 = atan2(x3_in_2(3), x3_in_2(1)); % angle between x2 and x3 about z3

%% 7. q5 (forearm rotation) -- paper Eq. (36, 67)
p4 = p2 + d3*z3 + a4*x3;
z5 = R6 * [sin(q6); cos(q6); 0];
vec_HO4 = p4 - (p6 - d5*z5);

R5_6 = rot_x(pi/2) * rot_z(q6); % frame 5 -> 6 transition (paper Eq. (36))
vec_5_HO4 = R5_6 * R6' * vec_HO4; % note: R5_6 is not transposed
q5 = -atan2(vec_5_HO4(2), vec_5_HO4(1));

%% OUTPUT
fprintf('--- IK result [%s-frame] (in degrees) ---\n', frame_mode);
res = rad2deg([q1, q2, q3, q4, q5, q6, qa(7)]);
labels = {'q1','q2','q3','q4','q5','q6','q7'};
for i = 1:7
    fprintf('%s: %8.2f (target: %8.2f)\n', labels{i}, res(i), rad2deg(qa(i)));
end

%% --- HELPER FUNCTIONS ---

function R = rot_y(phi)
    R = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
end

function out = clamp(in, low, high)
    out = max(min(in, high), low);
end

% x5-axis in world coordinates via Modified DH -- branch selection B (paper Eq. (83))
function x5 = get_x5_world(q)
    d1 = 0.333; d3 = 0.316; a4 = 0.0825; a5 = -0.0825; d5 = 0.384;
    MDH = @(al, a, d, th) [cos(th), -sin(th), 0, a;
                            sin(th)*cos(al), cos(th)*cos(al), -sin(al), -d*sin(al);
                            sin(th)*sin(al), cos(th)*sin(al),  cos(al),  d*cos(al);
                            0, 0, 0, 1];
    T05 = MDH(0,0,d1,q(1)) * MDH(-pi/2,0,0,q(2)) * MDH(pi/2,0,d3,q(3)) * ...
          MDH(pi/2,a4,0,q(4)) * MDH(-pi/2,a5,d5,q(5));
    x5 = T05(1:3, 1);
end
