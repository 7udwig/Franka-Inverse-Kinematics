% forward_kinematics_modDH.m  --  Franka Emika Panda Forward Kinematics (Modified DH)
%
% Computes the end-effector (flange) pose from given joint angles using the
% Modified Denavit-Hartenberg (MDH) parameterisation.
% Kinematic chain: T_flange = T1 * T2 * ... * T7 * T_flange_offset
%
% Frame convention : Flange frame  (d_f = 0.107 m offset appended at the end).
%   Results match IK_Franka.m / validate_IK.m in 'flange' mode.
%
% Inputs  (hard-coded in script, modify as needed):
%   q  -- (7x1) joint angles in degrees, converted internally to radians.
%
% Outputs (printed to console):
%   TCP position [mm] and Roll-Pitch-Yaw orientation [deg] (ZYX convention).
%
% Reference:
%   He & Liu, "Analytical Inverse Kinematics for Franka Emika Panda --
%   a Geometrical Solver for the 7-DOF Robot Arm", 2021.

% Joint angles in degrees
%q = [10; 16; 10; -100; 14; 10; 10];
%q = [10.00;   20.00;   10.00; -100.00;   10.00;  100.00;   0];
%q = [118.68;   35.20;    0.00; -103.83;  158.32; -123.24;   10.00];
q = [118.68;   35.20;    0.00; -103.83;  158.32; 123.24;   10.00];

% Convert to radians
q = deg2rad(q);

% Add theta offsets (as in original parameterisation)
A1 = q(1) + 0;
A2 = q(2) + 0;
A3 = q(3) + 0;
A4 = q(4) + 0;
A5 = q(5) + 0;
A6 = q(6) + 0;
A7 = q(7) + 0;

% Modified DH parameters (d, a, alpha, theta)
%              1        2      3       4       5       6     7
d     = [333.0;    0;     316;    0;     384;     0; 0];
a     = [    0;    0;     0;     82.5;      -82.5;        0; 88];
alpha = [ 0;     -pi/2;    pi/2;  pi/2;   -pi/2;     pi/2; pi/2];
theta = [ A1; A2; A3; A4; A5; A6; A7];

% Modified DH transformation matrices
T1 = modifiedDHTransform(a(1), alpha(1), d(1), theta(1));
T2 = modifiedDHTransform(a(2), alpha(2), d(2), theta(2));
T3 = modifiedDHTransform(a(3), alpha(3), d(3), theta(3));
T4 = modifiedDHTransform(a(4), alpha(4), d(4), theta(4));
T5 = modifiedDHTransform(a(5), alpha(5), d(5), theta(5));
T6 = modifiedDHTransform(a(6), alpha(6), d(6), theta(6));
T7 = modifiedDHTransform(a(7), alpha(7), d(7), theta(7));

% Final transform: flange offset (d_f = 0.107 m along z7)
T8 = modifiedDHTransform(0, 0, 107, 0);


% Total TCP pose
TCP_Pose = T1*T2*T3*T4*T5*T6*T7*T8;

% Extract position
position = TCP_Pose(1:3,4);
fprintf('\n--- MDH TCP position ---\n');
fprintf('X: %.2f mm\n', position(1));
fprintf('Y: %.2f mm\n', position(2));
fprintf('Z: %.2f mm\n', position(3));

% Extract rotation matrix (upper-left 3x3)
R = TCP_Pose(1:3, 1:3);

% Compute Roll-Pitch-Yaw angles (ZYX convention)
% pitch (beta)
pitch = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));

% yaw (alpha)
yaw = atan2(R(2,1), R(1,1));

% roll (gamma)
roll = atan2(R(3,2), R(3,3));

% Convert to degrees for readability
angles_deg = rad2deg([roll; pitch; yaw]);

% Print orientation
fprintf('\n--- MDH TCP orientation (Roll-Pitch-Yaw) ---\n');
fprintf('Roll  (about X): %.2f deg\n', angles_deg(1));
fprintf('Pitch (about Y): %.2f deg\n', angles_deg(2));
fprintf('Yaw   (about Z): %.2f deg\n', angles_deg(3));

% --- modifiedDHTransform ---
% Homogeneous transformation matrix for one Modified DH joint.
%   a     -- link length  [mm]
%   alpha -- link twist   [rad]
%   d     -- link offset  [mm]
%   theta -- joint angle  [rad]
%   Returns (4x4) homogeneous transformation matrix T.
function [T] = modifiedDHTransform(a, alpha, d, theta)
    T = [cos(theta),             -sin(theta),             0,              a;
         sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -sin(alpha),  -sin(alpha)*d;
         sin(theta)*sin(alpha),  cos(theta)*sin(alpha),   cos(alpha),   cos(alpha)*d;
         0,                      0,                       0,              1];
end
