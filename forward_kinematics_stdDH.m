% forward_kinematics_stdDH.m  --  Franka Emika Panda Forward Kinematics (Standard DH)
%
% Computes the end-effector (flange) pose from given joint angles using the
% standard Denavit-Hartenberg (DH) parameterisation.
% Kinematic chain: T_flange = T1 * T2 * ... * T7 * T_flange_offset
%
% Frame convention : Flange frame  (d_f = 0.107 m offset appended at the end).
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
%q = [10; 16; 10; -30; 14; 10; 10];
%q = [10.00;   20.00;   10.00; -100.00;   10.00;  100.00;   0];
%q = [118.68;   35.20;    0.00; -103.83;  158.32; -123.24;   10.00];
q = [118.68;   35.20;    0.00; -103.83;  158.32; 123.24;   10.00];

for i=1:7
    q(i,1) = deg2rad(q(i,1));
end

A1 = q(1,1) + 0;
A2 = q(2,1) + 0;
A3 = q(3,1) + 0;
A4 = q(4,1) + 0;
A5 = q(5,1) + 0;
A6 = q(6,1) + 0;
A7 = q(7,1) + 0;


% Standard DH parameters
%                 1           2        3       4       5        6    7
d     = [333.0;    0;     316;    0;     384;     0; 0];
a     = [    0;     0;     82.5;      -82.5;        0; 88; 0];
alpha = [-pi/2;    pi/2;  pi/2;   -pi/2;     pi/2; pi/2;  0;];
theta = [ A1; A2; A3; A4; A5; A6; A7];

% Compute homogeneous transformation matrices for each joint
T1 = standardDHTransform(d(1), A1, a(1), alpha(1));
T2 = standardDHTransform(d(2), A2, a(2), alpha(2));
T3 = standardDHTransform(d(3), A3, a(3), alpha(3));
T4 = standardDHTransform(d(4), A4, a(4), alpha(4));
T5 = standardDHTransform(d(5), A5, a(5), alpha(5));
T6 = standardDHTransform(d(6), A6, a(6), alpha(6));
T7 = standardDHTransform(d(7), A7, a(7), alpha(7));
T8 = standardDHTransform(107, 0, 0, 0); % Flange offset


% TCP transformation
TCP_Pose = T1*T2*T3*T4*T5*T6*T7*T8;

position = TCP_Pose(1:3,4);
fprintf('\n--- DH TCP position ---\n');
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
fprintf('\n--- DH TCP orientation (Roll-Pitch-Yaw) ---\n');
fprintf('Roll  (about X): %.2f deg\n', angles_deg(1));
fprintf('Pitch (about Y): %.2f deg\n', angles_deg(2));
fprintf('Yaw   (about Z): %.2f deg\n', angles_deg(3));

% --- standardDHTransform ---
% Homogeneous transformation matrix for one standard DH joint.
%   d      -- link offset  [mm]
%   thetar -- joint angle  [rad]
%   a      -- link length  [mm]
%   alpha  -- link twist   [rad]
%   Returns (4x4) homogeneous transformation matrix DHT.
function [DHT] = standardDHTransform(d, thetar, a, alpha)

    DHT = [cos(thetar), -sin(thetar)*cos(alpha),  sin(thetar)*sin(alpha), a*cos(thetar)       ;
               sin(thetar),  cos(thetar)*cos(alpha), -cos(thetar)*sin(alpha), a*sin(thetar);
               0,              sin(alpha),                cos(alpha),                     d;
               0,              0,                            0,                          1];
end
