global J m T grv
%system parameters
grv = [0;0;9.81];
m = 4; % mass of the drone
J = diag([0.006666666667, 0.03333333333, 0.03333333333]); % inertia matrix of the drone

% thurster positions in body frame
RR = [10	0	5
    -10	0	5
    10	5	0
    -10	5	0
    10	0	-5
    -10	0	-5
    10	-5	0
    -10	-5	0
    15	0	0
    -15	0	0]'*1e-2; %each row is the position of a propeller

%thurster directions in the body frame
U = [0	0	1
    0	0	1
    0	1	0
    0	1	0
    0	0	-1
    0	0	-1
    0	-1	0
    0	-1	0
    1	0	0
    -1	0	0]'; % each row is the direction of the thrust
direction = [1 1 1 1 1 1 1 1 1 -1]; % directio of the moments
kf = 1; %rho*D^4*CT, thrust coeff
km = 1; %rho*D^5*CM, torque coeff
T = []; % thurst to force and torque transformation

% compute the force and torque produced by each thurster 
for I = 1:width(RR)
    aux = [kf*U(:,I); 
           kf*skew(RR(:,I))*U(:,I)-km*U(:,I)*direction(I)];
    T = [T aux];
end 

%simulation parameters
TSPAN = [0 100];
JSPAN = [0 10];
rule = 1; %priority for jumps
p0 = [0; 0; 0]; % initial position
v0 = [0; 0; 0]; % initial velocity
R0 = angle2dcm(pi/4,0,0); % initial orientation
w0 = [0; 0; 0]; % initial angular velocity
t0 = 0;
xi0 = [p0; v0; reshape(R0, 9, 1); w0;t0]; % initial state vector

[t,j,xi] = HyEQsolver(@F,@G,@C,@D,xi0, TSPAN, JSPAN,rule);

% Extract unit-quaternions from R and plot them
quaternions = zeros(length(t), 4);
for i = 1:length(t)
    R = reshape(xi(i, 7:15), 3, 3);
    q = rotm2quat(R); % Convert rotation matrix to quaternion
    quaternions(i, :) = q;
end

figure;
subplot(2,1,1)
plot(t, quaternions);
xlabel('Time (s)');
ylabel('Quaternion Components');
legend('q_0', 'q_1', 'q_2', 'q_3');
title('Unit-Quaternions Over Time');
subplot(2,1,2)
rpm = zeros(length(t), width(RR)); % Initialize RPMs matrix
for I = 1:numel(t)
    [~,aux] = F(xi(I,:)');
    rpm(I,:) = aux(1:10); % Extract RPMs for each time step
end
plot(t, rpm);
xlabel('Time (s)');
ylabel('Motor RPMs');

%% 3D Animation of 2U Satellite
% Define 2U satellite dimensions (in meters)
length_2U = 0.2; % 2U CubeSat: 20cm long
width_2U = 0.1;  % 10cm wide
height_2U = 0.1; % 10cm high

% Define vertices of the box centered at origin
verts = 0.5 * [
    -length_2U, -width_2U, -height_2U;
    length_2U, -width_2U, -height_2U;
    length_2U,  width_2U, -height_2U;
    -length_2U,  width_2U, -height_2U;
    -length_2U, -width_2U,  height_2U;
    length_2U, -width_2U,  height_2U;
    length_2U,  width_2U,  height_2U;
    -length_2U,  width_2U,  height_2U
];

faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
figure;
axis equal;
grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('2U Satellite 3D Animation');
view(3);
hold on;

% Set axis limits based on trajectory
margin = 0.1;
traj = xi(:,1:3);
min_traj = min(traj) - margin;
max_traj = max(traj) + margin;
axis([min_traj(1) max_traj(1) min_traj(2) max_traj(2) min_traj(3) max_traj(3)]);

h_patch = patch('Vertices', verts, 'Faces', faces, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.7);

for i = 1:10:length(t)
    % Get position and orientation
    pos = xi(i,1:3)';
    R = reshape(xi(i,7:15),3,3);
    % Transform vertices
    verts_trans = (R * verts')' + pos';
    set(h_patch, 'Vertices', verts_trans);
    % Optionally, plot trajectory so far
    plot3(traj(1:i,1), traj(1:i,2), traj(1:i,3), 'k-', 'LineWidth', 1);
    drawnow;
    if i == 1
        pause(1); % Pause at start
    else
        pause(0.01);
    end
end
hold off;

function [pd,vd,ad,Rd,wd,dwd] = trajectory(t)
    % Define a simple trajectory for the drone
    pd = [0; 0; 0]; % position
    vd = [0; 0; 0]; % velocity
    ad = [0; 0; 0]; % acceleration
    wd = [1; 0; 0]; % angular velocity
    Rd = expm(skew(wd)*t); % orientation (identity matrix)
    dwd = [0; 0; 0]; % angular acceleration
end

function [dxi,rpm] = F(xi)
    global J m T grv
    p = xi(1:3);
    v = xi(4:6);
    R = reshape(xi(7:15), 3, 3);
    w = xi(16:18);
    t = xi(19); % time
    % tracking errors
    [pd,vd,ad,Rd,wd,dwd] = trajectory(t);
    ep = p - pd; % position error
    ev = v - vd; % velocity error
    eR = Rd' * R; % orientation error (R * Rd' gives the relative rotation)
    ew = w - eR'*wd; % angular velocity error
    %control commands
    Fctrl = ad-m*grv-ev -ep; % simple PD control
    Mctrl = J*eR'*dwd+skew(eR'*wd)*J*eR'*wd-inverse_skew(eR-eR')-ew; 
    % individual motors rpm >> TODO: replace rpm with thurst
    rpm = T \ [Fctrl;Mctrl]; % solve for motor RPMs
    % rpm to forces and moments
    aux = T * rpm; % forces and moments from RPMs
    F = aux(1:3); % forces
    M = aux(4:6); % moments
    %dynamics
    dp = v;
    dv = grv+F/m;
    dR = R * skew(w);
    dw = J \ (M - skew(w) * J * w);
    dxi = [dp; dv; reshape(dR, 9, 1); dw; 1];
end

function out = G(xi)
    % Jump map
    out = xi; 
end

function out = C(xi)
    % Flow set
    out = 1; %always in the flow set
end

function out = D(xi)
    % Jump set
    out = 0; %never in the jump set
end

function out = inverse_skew(S)
    % Computes the inverse of a skew-symmetric matrix
    out = [S(3,2); S(1,3); S(2,1)];
end

function S = skew(v)
    % Computes the skew-symmetric matrix of a vector
    S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

