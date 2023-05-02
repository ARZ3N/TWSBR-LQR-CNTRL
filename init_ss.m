%% TWSBR: INITIAL STATE SPACE MODEL
% TEAM 1-
% ANIRUDDHA BHATTACHARJEE
% POORNESH RAVAL
% VAIBHAV RAJSEWATI

%%% System Description Notes
%----------------------------------------------------
%  State vector: [x, v, theta, omega, del, del-dot]
%  Original control var: [Cl, Cr]'
%  Decoupled control vars: [Ctheta, Cdel]'
%  Decoupled system 1 states: [x, v, theta, omega]
%  Decoupled system 2 states: [del, del-dot]
%  Decoupled system 1 controls: [ctheta]
%  Decoupled system 2 controls: [Cdel]
% del: realted to yawing motion
% theta: realted to tilting motion

%% KINEMATIC VARIABLES
M = 1*10; % Mass of the Chassis
m = 0.012*10; % Mass of each wheel
Jz = 0.0144*100; % Moment of Inertia of Chassis about Z-axis
Jy = 0.00838*10; % Moment of Inertia of body about Y-axis
Jwh = 0.0000054*10; % Moment of Inertia of each wheel
L = 0.12; % Distance b/w CG and wheel-axis
R = 0.06; % Wheel radius
D = 0.12; % Distance between two wheels
g = 9.8; % Acceleration due to Gravity

%%% State Space Model
a23 = -1.0 * ( ((M^2) * (L^2) * g) / (M*Jz + 2*(Jz + M*(L^2) * (m + (Jwh/(R^2))) ) ))

a43 = ( ((M^2)*g*L) + 2*M*g*L*(m + (Jwh/(R^2))) ) / ( M*Jz + 2*(Jz + M*(L^2))*(m + (Jwh/(R^2))) )

b21 = ( ((Jz + M*(L^2))/R) + (M*L) ) / ( M*Jz + 2*(Jz + M*(L^2))*(m + (Jwh/(R^2))) )

b41 = -1.0 * ( (((R + L)*M)/R) - 2*(m + (Jwh/(R^2))) ) / ( M*Jz + 2*(Jz + M*(L^2))*(m + (Jwh/(R^2))) )

b61 = (D / (2*R)) / ( Jy + ((D^2)/(2*R))*(m*R + (Jwh/R)) )

%%% Total State Space (SS)
A = [0 1 0 0 0 0;
     0 0 a23 0 0 0;
     0 0 0 1 0 0;
     0 0 a43 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0];

B = [0 0;
     b21 b21;
     0 0;
     b41 b41;
     0 0;
     b61 -1*b61];

C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];

D = [0];

% Input Decoupling Matrix
Udcpl = [0.5 0.5;
         0.5 -0.5];

sys = ss(A, B, C, D);

%% State Space: Decoupled Sys 1: Pitch Dynamics
Asys1 = [0 1 0 0;
         0 0 a23 0;
         0 0 0 1;
         0 0 a43 0;];

Bsys1 = [0;
         b21;
         0;
         b41];

Bsys10 = [0;
         b21;
         0;
         b41];

Csys1 = [1 0 0 0]
         %0 0 1 0];    % Selecting [x, theta]

Dsys1 = [0];

sys1 = ss(Asys1, Bsys10, Csys1, Dsys1);

%% State Space: Decoupled Sys 2: Yaw Dynamics
Asys2 = [0 1;
         0 0];    % [del, del-dot]

Bsys2 = [0;
         b61];

Bsys20 = [0;
          0];

Csys2 = [1 0];    % Selecting [del]

Dsys2 = [0];

sys2 = ss(Asys2, Bsys2, Csys2, Dsys2);

%% LQR Gain for Sys-1
% Circular bounding
% q1: x, q2: v, q3: theta, q4: omega
radius = 10.0;
center = 1.0;
q1 = 1.0;
q2 = 1.0;
Q = {};
for k = 1:10
    q3 = circular(q1, center, radius);
    q4 = circular(q2, center, radius);
    fprintf("q1: %d, q2: %d, q3: %d, q4: %d \n", q1, q2, q3, q4);
    Q{k} = [q1 0 0 0; 0 q2 0 0; 0 0 q3 0; 0 0 0 q4];
    q1 = q1 + 1.0;
    q2 = q2 + 1.0;
end
R = [center];

fprintf("LQR gain matrix: \n");
K = {};
fd_sys = {};
fd_sys_fs = {};
Kp = [2];

for k = 1:10
    fprintf("pass k: %d \n", k);
    K{k} = lqr(Asys1, Bsys1, Q{k}, R);
    %K{k} = lqr(Asys1 - Bsys1*K0 - Bsys1*Kp*Csys1, Bsys1*Kp, Q{k}, R);
    fd_sys_fs{k} = ss(Asys1 - Bsys1*K{k}, Bsys1, Csys1, Dsys1); %#ok<SAGROW> 
    fd_sys{k} = ss(Asys1 - Bsys1*K{k} + Bsys1*Kp*Csys1, -Bsys1*Kp, Csys1, Dsys1);

end
sys1_dcg = dcgain(fd_sys{10});

%% LQR Gain for Sys-2
rc2 = 1.0;
qc1 = 1.0;
Q2 = {};
K2 = {};
for k = 1:10
    qc2 = circular(qc1, rc2, 10.0);
    Q2{k} = [qc1 0; 0 qc2];
    qc1 = qc1 + 1.0;
end
R2 = [rc2];

fd_sys2 = {};
fd_sys2_fs = {};
Kp2 = [10.0];
K0 = lqr(Asys2, Bsys2, eye(2), eye(1));
for k = 1:10
    fprintf("pass k: %d \n", k);
    %K2{k} = lqr(Asys2, Bsys2, Q2{k}, R);
    K2{k} = lqr(Asys2 - Bsys2*K0 - Bsys2*Kp2*Csys2, Bsys2*Kp2, Q2{k}, R2);
    fd_sys2_fs{k} = ss(Asys2 - Bsys2*K2{k}, Bsys2, Csys2, Dsys2);
    fd_sys2{k} = ss(Asys2 - Bsys2*K2{k} - Bsys2*Kp2*Csys2, Bsys2*Kp2, Csys2, Dsys2);
end


%% Simulation Runs and Graphs

% System 1 plots
t = 0:0.1:10;
[i,j] = size(t);
u = [1];
sig = 0;
for k=2:j
    if mod(k, 3) == 0
        %fprintf("hit mod(k,3), k = %d \n", k);
        sig = sig + 0.5;
    else 
        sig = sig - 0.5;
    end
    u(k) =  sig;
end
u = u.*-1;
plot(t, u, 'LineWidth', 1.4);
grid on
title('Reference Signal: r (x desired)')
xlabel("Time (sec)")
ylabel("Magnitude")
%pause;

x0 = [0.0, 0.0, 0.1, 0.0];

u = t.^3;


%{
input_sig = zeros(size(u));
len = size(t);
full_state = zeros(len(2));
for it = 1: len(2)
    full_state(it) = K{10} * (x(it, :).');
end
for it = 1: len(2)
    input_sig(it) = -Kp*(u(it) - y(it)) - full_state;
end
%}
%{ 
% Uncomment to plot Response for Full-state and Error Term 
--------------- FS + ER ----------------------------------------
[y, t, x] = lsim(fd_sys{10}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fd-10 x', 'LineWidth', 1.4);
grid on
hold on
plot(t, x(:, 3), 'DisplayName', 'fd-10 theta', 'LineWidth', 1.4);
hold on

[y, t, x] = lsim(fd_sys{5}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fd-5 x', 'LineWidth', 1.4);
hold on
plot(t, x(:, 3), 'DisplayName', 'fd-5 theta', 'LineWidth', 1.4);

[y, t, x] = lsim(fd_sys{1}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fd-1 x', 'LineWidth', 1.4);
hold on;
plot(t, x(:, 3), 'DisplayName', 'fd-1 theta', 'LineWidth', 1.4);
hold off
legend();
title('Sys-1 States(x, theta) with FS + ER FDBK, Q-R varying ')
xlabel("Time (sec)")
ylabel("States")
%pause;
%}
%{
% Uncomment to plot responses for Sys-1 only Full-State
------------------ FS -------------------------------------------
[y_fs, t_fs, x_fs] = lsim(fd_sys_fs{10}, u, t, x0);
plot(t_fs, x_fs(:, 1), 'DisplayName', 'fd-10 x', 'LineWidth', 1.4);
grid on
hold on
plot(t_fs, x_fs(:, 3), 'DisplayName', 'fd-10 theta', 'LineWidth', 1.4);
hold on

[y_fs, t_fs, x_fs] = lsim(fd_sys_fs{5}, u, t, x0);
plot(t_fs, x_fs(:, 1), 'DisplayName', 'fd-5 x', 'LineWidth', 1.4);
hold on
plot(t_fs, x_fs(:, 3), 'DisplayName', 'fd-5 theta', 'LineWidth', 1.4);
hold on

[y_fs, t_fs, x_fs] = lsim(fd_sys_fs{1}, u, t, x0);
plot(t_fs, x_fs(:, 1), 'DisplayName', 'fd-1 x', 'LineWidth', 1.4);
hold on;
plot(t_fs, x_fs(:, 3), 'DisplayName', 'fd-1 theta', 'LineWidth', 1.4);
hold off
hold off
legend();
title('Sys-1 States(x, theta) for SS with FS FDBK, Q-R varying ')
xlabel("Time (sec)")
ylabel("States")
%pause; 

%-------------- Combined ---------------------------------------------
[y, t, x] = lsim(fd_sys{10}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fd-10 x', 'LineWidth', 1.4);
grid on
hold on
plot(t, x(:, 3), 'DisplayName', 'fd-10 theta', 'LineWidth', 1.4);
hold on
[y_fs, t_fs, x_fs] = lsim(fd_sys_fs{10}, u, t, x0);
plot(t, x_fs(:, 1), 'DisplayName', 'fd-fs-10 x', 'LineWidth', 1.4);
hold on
plot(t_fs, x_fs(:, 3), 'DisplayName', 'fd-fs-10 theta', 'LineWidth', 1.4);
hold off
legend();
title('Comparison of Tracking Control for FS, FS + ER')
xlabel("Time (sec)")
ylabel("States")
%}


%--------System 2 --------------------------------
u = ones(size(t));
x0 = [0.1, 0.0];
[y, t, x] = lsim(fd_sys2{10}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fdsys2-10 del', 'LineWidth', 1.4);
hold on
grid on
plot(t, x(:, 2), 'DisplayName', 'fdsys2-10 del-dot', 'LineWidth', 1.4);
%-------------------------------------------------
[y, t, x] = lsim(fd_sys2{5}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fdsys2-5 del', 'LineWidth', 1.4);
hold on
grid on
plot(t, x(:, 2), 'DisplayName', 'fdsys2-5 del-dot', 'LineWidth', 1.4);
%-------------------------------------------------
[y, t, x] = lsim(fd_sys2{1}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fdsys2-1 del', 'LineWidth', 1.4);
hold on
grid on
plot(t, x(:, 2), 'DisplayName', 'fdsys2-1 del-dot', 'LineWidth', 1.4);
hold off
legend();
title('System 2 states with FS FDBK, Q, R varying')
xlabel("Time (sec)")
ylabel("States")
%--------------------------------------------------

%{
%----------------- sys-2 combined -------------------------
u = ones(size(t));
x0 = [0.1, 0.0];
[y, t, x] = lsim(fd_sys2{10}, u, t, x0);
plot(t, x(:, 1), 'DisplayName', 'fdsys2-10 del', 'LineWidth', 1.4);
hold on
grid on
plot(t, x(:, 2), 'DisplayName', 'fdsys2-10 del-dot', 'LineWidth', 1.4);
[y_fs, t_fs, x_fs] = lsim(fd_sys2_fs{10}, u, t, x0);
plot(t_fs, x_fs(:, 1), 'DisplayName', 'fdsys2-fs-10 del', 'LineWidth', 1.4);
hold on
plot(t_fs, x_fs(:, 2), 'DisplayName', 'fdsys2-fs-10 del-dot', 'LineWidth', 1.4);
hold off
legend();
title("Sys-2 Comparison of FS & FS + ER FDBK")
xlabel('Time (sec)')
ylabel('States')
%}

%% Bounding Functions
function [ret] = circular(x, cen, rad)
    ret = cen + sqrt((rad^2) - ((x - cen)^2));
end

function [ret] = elliptical(x, cen, a, b)
    ret = (b * sqrt(1 - (((x-cen)/a)^2) ) )+ cen;
end
function [ret] = parabolic(x, a, b)
    ret = 4 * a * (x - b)^2;
end

function plot_bode()
    for it = 1:10
    bode(fd_sys{it})
    grid on
    hold on
    end
    legend()
    hold on
    grid on
    OL_sys = ss(Asys1, Bsys1, Csys1, Dsys1);
    bode(OL_sys)
    hold off
end