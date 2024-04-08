%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Upstairs NLP Formulation   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

addpath(genpath('./poly5'));

%% parameters
dt = 0.004;  % control period: 4 ms
T_fly = 0.8; % swing time

% boundary constraints —— given by hands —— foot placements determine
param.x0 = 0.00;      param.dx0 = 0;      param.ddx0 = 0;
param.x3 = 0.30;      param.dx3 = 0;      param.ddx3 = 0;

param.z0 = 0.08;      param.dz0 = 0;      param.ddz0 = 0;
param.z3 = 0.30;      param.dz3 = 0;      param.ddz3 = 0;

% total time and time interval 
param.T = T_fly;
param.dt = T_fly/12;  % 12 just for simplify, of course, the value is also flexible.

% 变量顺序： x1   dx1  ddx1       x2      dx2  ddx2         z1    dz1  ddz1        z2    dz2   ddz2
x0 =     [-0.05,  0,    0,  param.x3,   -0.2,    0,   param.z3,    0,    0, param.z3,  -0.2,     0]';

%% problem setup
problem.guess.x = x0;
problem.options.nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-6,...
    'FinDiffType', 'forward', ...
    'Algorithm', 'interior-point', ... % interior-point
    'MaxFunEvals',1e3,...
    'MaxIter', 100);   %options for fmincon

problem.options.nlpOpt.GradConstr = 'on';
problem.options.nlpOpt.GradObj = 'on';
problem.options.nlpOpt.DerivativeCheck = 'off'; % 'on'  optimset: 'FinDiffType', 'forward'

% KEY LINE: Solve the problem
soln = MX_footTra_Optimizer(problem, param);

%% plot the results at each control period
[foot_x, foot_z] = getFootTraFromSoln_grad(soln.x, dt, param);

figure(11); clf;
plot(foot_x, foot_z, 'LineWidth', 3)
hold on
plot([soln.x(1), soln.x(4)], [soln.x(7), soln.x(10)], 'b*', 'MarkerSize', 16, 'LineWidth', 3);
plot([param.x0, param.x0, param.x3], [param.z0, param.z3, param.z3], 'r-', 'LineWidth', 5)

[foot_x_node, foot_z_node] = getFootTraFromSoln_grad(soln.x, param.dt, param);

%% animation
figure(12); clf;
for i = 1:length(foot_x)
    hold off
    plot(foot_x, foot_z, 'LineWidth', 3)
    hold on    
    plot([param.x0, param.x0, param.x3], [param.z0, param.z3, param.z3], 'k-', 'LineWidth', 5)
    plot(foot_x(i), foot_z(i), 'r.', 'MarkerSize', 30)
    drawnow;
    pause(dt);
end

%% compare with the results from cpp-ipopt
copyfile ..\test_ipopt_cmake\build\xSoln_cpp.dat
x_ipopt = load('xSoln_cpp.dat');
[foot_x1, foot_z1] = getFootTraFromSoln_grad(x_ipopt, dt, param);

figure(111); clf;
p1 = plot(foot_x, foot_z, 'LineWidth', 3);
hold on
p2 = plot(foot_x1, foot_z1, 'LineWidth', 3);
plot([soln.x(1), soln.x(4)], [soln.x(7), soln.x(10)], 'b*', 'MarkerSize', 16, 'LineWidth', 3);
plot([x_ipopt(1), x_ipopt(4)], [x_ipopt(7), x_ipopt(10)], 'r*', 'MarkerSize', 16, 'LineWidth', 3);
p3 = plot([param.x0, param.x0, param.x3], [param.z0, param.z3, param.z3], 'k-', 'LineWidth', 5);
legend([p1, p2, p3], {'fmincon', 'ipopt', 'stair'})

[foot_x1_node, foot_z1_node] = getFootTraFromSoln_grad(x_ipopt, param.dt, param);

%% remove the poly5 from path
rmpath(genpath('./poly5'));

%% sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    Sub-functions     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function soln = MX_footTra_Optimizer(problem, param)

G   = problem.guess;
Opt = problem.options;

P.objective = @(z)( myObjective_foot_grad(z) );

P.nonlcon = @(z)( myConstraint_lite_grad(z, param) );

P.x0 = G.x;
P.lb = [];
P.ub = [];
P.Aineq = []; P.bineq = [];
P.Aeq = []; P.beq = [];
P.options = Opt.nlpOpt;
P.solver = 'fmincon';

% solve
tic; 
[zSoln, ~, exitFlag, output] = fmincon(P);
nlpTime = toc;

soln.x = zSoln;

soln.nlpTime = nlpTime;

soln.info.exitFlag = exitFlag;
soln.info.output = output;

soln.problem = problem;

end

function [cost, costGrad] = myObjective_foot_grad(z)
% I do not set an obejctive function!!!
cost = 0;

costGrad = zeros(1, 12);
end

function [c, ceq, cGrad, ceqGrad] = myConstraint_lite_grad(z, mx)

[foot_x, foot_z, foot_x_grad, foot_z_grad] = getFootTraFromSoln_grad(z, mx.dt, mx);
nt = length(foot_x);

%% foot trajectory constraints
nn = floor(nt/3);
% -0.08 <= x <= 0.0;     0.1 <= z <= 0.25
idx_temp = 2:(nn+1); % 2:5
idx_n = length(idx_temp);
c_x1 = addBoundBox(foot_x(idx_temp)',   -0.05*ones(idx_n, 1),           0.00*ones(idx_n, 1));
c_z1 = addBoundBox(foot_z(idx_temp)',   mx.z0*ones(idx_n, 1), (mx.z3 + 0.05)*ones(idx_n, 1));

c_x1_grad = zeros(2 * idx_n, 12);
c_x1_grad(1:idx_n, 1:6) = - foot_x_grad(idx_temp, :);
c_x1_grad(idx_n + (1:idx_n), 1:6) = foot_x_grad(idx_temp, :);

c_z1_grad = zeros(2 * idx_n, 12);
c_z1_grad(1:idx_n, 7:12) = - foot_z_grad(idx_temp, :);
c_z1_grad(idx_n + (1:idx_n), 7:12) = foot_z_grad(idx_temp, :);

% 0.0 <= x <= 0.24;     0.21 <= z <= 0.24
idx_temp = (nn+2):(nn+nn+1); % 6:9
idx_n = length(idx_temp);
c_x2 = addBoundBox(foot_x(idx_temp)',            0.00*ones(idx_n, 1), (mx.x3 - 0.01)*ones(idx_n, 1));
c_z2 = addBoundBox(foot_z(idx_temp)',  (mx.z3 + 0.01)*ones(idx_n, 1), (mx.z3 + 0.04)*ones(idx_n, 1));

c_x2_grad = zeros(2 * idx_n, 12);
c_x2_grad(1:idx_n, 1:6) = - foot_x_grad(idx_temp, :);
c_x2_grad(idx_n + (1:idx_n), 1:6) = foot_x_grad(idx_temp, :);

c_z2_grad = zeros(2 * idx_n, 12);
c_z2_grad(1:idx_n, 7:12) = - foot_z_grad(idx_temp, :);
c_z2_grad(idx_n + (1:idx_n), 7:12) = foot_z_grad(idx_temp, :);

% 0.24 <= x <= 0.26;     0.2 <= z <= 0.23
idx_temp = (nn+nn+2):(nt-1); % 10:12
idx_n = length(idx_temp);
c_x3 = addBoundBox(foot_x(idx_temp)',  (mx.x3 - 0.01)*ones(idx_n, 1),           mx.x3*ones(idx_n, 1));
c_z3 = addBoundBox(foot_z(idx_temp)',           mx.z3*ones(idx_n, 1),  (mx.z3 + 0.03)*ones(idx_n, 1));

c_x3_grad = zeros(2 * idx_n, 12);
c_x3_grad(1:idx_n, 1:6) = - foot_x_grad(idx_temp, :);
c_x3_grad(idx_n + (1:idx_n), 1:6) = foot_x_grad(idx_temp, :);

c_z3_grad = zeros(2 * idx_n, 12);
c_z3_grad(1:idx_n, 7:12) = - foot_z_grad(idx_temp, :);
c_z3_grad(idx_n + (1:idx_n), 7:12) = foot_z_grad(idx_temp, :);

% 从第6个开始，x方向只增不减
idx_upp = (nn+2):nt;
idx_low = idx_upp - 1;
c_x4 = foot_x(idx_upp) - foot_x(idx_low);

idx_n = length(idx_upp);
c_x4_grad = zeros(idx_n, 12);

c_x4_grad(:, 1:6) = foot_x_grad(idx_upp, :) - foot_x_grad(idx_low, :);

% 从第5个开始，z方向只增不减
idx_upp = (nn+1):nt;
idx_low = idx_upp - 1;
c_z4 = foot_z(idx_upp) - foot_z(idx_low);

idx_n = length(idx_upp);
c_z4_grad = zeros(idx_n, 12);

c_z4_grad(:, 7:12) = foot_z_grad(idx_upp, :) - foot_z_grad(idx_low, :);

c = [c_x1; c_x2; c_x3; c_z1; c_z2; c_z3; -c_x4'; c_z4'];
ceq = [];

cGrad = [c_x1_grad; c_x2_grad; c_x3_grad; c_z1_grad; c_z2_grad; c_z3_grad; -c_x4_grad; c_z4_grad]';
ceqGrad = [];

end

function [foot_x, foot_z, foot_x_grad, foot_z_grad] = getFootTraFromSoln_grad(z, dt, mx)

x1 = z(1);      dx1 = z(2);     ddx1 = z(3);
x2 = z(4);      dx2 = z(5);     ddx2 = z(6);
z1 = z(7);      dz1 = z(8);     ddz1 = z(9);
z2 = z(10);     dz2 = z(11);    ddz2 = z(12);

T = mx.T;

% x 
S_x1 = getFivePolyCoeff(T/3, mx.x0, mx.dx0,  mx.ddx0,    x1,    dx1,    ddx1);
S_x2 = getFivePolyCoeff(T/3,    x1,    dx1,     ddx1,    x2,    dx2,    ddx2);
S_x3 = getFivePolyCoeff(T/3,    x2,    dx2,     ddx2, mx.x3, mx.dx3, mx.ddx3);

% z 
S_z1 = getFivePolyCoeff(T/3, mx.z0, mx.dz0,  mx.ddz0,    z1,    dz1,    ddz1);
S_z2 = getFivePolyCoeff(T/3,    z1,    dz1,     ddz1,    z2,    dz2,    ddz2);
S_z3 = getFivePolyCoeff(T/3,    z2,    dz2,     ddz2, mx.z3, mx.dz3, mx.ddz3);

tspan = 0:dt:T;
nt = length(tspan);

foot_x = zeros(1, nt);
foot_z = zeros(1, nt);

for i = 1:nt
    tt = tspan(i);

    if tt <= T/3
        foot_x(i) = getFiveOrderPoly(S_x1, tt);
        foot_z(i) = getFiveOrderPoly(S_z1, tt);
    elseif tt > T/3 && tt <= 2 * T/3
        foot_x(i) = getFiveOrderPoly(S_x2, tt - T/3);
        foot_z(i) = getFiveOrderPoly(S_z2, tt - T/3);
    else
        foot_x(i) = getFiveOrderPoly(S_x3, tt - 2 * T/3);
        foot_z(i) = getFiveOrderPoly(S_z3, tt - 2 * T/3);
    end
end

if nargout > 2
    foot_x_grad = zeros(nt, 6);
    foot_z_grad = zeros(nt, 6);
    for i = 1:nt
        tt = tspan(i);

        if tt <= T/3
%             foot_x(i) = getFiveOrderPoly(S_x1, tt);
            foot_x_grad(i, 1:3) = autoGen_jacobianRight3(tt, T/3);

%             foot_z(i) = getFiveOrderPoly(S_z1, tt);
            foot_z_grad(i, 1:3) = autoGen_jacobianRight3(tt, T/3);
        elseif tt > T/3 && tt <= 2 * T/3
%             foot_x(i) = getFiveOrderPoly(S_x2, tt - T/3);
            foot_x_grad(i, :) = autoGen_jacobianAll6(tt - T/3, T/3);        

%             foot_z(i) = getFiveOrderPoly(S_z2, tt - T/3);
            foot_z_grad(i, :) = autoGen_jacobianAll6(tt - T/3, T/3);
        else
%             foot_x(i) = getFiveOrderPoly(S_x3, tt - 2 * T/3);
            foot_x_grad(i, 4:6) = autoGen_jacobianLeft3(tt - 2 * T/3, T/3); 

%             foot_z(i) = getFiveOrderPoly(S_z3, tt - 2 * T/3);
            foot_z_grad(i, 4:6) = autoGen_jacobianLeft3(tt - 2 * T/3, T/3);
        end
    end
end

end

function c = addBoundBox(x, xLow, xUpp)

nx = numel(x);
c = zeros(2*nx, 1);
idxLeft = 1:nx;
c(idxLeft) = -x + xLow;
idxRight = nx + (1:nx);
c(idxRight) = x - xUpp;

end