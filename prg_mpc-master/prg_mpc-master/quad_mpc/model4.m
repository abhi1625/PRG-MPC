clear;

BEGIN_ACADO

acadoSet('problemname','quad_mpc_4');
% System variables
DifferentialState     p_x p_y p_z;
DifferentialState     q_w q_x q_y q_z;
DifferentialState     v_x v_y v_z;
Control               T w_x w_y w_z;
f = acado.DifferentialEquation();

input1 = acado.MexInput;
input2 = acado.MexInputVector;
input3 = acado.MexInputVector;
input4 = acado.MexInputMatrix;


% Parameters with exemplary values. These are set/overwritten at runtime.
t_start = 0.0;     % Initial time [s]
t_end = 0.2;       % Time horizon [s]
dt = 0.1;          % Discretization time [s]
N = round(t_end/dt);  % Number of nodes
g_z = 9.8066;      % Gravity is everywhere [m/s^2]
w_max_yaw = 1;     % Maximal yaw rate [rad/s]
w_max_xy = 3;      % Maximal pitch and roll rate [rad/s]
T_min = 2;         % Minimal thrust [N]
T_max = 20;        % Maximal thrust [N]

% Bias to prevent division by zero.
epsilon = 0.1;     % Camera projection recover bias [m]

% System Dynamics
f.add(dot(p_x) ==  v_x);
f.add(dot(p_y) ==  v_y);
f.add(dot(p_z) ==  v_z);
f.add(dot(q_w) ==  0.5 * ( - w_x * q_x - w_y * q_y - w_z * q_z));
f.add(dot(q_x) ==  0.5 * ( w_x * q_w + w_z * q_y - w_y * q_z));
f.add(dot(q_y) ==  0.5 * ( w_y * q_w - w_z * q_x + w_x * q_z));
f.add(dot(q_z) ==  0.5 * ( w_z * q_w + w_y * q_x + w_z * q_y));
f.add(dot(v_x) ==  2 * ( q_w * q_y + q_x * q_z ) * T);
f.add(dot(v_y) ==  2 * ( q_y * q_z - q_w * q_x ) * T);
f.add(dot(v_z) ==  ( 1 - 2 * q_x * q_x - 2 * q_y * q_y ) * T - g_z);


h = {p_x, p_y, p_z,...
     q_w, q_x, q_y, q_z,...
     v_x, v_y, v_z,...
     T, w_x, w_y, w_z};
 
 % for final(Nth) step, does not include control terms
 hN = {p_x, p_y, p_z,...
       q_w, q_x, q_y, q_z,...
       v_x, v_y, v_z};

% reference 
r = input2;  

rN = input2;

% weights (variables in same order as h )

Q = diag([100;100;100;...
          100;100;100;100;...
          10;10;10;
          1;1;1;1]);

QN = Q(1:10,1:10);    % control not included
        
  
ocp = acado.OCP( t_start, t_end, N );


ocp.minimizeLSQ( Q, h, r );
ocp.minimizeLSQEndTerm( QN, hN, rN );


% Add system dynamics
ocp.subjectTo( f );
% Add constraints
ocp.subjectTo(-w_max_xy <= w_x <= w_max_xy);
ocp.subjectTo(-w_max_xy <= w_y <= w_max_xy);
ocp.subjectTo(-w_max_yaw <= w_z <= w_max_yaw);
ocp.subjectTo(T_min <= T <= T_max);

ocp.setNOD(10);


% %     % Define an algorithm to solve it.
%     algorithm = acado.OptimizationAlgorithm(ocp);
%     algorithm.set( 'INTEGRATOR_TOLERANCE', 1e-6 );
%     algorithm.set( 'KKT_TOLERANCE', 1e-3 );
%     algorithm.set('HESSIAN_APPROXIMATION',  'GAUSS_NEWTON');        % is robust, stable
%     algorithm.set('DISCRETIZATION_TYPE',    'MULTIPLE_SHOOTING');   % good convergence

% SETTING UP THE (SIMULATED) PROCESS
identity = acado.OutputFcn();
dynamicSystem = acado.DynamicSystem(f, identity);
% A dynamic system with the differential equation and an output function
process = acado.Process(dynamicSystem, 'INT_RK45');
% Simulate proces based on a dynamic model.

% SETTING UP THE MPC CONTROLLER:
algo = acado.RealTimeAlgorithm(ocp, 0.1);
% The class RealTimeAlgorithm serves as a user?interface to formulate and solve model predictive control problems.
algo.set('MAX_NUM_ITERATIONS', 3);
% Set some algorithm parameters

algo.set( 'INTEGRATOR_TOLERANCE', 1e-6 );
algo.set( 'KKT_TOLERANCE', 1e-3 );
algo.set('HESSIAN_APPROXIMATION',  'GAUSS_NEWTON');        % is robust, stable
algo.set('DISCRETIZATION_TYPE',    'MULTIPLE_SHOOTING');   % good convergence

% ref = acado.PeriodicReferenceTrajectory([0,[0,0,2,0.99,0,0,0,0,0,0,g_z,0,0,0];...
%                                          1,[0,0,2,0.99,0,0,0,0,0,0,g_z,0,0,0]])
ref = acado.StaticReferenceTrajectory(input4);
% Static reference trajectory that the ControlLaw aims to track.

controller = acado.Controller( algo,ref );
% Online control law for obtaining the control inputs of a process

% SETTING UP THE SIMULATION ENVIRONMENT, RUN THE EXAMPLE..
sim = acado.SimulationEnvironment( 0.0,5.0,process,controller );

sim.init(input3);
% 
% controller.init(input1, input3);
% controller.step(input1, input3);

END_ACADO

out = quad_mpc_4_RUN(1,[0,0,2,1,0,0,0,0,0,0,10,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,2,1,0,0,0,0,0,0,10,0,0,0])