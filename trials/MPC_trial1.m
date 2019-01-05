%% Trial 1 quad MPC

clear;

BEGIN_ACADO;

    acadoSet('problemname','MPC_trial1');

    DifferentialState   p_x p_y p_z;
    DifferentialState   q_w q_x q_y q_z;
    DifferentialState   v_x v_y v_z;
    Control             T w_x w_y w_z;

    input_time = acado.MexInput;
    input_ref  = acado.MexInputVector;
    input_x    = acado.MexInputVector;
    input_sim  = acado.MexInputMatrix;

    % Parameters
    t_start = 0.0;     % Initial time [s]
    t_end = 1.0;       % Time horizon [s]
    dt = 0.1;          % Discretization time [s]
    N = round(t_end/dt);  % Number of nodes
    g_z = 9.8066;      % Gravity is everywhere [m/s^2]
    w_max_yaw = 1;     % Maximal yaw rate [rad/s]
    w_max_xy = 3;      % Maximal pitch and roll rate [rad/s]
    T_min = 2;         % Minimal thrust [N]
    T_max = 20;        % Maximal thrust [N]

    % System Dynamics
    f = acado.DifferentialEquation();

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

    %% Optimization Problem
    ocp = acado.OCP(t_start,t_end,N);

    h = {p_x, p_y, p_z,...
         q_w, q_x, q_y, q_z,...
         v_x, v_y, v_z,...
         T, w_x, w_y, w_z};
    % for final(Nth) step, does not include control terms

    hN = {p_x, p_y, p_z,...
          q_w, q_x, q_y, q_z,...
          v_x, v_y, v_z};

    r = input_ref;
    rN = input_ref;

    % weights (variables in same order as h )

    Q = diag([100;100;100;...
              100;100;100;100;...
              10;10;10;
              1;1;1;1]);
    QN = Q(1:10,1:10);    % control not included

    ocp.minimizeLSQ( Q, h, r );
    ocp.minimizeLSQEndTerm( QN, hN, rN );

    % Add system dynamics
    ocp.subjectTo( f );

    % Add constraints
    ocp.subjectTo('AT_START',p_x == 0);
    ocp.subjectTo('AT_START',p_y == 0);
    ocp.subjectTo('AT_START',q_w == 0);
    ocp.subjectTo('AT_START',q_x == 0);
    ocp.subjectTo('AT_START',q_y == 0);
    ocp.subjectTo('AT_START',q_z == 0);
    ocp.subjectTo('AT_START',v_x == 0);
    ocp.subjectTo('AT_START',v_y == 0);
    ocp.subjectTo('AT_START',v_z == 0);
    
    ocp.subjectTo(-w_max_xy <= w_x <= w_max_xy);
    ocp.subjectTo(-w_max_xy <= w_y <= w_max_xy);
    ocp.subjectTo(-w_max_yaw <= w_z <= w_max_yaw);
    ocp.subjectTo(T_min <= T <= T_max);

    %% Setting up the MPC Controller

    algo = acado.RealTimeAlgorithm(ocp, 0.1);
    algo.set('MAX_NUM_ITERATIONS', 3);  % number of optimization iterations per cycle

    % Set some algorithm parameters
    algo.set( 'INTEGRATOR_TOLERANCE', 1e-6 );
    algo.set( 'KKT_TOLERANCE', 1e-3 );
    algo.set('HESSIAN_APPROXIMATION',  'GAUSS_NEWTON');        % is robust, stable
    algo.set('DISCRETIZATION_TYPE',    'MULTIPLE_SHOOTING');   % good convergence
    %algo.initializeDifferentialStates([0,0,0,0,1,0,0,0,0,0,0])
    ref = acado.StaticReferenceTrajectory(input_sim);

    controller = acado.Controller( algo,ref );

    % SETTING UP THE (SIMULATED) PROCESS
    identity = acado.OutputFcn();
    dynamicSystem = acado.DynamicSystem(f, identity);
    process = acado.Process(dynamicSystem, 'INT_RK45');
    % SETTING UP THE SIMULATION ENVIRONMENT, RUN THE EXAMPLE..
    sim = acado.SimulationEnvironment( 0.0,5.0,process,controller );
    sim.init(input_x);
   
 END_ACADO;

 out = MPC_trial1_RUN(1,[0,0,2,1,0,0,0,0,0,0,10,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,2,1,0,0,0,0,0,0,10,0,0,0])
 figure;
 plot3(out.STATES(:,1),out.STATES(:,2),out.STATES(:,3))

 
