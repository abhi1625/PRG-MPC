
digits(3)
% global iter
% global tt yy uu

params.g = 10;
params.Q = diag([100,100,100,1,1,1,1,10,10,10]);
params.R = diag([1,5,5,0.1]);

params.u0 = [0;0;0;10];
params.x0 = [0;0;0;1;0;0;0;0;0;0];
params.xr = [0;0;2;1;0;0;0;0;0;0;9.8;0;0;0];

% uu = [params.u0.'];

y0 = params.x0;
dyn = @(t,y)quad_dyn_mpc(t,y,params);
op = @(t,y,flag)idlk(t,y,flag,params);
options = odeset('OutputFcn',op);

[t,y] = ode45(dyn,[0 10],y0);%,options);