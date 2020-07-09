%% Trajectory stuff
x= 2;
T_max = 160000;     % max thrust (N)
T = 0.2*T_max;      % Thrust
T_no_u = 0;         % Thrust with no input
A = 20;            % Angle of attack (use -0.9 for engine failure)
L = 1;              % Lift starts at 100%
si = [-0.07, 0, 0.07];      % side information +/- 10
%si = 0;
len_si = length(si);
y0 = [95 2.2 120];  % initial values of velocity (m/s), flight path angle (deg), altitude (m)
failure = 'engine'; % System failure type: options are 'wing' or 'engine'

step = 0.05;          % Time step size (s) for solving ODE
tspan = 0:step:200;   % Time span to solve over
t_f = 30;             % Time of system failure (s) (should be a multiple of step) use 0.5 for engine failure
dur = 1;             % Duration of data measurement after failure (s) (multiple of step) use 14.5 for engine failure
ex_dur = 100;         % Extrapolation time after data stops

[t,y] = ode45(@(t,y) true_sys(y,T,A,L), tspan(1:find(tspan==t_f)), y0); 
S = [t(1:end-1),y(1:end-1,:)]; % Store history of states in S (except last state)
y_end = y(end,:);              % Last state will be used as initial cond. for next ode45

[t,y] = ode45(@(t,y) true_sys(y,0.6*T_max+T,A,L), tspan(1:find(tspan==t_f)), y0); 
SS = [t(1:end-1),y(1:end-1,:)]; % Store history of states in S (except last state)
Sy_end = y(end,:);              % Last state will be used as initial cond. for next ode45

A = 1;
[t,y] = ode45(@(t,y) true_sys(y,T,A,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end); 
S2 = [t,y]; % Store history of states in S (except last state)

[t,y] = ode45(@(t,y) true_sys(y,0.6*T_max+T,A,L), tspan(find(tspan==t_f):find(tspan==(t_f+dur))), Sy_end); 
SS2 = [t,y]; % Store history of states in S (except last state)

figure
hold on
axis tight
plot(S(:,1),S(:,x),'b',S2(:,1),S2(:,x),'k');
plot(SS(:,1),SS(:,x),'k',SS2(:,1),SS2(:,x),'b');
figure
axis tight
plot(S(:,1),SS(:,x)-S(:,x),'r',S2(:,1),SS2(:,x)-S2(:,x),'m')

