% Ahmadi_sim.m replicates the results presented in a paper by
% M. Ahmadi et al. titled "Safe Controller Synthesis for
% Data-Driven Differential Inclusions"
% Updated 7/08/20 by Sequoyah Walters

%% --- Parameter Panel --- %%
Th_max = 160000;     % max thrust (N)
Th = 0.2*Th_max;      % Thrust
Th_no_u = 0;         % Thrust with no input
A = -.9;            % Angle of attack (use -0.9 for engine failure)
L = 1;              % Lift starts at 100%
u1a = 0.3*Th;        % 1 corresponds to thrust
u2a = 0.8*A;        % 2 corresponds to ang of attack
u1b = Th;            % b corresponds to final trajectory in the 3-traj approx
u2b = A;
%si = [0, -0.07, 0.07];    % side information
si = 0;
len_si = length(si);
y0 = [95 2.2 120];    % initial values of velocity (m/s), flight path angle (deg), altitude (m)
failure = 'engine';   % System failure type: options are 'wing' or 'engine'

step = 0.05;          % Time step size (s) for solving ODE
end_time = 200;
tspan = 0:step:end_time;   % Time span to solve over
t_f = 0.5;            % Time of system failure (s) (should be a multiple of step) use 0.5 for engine failure
dur = 14.5;           % Duration of data measurement after failure (s) (multiple of step) use 14.5 for engine failure
wait = 1; % wait 1 sec after approx is started to engage controller
ex_dur = 100;         % Extrapolation time after data stops
k = 3;                % Spline order
ex_k = 3;             % Order of extrapolation

%% --- Aircraft Flight and Failure --- %%

% Solve until t_f (system failure) %
[t,y] = ode45(@(t,y) true_sys(y,Th,A,L), tspan(1:find(tspan==t_f)), y0); 
S = [t(1:end-1),y(1:end-1,:)]; % Store history of states in S (except last state)
y_end = y(end,:);              % Last state will be used as initial cond. for next ode45

% System Failure %
if strcmp(failure,'wing') == 1
    % Wing failure:   100% Lift -> 20% Lift
    L = 0;
elseif strcmp(failure,'engine') == 1
    % Engine failure: 20% T_max -> 80% T_max (interpreted as 60% max thrust when no input is given)
    Th_no_u = 0.6*Th_max;
end
[t,y] = ode45(@(t,y) true_sys(y,Th_no_u,0,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end);
S_no_u = [t,y];
[t,y] = ode45(@(t,y) true_sys(y,u1a + Th_no_u,u2a,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end); 
S_a = [t,y];

% --- True flight system --- %
[t,y] = ode45(@(t,y) true_sys(y,u1b + Th_no_u,u2b,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end);
S_true_pre = [t,y];
S_true = [S;t,y];
y_end_data = y(end,:);
% The rest of the flight if there were no controller
[t,y] = ode45(@(t,y) true_sys(y,u1b + Th_no_u,u2b,L), tspan(find(tspan==t_f+dur):end), y_end_data);
S_true_no = [S_true;t,y];
% Flight until controller
T   = t_f + dur + wait;    % time to start controller
[t,y] = ode45(@(t,y) true_sys(y,u1b + Th_no_u,u2b,L), tspan(find(tspan==t_f+dur):find(tspan==T)), y_end_data);
S_true_c = [S_true;t,y];
S_true_c = S_true_c(1:end-1,:);
y_control = y(end,:);

times = S_no_u(:,1);  % all times are the same for S_no_u, S_u1 and S_u2
ex_times = times(1):step:times(end) + ex_dur; % Time span plus extrap time 

%% --- Approximate System Model --- %%

% Coefficients for Velocity
[V_C_coefs,V_u_coefs] = spline_coefs('v',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'v' for velocity

% Coefficients for Flight Path Angle
[P_C_coefs,P_u_coefs] = spline_coefs('f',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'f' for FPA

% Coefficients for Altitude
[A_C_coefs,A_u_coefs] = spline_coefs('a',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'a' for altitude

system_coefs = {V_C_coefs,V_u_coefs; P_C_coefs,P_u_coefs; A_C_coefs,A_u_coefs};

% SYSTEM MODEL
eval_times = (times(end)+step:step:times(end)+step+ex_dur)';
S_approx = cell(len_si,len_si);
for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        [t,x] = ode45(@(t,x) approx_sys(t,ex_times,k,system_coefs,u1b,u2b,si(si1_ind),si(si2_ind)), eval_times, y_end_data); 
        S_approx(si1_ind,si2_ind) = {[t,x]};
    end
end

%% --- Controller Synthesis and Usage --- %%

m = 60000;  % mass is 60,000 kg
grav = 9.8;    % accel due to grav is 9.8 ms^(-2)
v_safe = 86;
M = 10;
PI = 3.14159;

t_N   = t_f;                   % time t_N
v_t_N = S_true(S_true==t_N,2); % velocity @ time t_N
p_t_N = S_true(S_true==t_N,3); % flight path angle @ time t_N
a_t_N = S_true(S_true==t_N,4); % altitude @ time t_N

% T defined above as t_f + t_dur + wait ... time to start controller
v_T = S_true_c(S_true_c==T,2); % velocity @ time T
p_T = S_true_c(S_true_c==T,3); % flight path angle @ time T
a_T = S_true_c(S_true_c==T,4); % altitude @ time T

syms v p a t c % vel, fpa, alt, time, and constant
vars = [v; p; a; t];
% Constructing the vector field dx/dt = f + g*u

f = [(-7.5125*v^2/m) - (m*grav*p*PI/180);
     (85.75*v/m) - (grav/v);
     (v*p*PI/180)];

g = [1/m , -32.34*v^2/m;
     0   , (PI*Th/(180*m*v)) + (288.12*v/m);
     0   , 0];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);
% =============================================
% Declare variables and functions
[prog,B(v,p,a,t)]  = sospolyvar(prog,monomials(vars,[0:4]),'wscoeff');
[prog,W(v,p,a,t)]  = sospolyvar(prog,monomials(vars,[0:4]),'wscoeff');
[prog,m1(v,p,a,t)] = sospolyvar(prog,monomials(vars,[0:2]),'wscoeff');
[prog,m2(v,p,a,t)] = sospolyvar(prog,monomials(vars,[0:2]),'wscoeff');
[prog,s(v,p,a,t)]  = sospolyvar(prog,monomials(vars,[0:2]),'wscoeff');
[prog,q2(v,p,a,t)] = sospolyvar(prog,monomials(vars,[0:2]),'wscoeff');
prog = sosdecvar(prog,[c]);
% =============================================
% Next, define SOSP constraints
prog = sosineq(prog,W(v,p,a,t));
prog = sosineq(prog,m1(v,p,a,t));
prog = sosineq(prog,m2(v,p,a,t));
prog = sosineq(prog,q2(v,p,a,t));
prog = sosineq(prog,s(v,p,a,t));
prog = sosineq(prog,c);

q1(v,p,a,t) = v^2 + p^2 + a^2 +t^2;

expr1 = B(v_T,p_T,a_T,T) - B(v_t_N,p_t_N,a_t_N,t_N) + s(v,p,a,t)*(v_safe-v) - c;
expr2 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)-M) + diff(B(v,p,a,t),p)*(f(2)-M) + diff(B(v,p,a,t),a)*(f(3)-M) + W(v,p,a,t) + m1(v,p,a,t)*(t-t_N)*(t-T))*v;
expr3 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)+M) + diff(B(v,p,a,t),p)*(f(2)+M) + diff(B(v,p,a,t),a)*(f(3)+M) + W(v,p,a,t) + m2(v,p,a,t)*(t-t_N)*(t-T))*v;
expr4 = -q1(v,p,a,t)*([diff(B(v,p,a,t),v), diff(B(v,p,a,t),p), diff(B(v,p,a,t),a)]*g) - q2(v,p,a,t);
expr4_1 = expr4(1);
expr4_2 = expr4(2)*v;

prog = sosineq(prog,expr1);
prog = sosineq(prog,expr2);
prog = sosineq(prog,expr3);
prog = sosineq(prog,expr4_1);
prog = sosineq(prog,expr4_2);

% Constrain outputs --> 0 <= Th <= 15;  -3 <= fpa <= 15
% dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% u = (g.'*dBdx(v,p,a,t)*W(v,p,a,t))/(dBdx(v,p,a,t).'*(g*g.')*dBdx(v,p,a,t));
% 
% prog = sosineq(prog,u(1));
% prog = sosineq(prog,Th_max - u(1));
% prog = sosineq(prog,u(2) + 3);
% prog = sosineq(prog,15 - u(2));
% =============================================
% And call solver
solver_opt.solver = 'sedumi';
echo on;
prog = sossolve(prog,solver_opt);
echo off;
% =============================================
% Finally, get solution
SOL_B = sosgetsol(prog,B);
SOL_W = sosgetsol(prog,W);

% Synthesize Controller
dBdx(v,p,a,t) = [diff(SOL_B,v); diff(SOL_B,p); diff(SOL_B,a)];
dBdx = dBdx(v,p,a,t);
u = (g.'*dBdx*SOL_W)/(dBdx.'*(g*g.')*dBdx);

% Use Controller
for TT = T:2*step:end_time/2
    V = y_control(1);
    P = y_control(2);
    A = y_control(3);
    U = u(V,P,A,TT);
    TT_nxt = TT+2*step;
    [t,y] = ode45(@(t,y) true_sys(y,U(1) + Th_no_u,U(2),L), TT:step:TT_nxt, y_control);
    t = [t(1);t(2)]; % remove TT+2*step data
    y = [y(1,:);y(2,:)];
    disp(y);
    S_true_c = [S_true_c;t,y];
    y_control = y(end,:);
end

%% --- Remove Negative Altitude Data --- %%
for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        S_approx{si1_ind,si2_ind} = above_zero_alt(S_approx{si1_ind,si2_ind},'approximate');
    end
end
S_true_no   = above_zero_alt(S_true_no,  '       true_no');
S_true_c   = above_zero_alt(S_true_c,  '       true_c');

%% --- Plot --- %%
figure
hold on
axis tight
fill([86,86,100,100],[-10,3.3,3.3,-10],[1 0.8 0.8])  % Make the unsafe region pink
%plot3(S_true_no(:,2),S_true_no(:,3),S_true_no(:,4),'k')
plot3(S_true_c(:,2),S_true_c(:,3),S_true_c(:,4),'b')

for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        color = 'b';
        if si(si1_ind) == 0 && si(si2_ind)==0
            color = 'm';
        end
        plot3(S_approx{si1_ind,si2_ind}(:,2),S_approx{si1_ind,si2_ind}(:,3),S_approx{si1_ind,si2_ind}(:,4),color)
    end
end

plot3(S_true(1,2),S_true(1,3),S_true(1,4),'ko') % start
plot3(S_true(S_true==t_f,2),S_true(S_true==t_f,3),S_true(S_true==t_f,4),'k*') % failure
plot3(S_approx{1,1}(1,2),S_approx{1,1}(1,3),S_approx{1,1}(1,4),'k^') % Start of approximation

xlabel('Velocity (m/s)')
ylabel('Flight path angle (deg)')
zlabel('Altitude (m)')
legend('Unsafe set','True System','Approx System','Start','Failure','Approx Starts','Location','northwest')
grid on
