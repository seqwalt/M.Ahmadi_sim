% Ahmadi_sim.m replicates the results presented in a paper by
% M. Ahmadi et al. titled "Safe Controller Synthesis for
% Data-Driven Differential Inclusions"
% Updated 7/08/20 by Sequoyah Walters
tic
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
k = 5;                % Spline order
ex_k = 5;             % Order of extrapolation

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
y_control = y(end,:);

times = S_no_u(:,1);  % all times are the same for S_no_u, S_u1 and S_u2
ex_times = times(1):step:times(end) + dur; % Time span plus extrap time 

%% --- Approximate System Model --- %%

% Coefficients for Velocity
[V_C_coefs,V_u_coefs] = spline_coefs('vel',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'vel' for velocity

% Coefficients for Flight Path Angle
[P_C_coefs,P_u_coefs] = spline_coefs('fpa',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'fpa' for FPA

% Coefficients for Altitude
[A_C_coefs,A_u_coefs] = spline_coefs('alt',times,ex_times,k,ex_k,S_no_u,S_a,S_true_pre,u1a,u2a,u1b,u2b); % 'alt' for altitude

system_coefs = {V_C_coefs,V_u_coefs; P_C_coefs,P_u_coefs; A_C_coefs,A_u_coefs};
[f,g] = F_and_G(t_f+dur, ex_times, k, system_coefs);

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

% Sythesize the controller "u" using MATLAB file synthesize_controller
synthesize_controller

% Use Controller
TT = T;
U_hist = [];
while y_control(3) > 0
    TT_nxt = TT+2*step;
    V = y_control(1);
    P = y_control(2);
    A = y_control(3);
    U = double(u(V,P,A,TT)); % u(v,p,a,t) is created in synthesize_controller
    U1 = U(1); U2 = U(2);

    U_hist = [U_hist;U.'];
    [t,y] = ode45(@(t,y) true_sys(y,U1 + Th_no_u,U2,L), TT:step:TT_nxt, y_control);
    S_true_c = [S_true_c;t,y];
    S_true_c = S_true_c(1:end-1,:);
    y_control = y(end,:);
    TT = TT_nxt;
end

%% --- Remove Negative Altitude Data --- %%
for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        S_approx{si1_ind,si2_ind} = above_zero_alt(S_approx{si1_ind,si2_ind},'approximate');
    end
end
S_true_no   = above_zero_alt(S_true_no,  '       true_no_c');
S_true_c   = above_zero_alt(S_true_c,  '       true_c');

%% --- Plot --- %%
figure
hold on
axis tight
fill([86,86,110,110],[-65,3.3,3.3,-65],[1 0.8 0.8])  % Make the unsafe region pink
plot3(S_true_no(:,2),S_true_no(:,3),S_true_no(:,4),'k')
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
legend('Unsafe set','True: No Control','True: Control','Approx System','Start','Failure','Approx Starts','Location','northwest')
grid on
toc
system('say simulation complete')