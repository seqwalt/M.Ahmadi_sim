% Ahmadi_sim.m replicates the results presented in a paper by
% M. Ahmadi et al. titled "Safe Controller Synthesis for
% Data-Driven Differential Inclusions"

% Updated 7/08/20 by Sequoyah Walters

%% --- Parameter Panel --- %%
T_max = 160000;     % max thrust (N)
T = 0.2*T_max;      % Thrust
T_no_u = 0;         % Thrust with no input
A = -.9;            % Angle of attack (use -0.9 for engine failure)
L = 1;              % Lift starts at 100%
%si = [-0.07, 0, 0.07];      % side information +/- 10
si = 0;
len_si = length(si);
y0 = [95 2.2 120];  % initial values of velocity (m/s), flight path angle (deg), altitude (m)
failure = 'engine'; % System failure type: options are 'wing' or 'engine'

step = 0.05;          % Time step size (s) for solving ODE
tspan = 0:step:200;   % Time span to solve over
t_f = 0.5;            % Time of system failure (s) (should be a multiple of step) use 0.5 for engine failure
dur = 14.5;           % Duration of data measurement after failure (s) (multiple of step) use 14.5 for engine failure
ex_dur = 100;         % Extrapolation time after data stops

%% --- Aircraft Flight and Failure --- %%

% Solve until t_f (system failure) %
[t,y] = ode45(@(t,y) true_sys(y,T,A,L), tspan(1:find(tspan==t_f)), y0); 
S = [t(1:end-1),y(1:end-1,:)]; % Store history of states in S (except last state)
y_end = y(end,:);              % Last state will be used as initial cond. for next ode45

% System Failure %
if strcmp(failure,'wing') == 1
    % Wing failure:   100% Lift -> 20% Lift
    L = 0;
elseif strcmp(failure,'engine') == 1
    % Engine failure: 20% T_max -> 80% T_max (interpreted as 60% max thrust when no input is given)
    T = 0.8*T_max;
    T_no_u = 0.6*T_max;
end
[t,y] = ode45(@(t,y) true_sys(y,T_no_u,0,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end);
S_no_u = [t,y];
[t,y] = ode45(@(t,y) true_sys(y,T,0,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end); 
S_u1 = [t,y];  % u1 corresponds to thrust. (Ang of attack is 0)
[t,y] = ode45(@(t,y) true_sys(y,T_no_u,A,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end); 
S_u2 = [t,y];  % u2 corresponds to ang of attack. (Thrust input is 0)

% --- True flight system --- %
[t,y] = ode45(@(t,y) true_sys(y,T,A,L), tspan(find(tspan==t_f):find(tspan==t_f+dur)), y_end);
S_true_pre = [t,y];
S_true = [S;t,y];
y_end_data = y(end,:);
    % The rest of the flight if there were no controller
[t,y] = ode45(@(t,y) true_sys(y,T,A,L), tspan(find(tspan==t_f+dur):end), y_end_data);
S_true = [S_true;t,y];

times = S_no_u(:,1);  % all times are the same for S_no_u, S_u1 and S_u2

%% --- Approximate System Model --- %%

u1 = T; u2 = A;

% Approximation for Velocity
[V1_no_u,V2_no_u,V1_u1,V2_u2] = control_affine_spline('v',times,ex_dur,step,S_true_pre,S_u1,S_u2,S_no_u,u1,u2); % 'v' for velocity

% Approximation for Flight Path Angle
[P1_no_u,P2_no_u,P1_u1,P2_u2] = control_affine_spline('f',times,ex_dur,step,S_true_pre,S_u1,S_u2,S_no_u,u1,u2); % 'f' for FPA

% Approximation for Altitude
[A1_no_u,A2_no_u,A1_u1,A2_u2] = control_affine_spline('a',times,ex_dur,step,S_true_pre,S_u1,S_u2,S_no_u,u1,u2); % 'a' for altitude

% Combine Splines
system = [V1_no_u,V2_no_u,V1_u1,V2_u2; P1_no_u,P2_no_u,P1_u1,P2_u2; A1_no_u,A2_no_u,A1_u1,A2_u2];

% SYSTEM MODEL
eval_times = (times(end)+step:step:times(end)+step+ex_dur)';
S_approx = cell(len_si,len_si);
for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        [t,x] = ode45(@(t,x) approx_sys(t,system,u1,u2,si(si1_ind),si(si2_ind)), eval_times, y_end_data); 
        S_approx(si1_ind,si2_ind) = {[t,x]};
    end
end

%% --- Controller Synthesis and Usage --- %%

% Find B (Lyapunov Function)

% Find W (Positive definite function)

%% --- Remove Negative Altitude Data --- %%
for si1_ind = 1:len_si
    for si2_ind = 1:len_si
        S_approx{si1_ind,si2_ind} = above_zero_alt(S_approx{si1_ind,si2_ind},'approximate');
    end
end
S_true   = above_zero_alt(S_true,  '       true');

%% --- Plot --- %%
figure
hold on
axis tight
fill([86,86,100,100],[-10,3.3,3.3,-10],[1 0.8 0.8])  % Make the unsafe region pink
plot3(S_true(:,2),S_true(:,3),S_true(:,4),'k')

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
plot3(S_approx{1,1}(1,2),S_approx{1,1}(1,3),S_approx{1,1}(1,4),'k^') % data record stops

xlabel('Velocity (m/s)')
ylabel('Flight path angle (deg)')
zlabel('Altitude (m)')
legend('Unsafe set','True System','Approx System','Start','Failure','Approx Starts','Location','northwest')
grid on