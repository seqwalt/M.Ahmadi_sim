m = 60000;   % mass is 60,000 kg
grav = 9.8;  % accel due to grav is 9.8 ms^(-2)
v_safe = 86;
M = 0;

t_N   = t_f;                   % time t_N
v_t_N = S_true(S_true==t_N,2); % velocity @ time t_N
p_t_N = S_true(S_true==t_N,3); % flight path angle @ time t_N
a_t_N = S_true(S_true==t_N,4); % altitude @ time t_N

% T defined as t_f + t_dur + wait ... time to start controller
v_T = y_control(1); % velocity @ time T
p_T = y_control(2); % flight path angle @ time T
a_T = y_control(3); % altitude @ time T

syms v p a t c % vel, fpa, alt, time, and constant
vars = [v; p; a; t];

%% 1) Find a reasonable candidate for B

% Initialize the sum of squares program
prog1 = sosprogram(vars);
% Declare variables and functions
[prog1,B(v,p,a,t)]  = sospolyvar(prog1,monomials(vars,0:4),'wscoeff');
[prog1,W(v,p,a,t)]  = sospolyvar(prog1,monomials(vars,0:4),'wscoeff');
[prog1,m1(v,p,a,t)] = sospolyvar(prog1,monomials(vars,0:2),'wscoeff');
[prog1,m2(v,p,a,t)] = sospolyvar(prog1,monomials(vars,0:2),'wscoeff');
[prog1,s(v,p,a,t)]  = sospolyvar(prog1,monomials(vars,0:2),'wscoeff');
[prog1,q2(v,p,a,t)] = sospolyvar(prog1,monomials(vars,0:2),'wscoeff');
prog1 = sosdecvar(prog1,c);
dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% Define SOSP constraints
prog1 = sosineq(prog1,W(v,p,a,t));
prog1 = sosineq(prog1,m1(v,p,a,t));
prog1 = sosineq(prog1,m2(v,p,a,t));
prog1 = sosineq(prog1,q2(v,p,a,t));
prog1 = sosineq(prog1,s(v,p,a,t));
prog1 = sosineq(prog1,c);

expr1 = B(v_T,p_T,a_T,T) - B(v_t_N,p_t_N,a_t_N,t_N) + s(v,p,a,t)*(v_safe-v) - c;
expr2 = -(diff(B(v,p,a,t),t) + dBdx(v,p,a,t).'*(f(t) - M) + W(v,p,a,t) + m1(v,p,a,t)*(t-t_N)*(t-T));
expr3 = -(diff(B(v,p,a,t),t) + dBdx(v,p,a,t).'*(f(t) + M) + W(v,p,a,t) + m2(v,p,a,t)*(t-t_N)*(t-T));
prog1 = sosineq(prog1,expr1);
prog1 = sosineq(prog1,expr2);
prog1 = sosineq(prog1,expr3);

% Constrain outputs --> 0 <= Th <= 15;  -3 <= fpa <= 15
% u = (g(t).'*dBdx(v,p,a,t)*W(v,p,a,t))/(dBdx(v,p,a,t).'*(g(t)*g(t).')*dBdx(v,p,a,t));
% prog1 = sosineq(prog1,u(1));
% prog1 = sosineq(prog1,Th_max - u(1));
% prog1 = sosineq(prog1,u(2) + 3);
% prog1 =  sosineq(prog1,15 - u(2));

% Solve for B and q2
solver_opt.solver = 'sedumi';
prog1 = sossolve(prog1,solver_opt);
B = sosgetsol(prog1,B);
q2 = sosgetsol(prog1,q2);

%% 2) Use the previously determined B to find q1

% Initialize the sum of squares program
prog2 = sosprogram(vars);
% Declare variables and functions
[prog2,q1(v,p,a,t)] = sospolyvar(prog2,monomials(vars,0:2),'wscoeff');
dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% Define SOSP constraints
expr4 = -q1(v,p,a,t)*(dBdx(v,p,a,t).'*g(t)) - q2(v,p,a,t);expr4_1 = expr4(1);
expr4_2 = expr4(2);
prog2 = sosineq(prog2,expr4_1);
prog2 = sosineq(prog2,expr4_2);
% Solve for q1
solver_opt.solver = 'sedumi';
prog2 = sossolve(prog2,solver_opt);
q1 = sosgetsol(prog2,q1);

%% 3) Use the previously found q1 to find a better B

% Initialize the sum of squares program
prog_fin = sosprogram(vars);
% Declare variables and functions
[prog_fin,B(v,p,a,t)]  = sospolyvar(prog_fin,monomials(vars,0:4),'wscoeff');
[prog_fin,W(v,p,a,t)]  = sospolyvar(prog_fin,monomials(vars,0:4),'wscoeff');
[prog_fin,m1(v,p,a,t)] = sospolyvar(prog_fin,monomials(vars,0:2),'wscoeff');
[prog_fin,m2(v,p,a,t)] = sospolyvar(prog_fin,monomials(vars,0:2),'wscoeff');
[prog_fin,s(v,p,a,t)]  = sospolyvar(prog_fin,monomials(vars,0:2),'wscoeff');
[prog_fin,q2(v,p,a,t)] = sospolyvar(prog_fin,monomials(vars,0:2),'wscoeff');
prog_fin = sosdecvar(prog_fin,c);
dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% Define SOSP constraints
prog_fin = sosineq(prog_fin,W(v,p,a,t));
prog_fin = sosineq(prog_fin,m1(v,p,a,t));
prog_fin = sosineq(prog_fin,m2(v,p,a,t));
prog_fin = sosineq(prog_fin,q2(v,p,a,t));
prog_fin = sosineq(prog_fin,s(v,p,a,t));
prog_fin = sosineq(prog_fin,c);

expr1 = B(v_T,p_T,a_T,T) - B(v_t_N,p_t_N,a_t_N,t_N) + s(v,p,a,t)*(v_safe-v) - c;
expr2 = -(diff(B(v,p,a,t),t) + dBdx(v,p,a,t).'*(f(t) - M) + W(v,p,a,t) + m1(v,p,a,t)*(t-t_N)*(t-T));
expr3 = -(diff(B(v,p,a,t),t) + dBdx(v,p,a,t).'*(f(t) + M) + W(v,p,a,t) + m2(v,p,a,t)*(t-t_N)*(t-T));
expr4 = -q1(v,p,a,t)*(dBdx(v,p,a,t).'*g(t)) - q2(v,p,a,t);
expr4_1 = expr4(1);
expr4_2 = expr4(2);

prog_fin = sosineq(prog_fin,expr1);
prog_fin = sosineq(prog_fin,expr2);
prog_fin = sosineq(prog_fin,expr3);
prog_fin = sosineq(prog_fin,expr4_1);
prog_fin = sosineq(prog_fin,expr4_2);

% Constrain outputs --> 0 <= Th <= 15;  -3 <= fpa <= 15
% u = (g(t).'*dBdx(v,p,a,t)*W(v,p,a,t))/(dBdx(v,p,a,t).'*(g(t)*g(t).')*dBdx(v,p,a,t));
% prog_fin = sosineq(prog_fin,u(1));
% prog_fin = sosineq(prog_fin,Th_max - u(1));
% prog_fin = sosineq(prog_fin,u(2) + 3);
% prog_fin =  sosineq(prog_fin,15 - u(2));

% Solve for B and W
solver_opt.solver = 'sedumi';
echo on;
prog_fin = sossolve(prog_fin,solver_opt);
echo off;
SOL_B = sosgetsol(prog_fin,B);
SOL_W = sosgetsol(prog_fin,W);

%% 4) Synthesize the controller
dBdx(v,p,a,t) = [diff(SOL_B,v); diff(SOL_B,p); diff(SOL_B,a)];
dBdx = dBdx(v,p,a,t);
u = (g(t).'*dBdx*SOL_W)/(dBdx.'*(g(t)*g(t).')*dBdx);
