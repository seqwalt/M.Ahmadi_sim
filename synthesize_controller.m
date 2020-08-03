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
v_T = y_control(1); % velocity @ time T
p_T = y_control(2); % flight path angle @ time T
a_T = y_control(3); % altitude @ time T

syms v p a t c % vel, fpa, alt, time, and constant
vars = [v; p; a; t];
% Constructing the vector field dx/dt = f + g*u

f = [(-7.5125*v^2/m) - (m*grav*p*PI/180);
     (85.75*v/m) - (grav/v);
     (v*p*PI/180)];

g = [1/m , -32.34*v^2/m;
     0   , (PI*Th/(180*m*v)) + (288.12*v/m);
     0   , 0];

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
% Define SOSP constraints
prog1 = sosineq(prog1,W(v,p,a,t));
prog1 = sosineq(prog1,m1(v,p,a,t));
prog1 = sosineq(prog1,m2(v,p,a,t));
prog1 = sosineq(prog1,q2(v,p,a,t));
prog1 = sosineq(prog1,s(v,p,a,t));
prog1 = sosineq(prog1,c);

expr1 = B(v_T,p_T,a_T,T) - B(v_t_N,p_t_N,a_t_N,t_N) + s(v,p,a,t)*(v_safe-v) - c;
expr2 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)-M) + diff(B(v,p,a,t),p)*(f(2)-M) + diff(B(v,p,a,t),a)*(f(3)-M) + W(v,p,a,t) + m1(v,p,a,t)*(t-t_N)*(t-T))*v; % need to multiply by v so there is no v in the denominator of SOSP
expr3 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)+M) + diff(B(v,p,a,t),p)*(f(2)+M) + diff(B(v,p,a,t),a)*(f(3)+M) + W(v,p,a,t) + m2(v,p,a,t)*(t-t_N)*(t-T))*v; % need to multiply by v so there is no v in the denominator of SOSP
prog1 = sosineq(prog1,expr1);
prog1 = sosineq(prog1,expr2);
prog1 = sosineq(prog1,expr3);

% Constrain outputs --> 0 <= Th <= 15;  -3 <= fpa <= 15
% dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% u = (g.'*dBdx(v,p,a,t)*W(v,p,a,t))/(dBdx(v,p,a,t).'*(g*g.')*dBdx(v,p,a,t));
% prog = sosineq(prog,u(1));
% prog = sosineq(prog,Th_max - u(1));
% prog = sosineq(prog,u(2) + 3);
% prog =  sosineq(prog,15 - u(2));

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
% Define SOSP constraints
expr4 = -q1(v,p,a,t)*([diff(B(v,p,a,t),v), diff(B(v,p,a,t),p), diff(B(v,p,a,t),a)]*g) - q2(v,p,a,t);
expr4_1 = expr4(1);
expr4_2 = expr4(2)*v;  % need to multiply by v so there is no v in the denominator of SOSP
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
% Define SOSP constraints
prog_fin = sosineq(prog_fin,W(v,p,a,t));
prog_fin = sosineq(prog_fin,m1(v,p,a,t));
prog_fin = sosineq(prog_fin,m2(v,p,a,t));
prog_fin = sosineq(prog_fin,q2(v,p,a,t));
prog_fin = sosineq(prog_fin,s(v,p,a,t));
prog_fin = sosineq(prog_fin,c);

expr1 = B(v_T,p_T,a_T,T) - B(v_t_N,p_t_N,a_t_N,t_N) + s(v,p,a,t)*(v_safe-v) - c;
expr2 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)-M) + diff(B(v,p,a,t),p)*(f(2)-M) + diff(B(v,p,a,t),a)*(f(3)-M) + W(v,p,a,t) + m1(v,p,a,t)*(t-t_N)*(t-T))*v; % need to multiply by v so there is no v in the denominator of SOSP
expr3 = -(diff(B(v,p,a,t),t) + diff(B(v,p,a,t),v)*(f(1)+M) + diff(B(v,p,a,t),p)*(f(2)+M) + diff(B(v,p,a,t),a)*(f(3)+M) + W(v,p,a,t) + m2(v,p,a,t)*(t-t_N)*(t-T))*v; % need to multiply by v so there is no v in the denominator of SOSP
expr4 = -q1(v,p,a,t)*([diff(B(v,p,a,t),v), diff(B(v,p,a,t),p), diff(B(v,p,a,t),a)]*g) - q2(v,p,a,t);
expr4_1 = expr4(1);
expr4_2 = expr4(2)*v;  % need to multiply by v so there is no v in the denominator of SOSP

prog_fin = sosineq(prog_fin,expr1);
prog_fin = sosineq(prog_fin,expr2);
prog_fin = sosineq(prog_fin,expr3);
prog_fin = sosineq(prog_fin,expr4_1);
prog_fin = sosineq(prog_fin,expr4_2);

% Constrain outputs --> 0 <= Th <= 15;  -3 <= fpa <= 15
% dBdx(v,p,a,t) = [diff(B(v,p,a,t),v); diff(B(v,p,a,t),p); diff(B(v,p,a,t),a)];
% u = (g.'*dBdx(v,p,a,t)*W(v,p,a,t))/(dBdx(v,p,a,t).'*(g*g.')*dBdx(v,p,a,t));
% prog = sosineq(prog,u(1));
% prog = sosineq(prog,Th_max - u(1));
% prog = sosineq(prog,u(2) + 3);
% prog =  sosineq(prog,15 - u(2));

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
u = (g.'*dBdx*SOL_W)/(dBdx.'*(g*g.')*dBdx);
