echo on;
syms v p a t;
vartable = [v,p,a,t];
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable); % No decision variables.
% =============================================
% Next, define the inequality
test_func = v^2 + p^2 + a^2 +t^2;
prog = sosineq(prog,test_func);
% =============================================
% And call solver
solver_opt.solver = 'sedumi';
[prog,info] = sossolve(prog,solver_opt);
% =============================================
% If program is feasible, W is an SOS.
echo off;