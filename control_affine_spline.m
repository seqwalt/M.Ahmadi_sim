% Finding coeffs c_{i,j} is easier if we say the final spline is a linear
% combination of the deg 1 spline and deg 2 spline, where
% Final_Spline = deg_1_sp + deg_2_sp
% r1*Final_Spline = deg_1_sp ,  r2*Final_Spline = deg_2_sp
% r1 + r2 = 1

% Find a_{i,j} by setting one of the inputs to 0, and the other non
% zero, then put sum into the following form:
% sum(c_{i,1}*a_{i,1}*Q_{i,1}(t)), which equals dS1 = (S_u1 - S_no_u)/u1.
% This allows for solving for a_{i,1} after finding the coeffs
% determined by spapi.  The coeffs will equal c_{i,1}*a_{i,1}, and we
% already know c_{i,1}.

% This function outputs the extrapolated versions of the splines, up until
% ex_dur after the data stops.

function [S1_no_u,S2_no_u,S1_u1,S2_u2] = control_affine_spline(output_type,times,ex_dur,step,S_true,S_u1,S_u2,S_no_u,u1,u2)

if strcmpi(output_type(1),'v') == 1
    column = 2;
elseif strcmpi(output_type(1),'f') == 1
    column = 3;
elseif strcmpi(output_type(1),'a') == 1
    column = 4;
end

dS1 = (S_u1 - S_no_u)/u1;
dS2 = (S_u2 - S_no_u)/u2;
r1 = 0.5;
r2 = 1 - r1;
k1 = 2;   % degree: 1 = k1 - 1
k2 = 3;   % degree: 2 = k2 - 1

ext_times = times(1):step:times(end) + ex_dur;

%% --- Extrapolate for coefficients when u_{j} = 0 --- %%
S1_no_u = fnxtr(spapi(k1,times,r1*S_no_u(:,column)),3);
ex1_no_u = fnval(S1_no_u,ext_times);
S1_no_u = spapi(k1,ext_times,ex1_no_u);

S2_no_u = fnxtr(spapi(k2,times,r2*S_no_u(:,column)),3);
ex2_no_u = fnval(S2_no_u,ext_times);
S2_no_u = spapi(k2,ext_times,ex2_no_u);

c1 = S1_no_u.coefs;
c2 = S2_no_u.coefs;

%% --- Extrapolate for coefficients when u_{2} = 0 --- %%
S1_u1 = fnxtr(spapi(k1,times,dS1(:,column)),3);
ex1_u1 = fnval(S1_u1,ext_times);
S1_u1 = spapi(k1,ext_times,ex1_u1);

ca1 = S1_u1.coefs;
a1 = ca1./c1;

%% --- Extrapolate for coefficients when u_{1} = 0 --- %%
S2_u2 = fnxtr(spapi(k2,times,dS2(:,column)),3);
ex2_u2 = fnval(S2_u2,ext_times);
S2_u2 = spapi(k2,ext_times,ex2_u2);

ca2 = S2_u2.coefs;
a2 = ca2./c2;

%% --- Extrapolate for coefficients c_{i,j} with actual flight inputs --- %%
S1_true = fnxtr(spapi(k1,times,r1*S_true(:,column)),3);
ex1_true = fnval(S1_true,ext_times);
S1_true = spapi(k1,ext_times,ex1_true);

S2_true = fnxtr(spapi(k2,times,r2*S_true(:,column)),3);
ex2_true = fnval(S2_true,ext_times);
S2_true = spapi(k2,ext_times,ex2_true);

% h_{i,j} = c_{i,j}*(1 + a_{i,j}*u_{j})
h1 = S1_true.coefs;
h2 = S2_true.coefs;
c1 = h1./(1 + a1*u1);
c2 = h2./(1 + a2*u2);

%% --- Correct the splines with the new c_{i,j} values --- %%
S1_no_u.coefs = c1;
S2_no_u.coefs = c2;
S1_u1.coefs = c1.*a1;
S2_u2.coefs = c2.*a2;