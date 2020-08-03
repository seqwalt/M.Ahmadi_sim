% This function outputs the extrapolated coefficients of the splines,
% that is, C_i, and the input weight coefs -- see the the Ahmadi_sim.pdf
% write-up.  Note: u_coefs = [u1_coefs; u2_coefs]

function [C_coefs, u_coefs] = spline_coefs(output_type,times,ex_times,order,ex_order,S_no_u,S_a,S_b,u1a,u2a,u1b,u2b)

if strcmpi(output_type(1),'vel') == 1
    column = 2;
elseif strcmpi(output_type(1),'fpa') == 1
    column = 3;
elseif strcmpi(output_type(1),'alt') == 1
    column = 4;
end

%% --- Extrapolate for coefficients when u_{j} = 0 (C)--- %%
exS_no_u = fnxtr(spapi(order,times,S_no_u(:,column)),ex_order); % exS_no_u is in ppform
exS_no_u_vals = fnval(exS_no_u,ex_times);
exS_no_u = spapi(order,ex_times,exS_no_u_vals);                 % now its in B-form (what we want)

C = exS_no_u.coefs;

%% --- Extrapolate for 2nd trajectory coefficients (A) --- %%
exS_a = fnxtr(spapi(order,times,S_a(:,column)),ex_order); % exS_a is in ppform
exS_a_vals = fnval(exS_a,ex_times);
exS_a = spapi(order,ex_times,exS_a_vals);                 % now in B-form

A = exS_a.coefs;

%% --- Extrapolate for 3rd (final) trajectory coefficients (B) --- %%
exS_b = fnxtr(spapi(order,times,S_b(:,column)),ex_order); % exS_b is in ppform
exS_b_vals = fnval(exS_b,ex_times);
exS_b = spapi(order,ex_times,exS_b_vals);                 % now in B-form

B = exS_b.coefs;

%% --- Calculate Spline Coefficients --- %%
u = [u1a, u2a; u1b, u2b];
C_coefs = C;
u_coefs = (1/det(u))*[u2b*(A-C) + u2a*(C-B); u1b*(C-A) + u1a*(B-C)];