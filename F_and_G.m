function [f,g] = F_and_G(aprx_strt_t,ex_times,k,coefs)

% Here are the coefficients
V_C_coefs = coefs{1,1};
V_u_coefs = coefs{1,2};
P_C_coefs = coefs{2,1};
P_u_coefs = coefs{2,2};
A_C_coefs = coefs{3,1};
A_u_coefs = coefs{3,2};

% Create splines for the coefficients, then put in polynomial form
template = spapi(k,ex_times,ones(1,length(ex_times)));

VC = template;   VC.coefs  = V_C_coefs;       VC  = fn2fm(VC,'pp');
Vu1 = template;  Vu1.coefs = V_u_coefs(1,:);  Vu1 = fn2fm(Vu1,'pp');
Vu2 = template;  Vu2.coefs = V_u_coefs(2,:);  Vu2 = fn2fm(Vu2,'pp');

PC = template;   PC.coefs  = P_C_coefs;       PC = fn2fm(PC,'pp');
Pu1 = template;  Pu1.coefs = P_u_coefs(1,:);  Pu1 = fn2fm(Pu1,'pp');
Pu2 = template;  Pu2.coefs = P_u_coefs(2,:);  Pu2 = fn2fm(Pu2,'pp');

AC = template;   AC.coefs  = A_C_coefs;       AC = fn2fm(AC,'pp');
Au1 = template;  Au1.coefs = A_u_coefs(1,:);  Au1 = fn2fm(Au1,'pp');
Au2 = template;  Au2.coefs = A_u_coefs(2,:);  Au2 = fn2fm(Au2,'pp');

% Turn the extrapolated splines into symbolic functions, and take deriv.
syms t
t_mons = monomials(t-aprx_strt_t,0:k-1);

f_V(t) = diff(flip(VC.coefs(end,:))*t_mons);
g_V1(t) = diff(flip(Vu1.coefs(end,:))*t_mons);
g_V2(t) = diff(flip(Vu2.coefs(end,:))*t_mons);

f_P(t) = diff(flip(PC.coefs(end,:))*t_mons);
g_P1(t) = diff(flip(Pu1.coefs(end,:))*t_mons);
g_P2(t) = diff(flip(Pu2.coefs(end,:))*t_mons);

f_A(t) = diff(flip(AC.coefs(end,:))*t_mons);
g_A1(t) = diff(flip(Au1.coefs(end,:))*t_mons);
g_A2(t) = diff(flip(Au2.coefs(end,:))*t_mons);

% Return f and g
f(t) = [f_V(t); f_P(t); f_A(t)];
g(t) = [g_V1(t),g_V2(t); g_P1(t),g_P2(t); g_A1(t),g_A2(t)];