% Note: S1_no_u = f(t), [S1_u1 S2_u2] = g(t)
% Form: dXdt = F(t) + G(t)*u
% F = [f_{vel} f_{FPA} f_{alt}]' ,  G = [g_{vel} g_{FPA} g_{alt}]
% T -> Thrust
% A -> Angle of attack

function dxdt = approx_sys(t,ex_times,k,coefs,T,A,si1,si2)
% Here are the coefficients
V_C_coefs = coefs{1,1};
V_u_coefs = coefs{1,2};
P_C_coefs = coefs{2,1};
P_u_coefs = coefs{2,2};
A_C_coefs = coefs{3,1};
A_u_coefs = coefs{3,2};

% Create splines for the coefficients
template = spapi(k,ex_times,ones(1,length(ex_times)));

VC = template; VC.coefs = V_C_coefs;
Vu1 = template; Vu1.coefs = V_u_coefs(1,:);
Vu2 = template; Vu2.coefs = V_u_coefs(2,:);

PC = template; PC.coefs = P_C_coefs;
Pu1 = template; Pu1.coefs = P_u_coefs(1,:);
Pu2 = template; Pu2.coefs = P_u_coefs(2,:);

AC = template; AC.coefs = A_C_coefs;
Au1 = template; Au1.coefs = A_u_coefs(1,:);
Au2 = template; Au2.coefs = A_u_coefs(2,:);

% Take derivatives of the splines and get values at time t
f_V = fnval(fnder(VC),t);
g_V1 = fnval(fnder(Vu1),t);
g_V2 = fnval(fnder(Vu2),t);

f_P = fnval(fnder(PC),t);
g_P1 = fnval(fnder(Pu1),t);
g_P2 = fnval(fnder(Pu2),t);

f_A = fnval(fnder(AC),t);
g_A1 = fnval(fnder(Au1),t);
g_A2 = fnval(fnder(Au2),t);

% Put the model together
dxdt = zeros(3,1);
dxdt(1) = f_V + [g_V1,g_V2]*[T;A] + si1;
dxdt(2) = f_P + [g_P1,g_P2]*[T;A] + si2;
dxdt(3) = f_A + [g_A1,g_A2]*[T;A];
