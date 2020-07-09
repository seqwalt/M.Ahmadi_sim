% Note: S1_no_u = f(t), [S1_u1 S2_u2] = g(t)
% Form: dXdt = F(t) + G(t)*u
% F = [f_{vel} f_{FPA} f_{alt}]' ,  G = [g_{vel} g_{FPA} g_{alt}]
% T -> Thrust
% A -> Angle of attack

function dxdt = approx_sys(t,system,T,A,si1,si2)

V1_no_u = system(1,1);
V2_no_u = system(1,2);
V1_u1   = system(1,3);
V2_u2   = system(1,4);
P1_no_u = system(2,1);
P2_no_u = system(2,2);
P1_u1   = system(2,3);
P2_u2   = system(2,4);
A1_no_u = system(3,1);
A2_no_u = system(3,2);
A1_u1   = system(3,3);
A2_u2   = system(3,4);

f_V1 = fnval(fnder(V1_no_u),t);
f_V2 = fnval(fnder(V2_no_u),t);   f_V = f_V1 + f_V2;
g_V1 = fnval(fnder(V1_u1),t);
g_V2 = fnval(fnder(V2_u2),t);
f_P1 = fnval(fnder(P1_no_u),t);
f_P2 = fnval(fnder(P2_no_u),t);   f_P = f_P1 + f_P2;
g_P1 = fnval(fnder(P1_u1),t);
g_P2 = fnval(fnder(P2_u2),t);
f_A1 = fnval(fnder(A1_no_u),t);
f_A2 = fnval(fnder(A2_no_u),t);   f_A = f_A1 + f_A2;
g_A1 = fnval(fnder(A1_u1),t);
g_A2 = fnval(fnder(A2_u2),t);

dxdt = zeros(3,1);
dxdt(1) = f_V + [g_V1,g_V2]*[T;A] + si1;
dxdt(2) = f_P + [g_P1,g_P2]*[T;A] + si2;
dxdt(3) = f_A + [g_A1,g_A2]*[T;A];
