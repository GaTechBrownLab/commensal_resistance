%% Benefits of antibiotic resistance in commensals and the scope for resistance optimization
% Kristofer Wollein Waldetoft, Sarah Sundius, Rachel Kuske, Sam P. Brown

% Code by: Sarah Sundius
% August 18, 2021


%% Qualitative Analysis
% code to complete qualitative (stability) analysis in SI: defines model,
% equilibrium points, and stability conditions, assuming presence of
% antibiotic (A > 0)

% variables:
% P = pathogen species absolute density
% C = commensal species absolute density

% parameters: 
% rp,rc = maximal growth rate
% kp,kc = carrying capacity
% apc, acp = interspecific interaction coeff
% x = maximal pathogen clearance rate
% A = antibiotic exposure (control parameter)
% f = commensal relative susceptibility (control parameter)

clear; clc;

% define model
syms P C rp rc kp kc apc acp f A x

dPdt = rp*P*(1 - (P+apc*C)/kp) - x*A*P;
dCdt = rc*C*(1 - (C+acp*P)/kc) - x*f*A*C;

% solve for equilibria
eqn_P = dPdt == 0;
eqn_C = dCdt == 0;

eqns = [eqn_P, eqn_C];
vars = [P, C];

[P_eq, C_eq] = vpasolve(eqns, vars);

eq1 = [P_eq(1), C_eq(1)];
eq2 = [P_eq(2), C_eq(2)];
eq3 = [P_eq(3), C_eq(3)];
eq4 = [P_eq(4), C_eq(4)];

% check equilibria - sub in each eq1-4, both should be 0
%dPdtcheck = simplify(subs(dPdt,[P,C],eq4))
%dCdtcheck = simplify(subs(dCdt,[P,C],eq4))

% compute the Jacobian
model = [dPdt; dCdt];
J = jacobian(model,vars);

% linearize system (compute Jacobian at equilibria)
J1 = subs(J,vars,eq1);
J2 = subs(J,vars,eq2);
J3 = subs(J,vars,eq3);
J4 = subs(J,vars,eq4);

% calculate eigenvalues/eigenvectors
[V1,D1] = eig(J1);
[V2,D2] = eig(J2);
[V3,D3] = eig(J3);
[V4,D4] = eig(J4);

% define/calculate linear stability conditions: tr(J) < 0, det(J) > 0 
cond1_1 = simplify(trace(J1));               
cond2_1 = simplify(det(J1));  

cond1_2 = simplify(trace(J2));                 
cond2_2 = simplify(det(J2));

cond1_3 = simplify(trace(J3));                    
cond2_3 = simplify(det(J3)); 

cond1_4 = simplify(trace(J4));                   
cond2_4 = simplify(det(J4));   



