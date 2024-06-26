function Jac2 = autoGen_jacobianAll6(tx,T)
%AUTOGEN_JACOBIANALL6
%    JAC2 = AUTOGEN_JACOBIANALL6(TX,T)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    28-Jun-2023 20:12:00

t2 = T.^2;
t3 = tx.^2;
t4 = tx.^3;
t5 = 1.0./T.^3;
t7 = 1.0./T.^5;
t8 = -tx;
t6 = 1.0./t2.^2;
t9 = t3.*6.0;
t10 = T+t8;
t11 = t10.^3;
Jac2 = [t7.*t11.*(t2+t9+T.*tx.*3.0),t6.*t11.*tx.*(T+tx.*3.0),(t3.*t5.*t11)./2.0,t4.*t7.*(t2.*1.0e+1+t9-T.*tx.*1.5e+1),-t4.*t6.*(t2.*4.0+t3.*3.0-T.*tx.*7.0),(t4.*t5.*t10.^2)./2.0];
