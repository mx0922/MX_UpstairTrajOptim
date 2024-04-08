function [a,b,c,d,e,f] = autoGen_fiveOrderPolyCoeff(T,p0,v0,a0,p1,v1,a1)
%AUTOGEN_FIVEORDERPOLYCOEFF
%    [A,B,C,D,E,F] = AUTOGEN_FIVEORDERPOLYCOEFF(T,P0,V0,A0,P1,V1,A1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    18-Apr-2022 15:00:49

a = p0;
if nargout > 1
    b = v0;
end
if nargout > 2
    c = a0./2.0;
end
if nargout > 3
    t2 = T.^2;
    t3 = a1.*t2;
    t4 = a0.*t2.*3.0;
    t5 = -t3;
    d = 1.0./T.^3.*(p0.*2.0e+1-p1.*2.0e+1+t4+t5+T.*v0.*1.2e+1+T.*v1.*8.0).*(-1.0./2.0);
end
if nargout > 4
    e = (1.0./t2.^2.*(p0.*3.0e+1-p1.*3.0e+1-t3.*2.0+t4+T.*v0.*1.6e+1+T.*v1.*1.4e+1))./2.0;
end
if nargout > 5
    f = 1.0./T.^5.*(p0.*1.2e+1-p1.*1.2e+1+t5+T.*v0.*6.0+T.*v1.*6.0+a0.*t2).*(-1.0./2.0);
end
