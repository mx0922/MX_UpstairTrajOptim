function [pos, vel, acc] = getFiveOrderPoly(S, t)

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;
f = S.f;

pos = a + b .* t +  c .* t.^2 +     d .* t.^3 +      e .* t.^4 +      f .* t.^5;
vel =          b + 2 * c .* t + 3 * d .* t.^2 +  4 * e .* t.^3 +  5 * f .* t.^4;
acc =                  2 * c  +    6 * d .* t + 12 * e .* t.^2 + 20 * f .* t.^3;

end