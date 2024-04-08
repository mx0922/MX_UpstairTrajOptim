#pragma once
#ifndef _MX_FIVEPOLY_H
#define _MX_FIVEPOLY_H

void mx_getFivePolyCoeff(double T, double p0, double v0, double a0, double p1, double v1, double a1, double S[6]);
void mx_getFivePoly(double S[6], double t, double& pos);

void mx_getJacobianRight3(double tx, double T, double Jac[]);
void mx_getJacobianAll6(double tx, double T, double Jac[]);
void mx_getJacobianLeft3(double tx, double T, double Jac[]);

#endif // !_MX_FIVEPOLY_H
