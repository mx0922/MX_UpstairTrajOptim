#include"mx_fivePoly.h"
#include<cmath>

void mx_getFivePolyCoeff(double T, double p0, double v0, double a0, double p1, double v1, double a1, double S[6])
{
	S[0] = p0;
	S[1] = v0;
	S[2] = a0 / 2.0;

	double t2 = T * T;
	double t3 = a1 * t2;
	double t4 = a0 * t2 * 3;
	double t5 = -t3;
	S[3] = 1.0 / pow(T, 3) * (p0 * 20 - p1 * 20 + t4 + t5 + T * v0 * 12 + T * v1 * 8) * (-1.0 / 2.0);

	S[4] = 1.0 / (t2 * t2) * (p0 * 30 - p1 * 30 - t3 * 2.0 + t4 + T * v0 * 16 + T * v1 * 14) / 2.0;

	S[5] = 1.0 / pow(T, 5) * (p0 * 12 - p1 * 12 + t5 + T * v0 * 6.0 + T * v1 * 6.0 + a0 * t2) * (-1.0 / 2.0);
}

void mx_getFivePoly(double S[6], double t, double &pos)
{
	pos = S[0] + S[1] * t + S[2] * t * t + S[3] * pow(t, 3) + S[4] * pow(t, 4) + S[5] * pow(t, 5);
}

void mx_getJacobianRight3(double tx, double T, double Jac[])
{
	double t2 = T * T;
	double t3 = tx * tx;
	double t4 = t3 * tx;

	Jac[0] = 1.0 / pow(T, 5) * t4 * (t2 * 10.0 + t3 * 6.0 - T * tx * 15.0);
	Jac[1] = -1.0 / (t2 * t2) * t4 * (t2 * 4.0 + t3 * 3.0 - T * tx * 7.0);
	Jac[2] = (1.0 / (t2 * T) * t4 * (T - tx) * (T - tx)) / 2.0;
}

void mx_getJacobianAll6(double tx, double T, double Jac[])
{
	double t2 = T * T;
	double t3 = tx * tx;
	double t4 = t3 * tx;
	double t5 = 1.0 / pow(T, 3);
	double t7 = 1.0 / pow(T, 5);
	double t8 = -tx;

	double t6 = 1.0 / (t2 * t2);
	double t9 = t3 * 6.0;
	double t10 = T + t8;
	double t11 = pow(t10, 3);

	Jac[0] = t7 * t11 * (t2 + t9 + T * tx * 3.0);
	Jac[1] = t6 * t11 * tx * (T + tx * 3.0);
	Jac[2] = (t3 * t5 * t11) / 2.0;
	Jac[3] = t4 * t7 * (t2 * 10.0 + t9 - T * tx * 15.0);
	Jac[4] = -t4 * t6 * (t2 * 4.0 + t3 * 3.0 - T * tx * 7.0);
	Jac[5] = (t4 * t5 * t10 * t10) / 2.0;
}

void mx_getJacobianLeft3(double tx, double T, double Jac[])
{
	double t2 = tx * tx;
	double t3 = -tx;
	double t4 = T + t3;
	double t5 = pow(t4, 3);

	Jac[0] = 1.0 / pow(T, 5) * t5 * (t2 * 6.0 + T * tx * 3.0 + T * T);
	Jac[1] = 1.0 / pow(T, 4) * t5 * tx * (T + tx * 3.0);
	Jac[2] = (1.0 / pow(T, 3) * t2 * t5) / 2.0;
}