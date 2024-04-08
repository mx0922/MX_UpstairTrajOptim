#include"IpIpoptApplication.hpp"
#include"IpSolveStatistics.hpp"

#include"mx_nlp_ankleTraj.h"
#include"mx_fivePoly.h"

#include<iostream>
#include<fstream>

using namespace std;
using namespace Ipopt;

// some parameters
double tCtrl = 0.004;  // control period: 4ms
double T_swg = 0.8;    // swing time

// 上高22cm、长30cm的台阶 -- 粗略的相对位置
// 真实的相对位置可以通过修改xf，zf测试

// boundary footholds (in sagittal plane)
double x0[3] = { 0.00, 0.0, 0.0 }; // pos vel acc (vel, acc -- default = 0.0)
double xf[3] = { 0.30, 0.0, 0.0 };

double z0[3] = { 0.08, 0.0, 0.0 };
double zf[3] = { 0.30, 0.0, 0.0 };

// number of sampling point (user-defined: do not suggest too many)
int N_sample = 12;
double DT = T_swg / N_sample; // sampling timing

// x guess
double xGuess[_NX] = { -0.05, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0, -0.2, 0.0};

// x_lb and x_ub -- necessary setting!
double x_lb[_NX] = { -0.20, -1.0, -2.0, 0.0, -1.0, -2.0, 0.0, -1.0, -10, 0.0, -1.0, -10 };
double x_ub[_NX] = {  0.05,  1.0,  2.0, 0.5,  1.0,  2.0, 0.5,  1.0,  10, 0.5,  1.0,  10 };

// save x solution file to show the results in matlab
ofstream xSoln;

// initia guess by the boundary footholds position
void initGuess()
{
	xGuess[3] = xf[0];
	xGuess[6] = zf[0];
	xGuess[9] = zf[0];
}

// functions
int main(
	int,
	char**
)
{
	// xGuess initialization
	initGuess();

	SmartPtr<TNLP> mxnlp1 = new MXNLP();

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	// // ===== only work in Release mode!!! =====
	app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetIntegerValue("print_level", 3); // the value smaller, the print output less

	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded)
	{
		cout << endl << endl << "*** Error during initialization!" << endl;
		return (int)status;
	}

	// solve the NLP using IPOPT solver
	status = app->OptimizeTNLP(mxnlp1);

	if(status == Solve_Succeeded)
		cout << "\n\n" << "*** The problem solved!" << endl;
	else
		cout << "\n\n" << "!!! The problem FAILED!" << endl;

	return (int)status;
}

int main_test33()
{

	double S[6] = { 0 };
	mx_getFivePolyCoeff(1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, S);

	for (int i = 0; i < 6; i++)
		cout << "S[" << i << "] = " << S[i] << endl;

	cout << "\n\n";

	double pos;
	mx_getFivePoly(S, 0.5, pos);

	cout << "pos = " << pos << endl;

	cout << "\n\n";

	double tx = 0.2;
	double T = 0.6;

	double Jac1[3], Jac2[6], Jac3[3];

	mx_getJacobianRight3(tx, T, Jac1);
	mx_getJacobianAll6(tx, T, Jac2);
	mx_getJacobianLeft3(tx, T, Jac3);

	for(int i = 0; i < 3; i++)
		cout << "Jac1[" << i << "] = " << Jac1[i] << endl;

	cout << "\n\n";

	for (int i = 0; i < 6; i++)
		cout << "Jac2[" << i << "] = " << Jac2[i] << endl;

	cout << "\n\n";

	for (int i = 0; i < 3; i++)
		cout << "Jac3[" << i << "] = " << Jac3[i] << endl;

	cout << "\n\n";

	return 0;
}
