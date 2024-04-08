#include "mx_nlp_ankleTraj.h"
#include "mx_fivePoly.h"

#include<cassert>
#include<iostream>
#include<fstream>
#include<io.h>

#include<vector>

#define MX_INF 2e19

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

using namespace Ipopt;
using namespace std;

extern double xGuess[_NX];

extern double tCtrl;
extern double T_swg;
extern double x0[3], xf[3], z0[3], zf[3];

extern int N_sample;
extern double DT;

extern double x_lb[_NX], x_ub[_NX];

extern ofstream xSoln;

void mx_getAnkleTrajNodes(const Number* z, double dt, double T, vector<double>& foot_x, vector<double>& foot_z)
{
	double S_x1[6], S_x2[6], S_x3[6];
	double S_z1[6], S_z2[6], S_z3[6];

	// x
	mx_getFivePolyCoeff(T / 3.0, x0[0], x0[1], x0[2], z[0], z[1], z[2], S_x1);
	mx_getFivePolyCoeff(T / 3.0,  z[0],  z[1],  z[2], z[3], z[4], z[5], S_x2);
	mx_getFivePolyCoeff(T / 3.0,  z[3],  z[4],  z[5], xf[0], xf[1], xf[2], S_x3);

	// z
	mx_getFivePolyCoeff(T / 3.0, z0[0], z0[1], z0[2], z[6], z[7], z[8], S_z1);
	mx_getFivePolyCoeff(T / 3.0,  z[6],  z[7],  z[8], z[9], z[10], z[11], S_z2);
	mx_getFivePolyCoeff(T / 3.0,  z[9], z[10], z[11], zf[0], zf[1], zf[2], S_z3);

	// x z -- nodes
	int NT = (int) round(T / dt) + 1;
	double tt = 0;
	for (int i = 0; i < NT; i++)
	{
		tt = i * dt;
		if (tt <= T / 3)
		{
			mx_getFivePoly(S_x1, tt, foot_x[i]);
			mx_getFivePoly(S_z1, tt, foot_z[i]);
		}
		else if (tt > T / 3 && tt <= 2 * T / 3)
		{
			mx_getFivePoly(S_x2, tt - T / 3, foot_x[i]);
			mx_getFivePoly(S_z2, tt - T / 3, foot_z[i]);
		}
		else
		{
			mx_getFivePoly(S_x3, tt - 2 * T / 3, foot_x[i]);
			mx_getFivePoly(S_z3, tt - 2 * T / 3, foot_z[i]);
		}
	}

	/*cout << "foot_x = " << endl;
	for (int i = 0; i < NT; i++)
		cout << foot_x[i] << endl;
	cout << "\n\n";
	cout << "foot_z = " << endl;
	for (int i = 0; i < NT; i++)
		cout << foot_z[i] << endl;
	cout << "\n\n";*/
}

void mx_getCstsLite(const Number* z, vector<double>& csts)
{
	vector<double> foot_x(N_sample + 1), foot_z(N_sample + 1);
	mx_getAnkleTrajNodes(z, DT, T_swg, foot_x, foot_z);

	int nn = (int) floor(N_sample / 3);

	int idx_n = 4;

	int idx_temp;

	// c_x1
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 1 + i;
		csts[i] = -foot_x[idx_temp] - 0.05;
		csts[idx_n + i] = foot_x[idx_temp] - 0.00;
	}
	
	// c_x2
	int idx_csts = 2 * idx_n;
	idx_n = 4;	
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		csts[idx_csts + i] = -foot_x[idx_temp] + 0.00;
		csts[idx_csts + idx_n + i] = foot_x[idx_temp] - (xf[0] - 0.01);
	}

	// c_x3
	idx_csts += 2 * idx_n;
	idx_n = 3;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + nn + 1 + i;
		csts[idx_csts + i] = -foot_x[idx_temp] + (xf[0] - 0.01);
		csts[idx_csts + idx_n + i] = foot_x[idx_temp] - xf[0];
	}

	// c_z1
	idx_csts += 2 * idx_n;
	idx_n = 4;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 1 + i;
		csts[idx_csts + i] = -foot_z[idx_temp] + z0[0];
		csts[idx_csts + idx_n + i] = foot_z[idx_temp] - (zf[0] + 0.05);
	}

	// c_z2
	idx_csts += 2 * idx_n;
	idx_n = 4;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		csts[idx_csts + i] = -foot_z[idx_temp] + (zf[0] + 0.01);
		csts[idx_csts + idx_n + i] = foot_z[idx_temp] - (zf[0] + 0.04);
	}

	// c_z3
	idx_csts += 2 * idx_n;
	idx_n = 3;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + nn + 1 + i;
		csts[idx_csts + i] = -foot_z[idx_temp] + zf[0];
		csts[idx_csts + idx_n + i] = foot_z[idx_temp] - (zf[0] + 0.03);
	}

	// c_x4 -- 从第nn+2个开始只增不减
	idx_csts += 2 * idx_n;
	idx_n = 8;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		csts[idx_csts + i] = foot_x[idx_temp - 1] - foot_x[idx_temp];
	}

	// c_z4 -- 从第5个开始只减不增
	idx_csts += idx_n;
	idx_n = 9;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 4 + i;
		csts[idx_csts + i] = foot_z[idx_temp] - foot_z[idx_temp - 1];
	}

	/*cout << "csts = " << endl;
	for (int i = 0; i < csts.size(); i++)
		cout << csts[i] << endl;*/
}

void mx_getAnkleTraNodes_grad(double dt, double T, vector<vector<double> >& foot_x_grad, vector<vector<double> >& foot_z_grad)
{
	int NT = (int)round(T / dt) + 1;
	double tt = 0;

	double Jac1[3], Jac2[6], Jac3[3];

	for (int i = 0; i < NT; i++)
	{
		tt = i * dt;
		if (tt <= T / 3)
		{
			mx_getJacobianRight3(tt, T / 3, Jac1);
			for (int j = 0; j < 3; j++)
			{
				foot_x_grad[i][j] = Jac1[j];
				foot_z_grad[i][j] = Jac1[j];
			}
		}
		else if (tt > T / 3 && tt <= 2 * T / 3)
		{
			mx_getJacobianAll6(tt - T / 3, T / 3, Jac2);
			for (int j = 0; j < 6; j++)
			{
				foot_x_grad[i][j] = Jac2[j];
				foot_z_grad[i][j] = Jac2[j];
			}
		}
		else
		{
			mx_getJacobianLeft3(tt - 2 * T / 3, T / 3, Jac3);
			for (int j = 0; j < 3; j++)
			{
				foot_x_grad[i][j + 3] = Jac3[j];
				foot_z_grad[i][j + 3] = Jac3[j];
			}
		}
	}
}

void mx_getCstsLite_grad(vector<vector<double> >& csts_grad)
{
	vector<vector<double> > foot_x_grad(N_sample + 1, vector<double>(6, 0.0)), foot_z_grad(N_sample + 1, vector<double>(6, 0.0));
	mx_getAnkleTraNodes_grad(DT, T_swg, foot_x_grad, foot_z_grad);

	// 开始！
	int nn = (int)floor(N_sample / 3);
	int idx_n = 4;
	int idx_temp = 0;

	// c_x1_grad
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[i][j] = - foot_x_grad[idx_temp][j];
			csts_grad[i + idx_n][j] = foot_x_grad[idx_temp][j];
		}
	}

	// c_x2_grad
	int idx_csts = 2 * idx_n;
	idx_n = 4;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[idx_csts + i][j] = -foot_x_grad[idx_temp][j];
			csts_grad[idx_csts + idx_n + i][j] = foot_x_grad[idx_temp][j];
		}
	}

	// c_x3_grad
	idx_csts += 2 * idx_n;
	idx_n = 3;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + nn + 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[idx_csts + i][j] = -foot_x_grad[idx_temp][j];
			csts_grad[idx_csts + idx_n + i][j] = foot_x_grad[idx_temp][j];
		}
	}

	// c_z1_grad
	idx_csts += 2 * idx_n;
	idx_n = 4;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[idx_csts + i][6 + j] = -foot_z_grad[idx_temp][j];
			csts_grad[idx_csts + idx_n + i][6 + j] = foot_z_grad[idx_temp][j];
		}
	}

	// c_z2_grad
	idx_csts += 2 * idx_n;
	idx_n = 4;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[idx_csts + i][6 + j] = -foot_z_grad[idx_temp][j];
			csts_grad[idx_csts + idx_n + i][6 + j] = foot_z_grad[idx_temp][j];
		}
	}

	// c_z3_grad
	idx_csts += 2 * idx_n;
	idx_n = 3;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + nn + 1 + i;
		for (int j = 0; j < 6; j++)
		{
			csts_grad[idx_csts + i][6 + j] = -foot_z_grad[idx_temp][j];
			csts_grad[idx_csts + idx_n + i][6 + j] = foot_z_grad[idx_temp][j];
		}
	}

	// c_x4_grad
	idx_csts += 2 * idx_n;
	idx_n = 8;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = nn + 1 + i;
		for(int j = 0; j < 6; j++)
			csts_grad[idx_csts + i][j] = foot_x_grad[idx_temp - 1][j] - foot_x_grad[idx_temp][j];
	}

	// c_z4_grad
	idx_csts += idx_n;
	idx_n = 9;
	for (int i = 0; i < idx_n; i++)
	{
		idx_temp = 4 + i;
		for (int j = 0; j < 6; j++)
			csts_grad[idx_csts + i][6 + j] = foot_z_grad[idx_temp][j] - foot_z_grad[idx_temp - 1][j];
	}

	/*cout << "csts_grad = " << endl;
	for (int i = 0; i < csts_grad.size(); i++)
	{
		for (int j = 0; j < csts_grad[i].size(); j++)
			cout << csts_grad[i][j] << "\t";
		cout << "\n";
	}*/

}

MXNLP::MXNLP()
{ }

MXNLP::~MXNLP()
{ }

bool MXNLP::get_nlp_info(
	Index& n,
	Index& m,
	Index& nnz_jac_g,
	Index& nnz_h_lag,
	IndexStyleEnum& index_style
)
{
	n = _NX;
	m = _NC;
	nnz_jac_g = n * m;
	nnz_h_lag = 0;
	index_style = C_STYLE; // from 0 start

	return true;
}

bool MXNLP::get_bounds_info(
	Index n,
	Number* x_l,
	Number* x_u,
	Index m,
	Number* g_l,
	Number* g_u
)
{
	assert(n == _NX);
	assert(m == _NC);

	for (Index i = 0; i < n; i++)
	{
		x_l[i] = x_lb[i];
		x_u[i] = x_ub[i];
	}

	for (Index i = 0; i < m; i++)
	{
		g_l[i] = -MX_INF;
		g_u[i] = 0;
	}

	return true;
}

bool MXNLP::get_starting_point(
	Index n,
	bool init_x,
	Number* x,
	bool init_z,
	Number* z_L,
	Number* z_U,
	Index m,
	bool init_lambda,
	Number* lambda
)
{
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	for (Index i = 0; i < n; i++)
	{
		x[i] = xGuess[i];
	}

	return true;
}

bool MXNLP::eval_f(
	Index n,
	const Number* x,
	bool new_x,
	Number& obj_value
)
{
	obj_value = 0.0;

	return true;
}

bool MXNLP::eval_grad_f(
	Index n,
	const Number* x,
	bool new_x,
	Number* grad_f
)
{
	assert(n == _NX);

	for (Index i = 0; i < n; i++)
	{
		grad_f[i] = 0.0;
	}

	return true;
}

bool MXNLP::eval_g(
	Index n,
	const Number* x,
	bool new_x,
	Index m,
	Number* g
)
{
	assert(n == _NX);
	assert(m == _NC);

	vector<double> csts(m, 0.0);

	mx_getCstsLite(x, csts);

	for (int i = 0; i < m; i++)
		g[i] = csts[i];

	return true;
}

bool MXNLP::eval_jac_g(
	Index n,
	const Number* x,
	bool new_x,
	Index m,
	Index nele_jac,
	Index* iRow,
	Index* jCol,
	Number* values
)
{
	assert(n == _NX);
	assert(m == _NC);

	if (values == NULL)
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				iRow[i * n + j] = i;
				jCol[i * n + j] = j;
			}
		}
	}
	else
	{
		vector<vector<double> > csts_grad(m, vector<double>(n, 0.0));
		mx_getCstsLite_grad(csts_grad);

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				values[i * n + j] = csts_grad[i][j];
	}

	return true;
}

bool MXNLP::eval_h(
	Index n,
	const Number* x,
	bool new_x,
	Number obj_factor,
	Index m,
	const Number* lambda,
	bool new_lambda,
	Index nele_hess,
	Index* iRow,
	Index* jCol,
	Number* values
)
{
	return true;
}

void MXNLP::finalize_solution(
	SolverReturn status,
	Index n,
	const Number* x,
	const Number* z_L,
	const Number* z_U,
	Index m,
	const Number* g,
	const Number* lambda,
	Number obj_value,
	const IpoptData* ip_data,
	IpoptCalculatedQuantities* ip_cq
)
{
	cout << endl << endl << "Solution of the primal variables, x" << endl;
	for (Index i = 0; i < n; i++)
	{
		cout << "x[" << i << "] = " << x[i] << endl;
	}

	/*if ((_access("xSoln_cpp.dat", 0)) == -1)
	{
		xSoln.open("xSoln_cpp.dat");
		for (int i = 0; i < n; i++)
		{
			xSoln << x[i] << endl;
		}
		xSoln.close();
	}*/

	// write the solutions into the file
	xSoln.open("xSoln_cpp.dat");
	for (int i = 0; i < n; i++)
	{
		xSoln << x[i] << endl;
	}
	xSoln.close();

	/*cout << endl << endl << "Solution of the bound multipliers, z_L and z_U" << endl;
	for (Index i = 0; i < n; i++)
	{
		cout << "z_L[" << i << "] = " << z_L[i] << endl;
	}
	for (Index i = 0; i < n; i++)
	{
		cout << "z_U[" << i << "] = " << z_U[i] << endl;
	}

	cout << endl << endl << "Objective value" << endl;
	cout << "f(x*) = " << obj_value << endl;

	cout << endl << endl << "Final value of the contraints:" << endl;
	for (Index i = 0; i < m; i++)
	{
		cout << "g(" << i << ") = " << g[i] << endl;
	}*/
}