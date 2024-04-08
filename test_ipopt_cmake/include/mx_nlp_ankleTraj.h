#pragma once
#ifndef __MXNLP_HPP__
#define __MXNLP_HPP__

#include "IpTNLP.hpp"

#include<vector>

#define _NX 12 // number of decision variables
#define _NC 61 // number of nonlinear csts

using namespace Ipopt;
using namespace std;

void mx_getAnkleTrajNodes(const Number* z, double dt, double T, vector<double>& foot_x, vector<double>& foot_z);
void mx_getCstsLite(const Number* z, vector<double>& csts);
void mx_getAnkleTraNodes_grad(double dt, double T, vector<vector<double> >& foot_x_grad, vector<vector<double> >& foot_z_grad);
void mx_getCstsLite_grad(vector<vector<double> >& csts_grad);

class MXNLP : public TNLP
{
public:
	// 构造函数
	MXNLP();
	// 析构函数
	virtual ~MXNLP();

	// get nlp info
	virtual bool get_nlp_info(
		Index& n,
		Index& m,
		Index& nnz_jac_g,
		Index& nnz_h_lag,
		IndexStyleEnum& index_style
	);

	// get bounds info
	virtual bool get_bounds_info(
		Index n,
		Number* x_l,
		Number* x_u,
		Index m,
		Number* g_l,
		Number* g_u
	);

	// get starting point
	virtual bool get_starting_point(
		Index n,
		bool init_x,
		Number* x,
		bool init_z,
		Number* z_L,
		Number* z_U,
		Index m,
		bool init_lamnda,
		Number* lambda
	);

	// eval f
	virtual bool eval_f(
		Index n,
		const Number* x,
		bool new_x,
		Number& obj_value
	);

	// eval grad f
	virtual bool eval_grad_f(
		Index n,
		const Number* x,
		bool new_x,
		Number* grad_f
	);

	// eval g
	virtual bool eval_g(
		Index n,
		const Number* x,
		bool new_x,
		Index m,
		Number* g
	);

	// eval_jac_g
	virtual bool eval_jac_g(
		Index n,
		const Number* x,
		bool new_x,
		Index m,
		Index nele_jac,
		Index* iRow,
		Index* jCol,
		Number* values
	);

	// eval_h
	virtual bool eval_h(
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
	);

	virtual void finalize_solution(
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
	);

private:
	//@{
	MXNLP(
		const MXNLP&
	);

	MXNLP& operator=(
		const MXNLP&
		);
	//@}
};

#endif