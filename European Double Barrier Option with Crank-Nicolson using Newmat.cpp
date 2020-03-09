//Group Members: Yuzheng Nan (ynan4), Mengyuan Chen (mchen100)
#include <cmath>
#include "newmat.h" // definitions for newmat matrix library
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

double EuDOputCK(double S, double K, double r, double T, double sigma, double Sb, double Smax, double dS, double dt)
{
	int M = (Smax - Sb) / dS;
	int N = T / dt;

	vector<double>vetS(M + 1);
	for (int m = 0; m <= M; m++)
		vetS[m] = Sb + m*dS;
	ColumnVector f(M - 1);
	ColumnVector f_next(M - 1);
	for (int i = 0; i < M - 1; i++)
		f.element(i) = max(K - vetS[i + 1], 0.0);//boundary condition on time
	vector<double>alpha(M + 1);
	vector<double>beta(M + 1);
	vector<double>gamma(M + 1);
	for (int m = 0; m <= M; m++)
	{
		int i = m + Sb / dS;
		alpha[m] = 0.25*dt*(pow(sigma, 2)*pow(i, 2) - r*i);
		beta[m] = -0.5*dt*(pow(sigma, 2)*pow(i, 2) + r);
		gamma[m] = 0.25*dt*(pow(sigma, 2)*pow(i, 2) + r*i);
	}
	BandMatrix M1(M - 1, 1, 1); M1 = 0.0;
	M1.element(0, 0) = 1 - beta[1]; M1.element(0, 1) = -gamma[1];
	M1.element(M - 2, M - 2) = 1 - beta[M - 1]; M1.element(M - 2, M - 3) = -alpha[M - 1];
	BandMatrix M2(M - 1, 1, 1); M2 = 0.0;
	M2.element(0, 0) = 1 + beta[1]; M2.element(0, 1) = gamma[1];
	M2.element(M - 2, M - 2) = 1 + beta[M - 1]; M2.element(M - 2, M - 3) = alpha[M - 1];
	for (int j = 1; j<M-2; ++j) 
	{
		M1.element(j, j - 1) = -alpha[j + 1];
		M1.element(j, j) = 1 - beta[j + 1];
		M1.element(j, j + 1) = -gamma[j + 1];
		M2.element(j, j - 1) = alpha[j + 1];
		M2.element(j, j) = 1 + beta[j + 1];
		M2.element(j, j + 1) = gamma[j + 1];
	}
	for (int n = N; n >= 1; n--)
	{
		f_next = M1.i()*M2*f;
		f = f_next;
	}
	vector<double>fval(M + 1);
	fval[0] = 0; fval[M] = 0;//These are the boundary conditions on stock price
	for (int m = 1; m < M; m++)
		fval[m] = f.element(m - 1);
	for (int m = 0; m <= M; m++)
		if (vetS[m] == S)
			return fval[m];
}

int main()
{
	double S = 50; double K = 50; double r = 0.1; double T = 5 / double(12); double sigma = 0.4; double Sb = 40; double Smax = 100;
	double dS = 0.5; double dt = 1 / double(1200);
	cout << "The Price of the European Down and Out Barrier Option using Crank-Nicolson Scheme is:" << endl;
	cout << EuDOputCK(S, K, r, T, sigma, Sb, Smax, dS, dt) << endl;
}