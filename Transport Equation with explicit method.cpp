//Group Members: Yuzheng Nan (ynan4), Mengyuan Chen (mchen100)
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

void transport_equation(double C, double dx, double dt, ofstream& outfile)
{
	cout << "When dx=" << dx << ":" << endl;
	outfile << "x,";
	int M = int(5 / dx);
	int steps = 1 / dt;
	double lamda = C*dt / dx;
	vector<double>u(M + 1);
	//boundary condition for t
	for (int i = 0; i <= M; i++)
	{
		if ((-2 + i*dx) < -1)
			u[i] = 0;
		if (((-2 + i*dx) >= -1) && ((-2 + i*dx) <= 0))
			u[i] = -1 + i*dx;
		if ((-2 + i*dx) > 0)
			u[i] = 1;
	}
	vector<double>u_next(M + 1);
	for (int i = 0; i <= M; i++)
	{
		outfile << -2 + i*dx << ",";
	}
	outfile << endl;
	for (int j = 1; j <= steps; j++)
	{
		outfile << dt*(j - 1) << ",";
		cout << "t=" << dt*(j - 1) << ": ";
		u_next[0] = 0; u_next[M] = 1;//boundary conditions for x
		for (int i = 1; i < M; i++)
			u_next[i] = u[i] + lamda*(u[i - 1] - u[i]);//implementing explicit scheme
		for (int i = 0; i <= M; i++)
		{
			cout << u[i] << ",";
			outfile<< u[i] << ",";
			u[i] = u_next[i];
		}
		cout << endl;
		outfile << endl;
	}
	outfile << 1 << ",";
	cout << "t=" << 1 << ": ";
	for (int i = 0; i <= M; i++)
	{
		cout << u[i] << ",";
		outfile << u[i] << ",";
	}
	cout << endl;
	outfile << endl;
}

int main()
{
	double C = 1; double dx1 = 0.05; double dx2 = 0.01; double dx3 = 0.005; double dt = 0.01;
	ofstream file1("file1.csv");
	transport_equation(C, dx1, dt, file1);
	ofstream file2("file2.csv");
	transport_equation(C, dx2, dt, file2);
	ofstream file3("file3.csv");
	transport_equation(C, dx3, dt, file3);

}