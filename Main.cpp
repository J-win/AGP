#include <iostream>
using namespace std;

double f(double x)
{
	return x * x;
}

void insertup(double *x, int k, double xk)
{
	bool l = true;
	int i;
	for (i = k - 1; (i >= 0) && (l); i--)
	{
		if (xk < x[i])
			x[i + 1] = x[i];
		else
			l = false;
	}
	if (l)
		x[0] = xk;
	else
		x[i + 2] = xk;
}

int main()
{
	double a, b, ee, m, r, mm, rr, maxrr, minf, minx;
	int k, n, t;
	r = 2.0;
	cout << "Write left and right side" << endl;
	cin >> a >> b;
	cout << "Write number iterations" << endl;
	cin >> n;
	cout << "Write inaccuracy" << endl;
	cin >> ee;
	double* x = new double[n + 3];
	x[0] = a;
	x[1] = b;
	k = 2;
	do
	{
		mm = 0.0;
		maxrr = -1000000.0;
		for (int i = 1; i < k; i++)
		{
			double max = fabs((f(x[i]) - f(x[i - 1])) / (x[i] - x[i - 1]));
			if (mm < max)
				mm = max;
		}
		if (mm > 0.0)
			m = r * mm;
		else
			m = 1.0;
		for (int i = 1; i < k; i++)
		{
			double rr = m * (x[i] - x[i - 1]) + (((f(x[i]) - f(x[i - 1])) * (f(x[i]) - f(x[i - 1]))) / (m * (x[i] - x[i - 1]))) - 2.0 * (f(x[i]) + f(x[i - 1]));
			if (maxrr < rr)
			{
				maxrr = rr;
				t = i;
			}
		}
		double xk = 0.5 * (x[t] + x[t - 1]) - ((f(x[t]) - f(x[t - 1])) / (2.0 * m));
		insertup(x, k, xk);
		k++;
	} while ((x[t] - x[t - 1] > ee) && (k - 2 <= n));
	minf = f(x[0]);
	minx = x[0];
	for (int i = 1; i < k; i++)
	{
		double ff = f(x[i]);
		if (minf > ff)
		{
			minf = ff;
			minx = x[i];
		}
	}
	cout << "Min f = " << minf << endl;
	cout << "Arg min f = " << minx << endl;
	delete[] x;
	cin >> k;
	return 0;
}