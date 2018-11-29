#include <iostream>
#include <list>
#include <queue>
#include <math.h>
#include <ctime>
#include <thread>
#include <mutex>
using namespace std;

#define pi 3.14159265358979323846

mutex m_lock;
mutex q_lock;
mutex p_lock;

//double f(double x)
//{
//	return sin(x) + sin(10.0 * x / 3.0);
//}

//double f(double x)
//{
//	double s = 0.0;
//	for (int k = 1; k <= 5; k++)
//	{
//		s += k * sin((k + 1) * x + k);
//	}
//	return -1.0 * s;
//}

//double f(double x)
//{
//	return (3.0 * x - 1.4) * sin(18.0 * x);
//}

//double f(double x)
//{
//	return -1.0 * (x + sin(x)) * exp(-1.0 * x * x);
//}

double f(double x)
{
	return sin(x) + sin(10.0 * x / 3.0) + log(x) - 0.84 * x + 3.0;
}

//double f(double x)
//{
//	return -1.0 * sin(2.0 * pi * x) * exp(-1.0 * x);
//}

//double f(double x)
//{
//	return (x * x - 5.0 * x + 6.0) / (x * x + 1.0);
//}

//double f(double x)
//{
//	return -1.0 * x + sin(3.0 * x) - 1.0;
//}

//double f(double x)
//{
//	return 2.0 * (x - 3.0) * (x - 3.0) + exp(x * x / 2.0);
//}

struct point
{
	double x;
	double z;
	point(double x_ = 0, double z_ = 0): x(x_), z(z_) {}
};

struct interval
{
	double R;
	point* lp;
	point* rp;
	interval(double R_ = 0.0, point* lp_ = NULL, point* rp_ = NULL): R(R_), lp(lp_), rp(rp_) {}
	interval(const interval& i)
	{
		R = i.R;
		lp = i.lp;
		rp = i.rp;
	}
	interval& operator=(const interval &i)
	{
		R = i.R;
		lp = i.lp;
		rp = i.rp;
		return *this;
	}
};

bool operator<(const interval& i1, const interval& i2) 
{
	return (i1.R < i2.R) ? true : false;
}

double Rfunc(const point &lp_, const point &rp_, double m)
{
	double dx = rp_.x - lp_.x;
	double dz = rp_.z - lp_.z;
	return (m * dx + dz * dz / (m * dx) - 2.0 * (rp_.z + lp_.z));
}

point* insertup_list(list<point> &p, point &xk)
{
	list<point>::iterator itl, itr;
	itl = itr = p.begin();
	while ((itr != p.end()) && (itr->x < xk.x))
	{
		itl = itr;
		itr++;
	}
	p.insert(itr, xk);
	itl++;
	return &(*itl);
}

void agp_potok(priority_queue<interval> &q, list<point> &p, interval &zk, double &m, bool &ff, double ee)
{
	list<point>::iterator itl, itr;
	double mold = m;
	double mm = 0.0;
	double r = 2.0;

	double xk = 0.5 * (zk.rp->x + zk.lp->x) - ((zk.rp->z - zk.lp->z) / (2.0 * m));
	point t;
	t.x = xk;
	t.z = f(xk);

	p_lock.lock();
	point *tt = insertup_list(p, t);
	p_lock.unlock();
	q_lock.lock();
	q.push(interval(Rfunc(*zk.lp, *tt, m), zk.lp, tt));
	q.push(interval(Rfunc(*tt, *zk.rp, m), tt, zk.rp));
	q_lock.unlock();

	ff = (zk.rp->x - zk.lp->x > ee);

	itr = itl = p.begin();
	itr++;
	while (itr != p.end())
	{
		double max = fabs((itr->z - itl->z) / (itr->x - itl->x));
		if (mm < max)
			mm = max;
		itr++;
		itl++;
	}
	m_lock.lock();
	if (mm > 0.0)
		m = r * mm;
	else
		m = 1.0;
	m_lock.unlock();

	if (mold != m)
	{
		q_lock.lock();
		while (!q.empty())
			q.pop();


		itr = itl = p.begin();
		itr++;
		while (itr != p.end())
		{
			q.push(interval(Rfunc(*itl, *itr, m), &(*itl), &(*itr)));
			itl++;
			itr++;
		}
		q_lock.unlock();
	}
}

int main()
{
	double ddd = f(0.0866);
	int v = 1;
	do
	{
		list<point> p;
		list<point>::iterator itl, itr;
		priority_queue<interval> q;
		interval zk1, zk2, zk3, zk4;

		double a, b, ee, m = -1.0, minf, minx;
		int k, n;
		int st, en;
		bool f1, f2, f3, f4;
		double mold = m;
		double mm = 0.0;
		double r = 2.0;

		cout << "Write left and right side" << endl;
		cin >> a >> b;
		cout << "Write number iterations" << endl;
		cin >> n;
		cout << "Write inaccuracy" << endl;
		cin >> ee;

		st = clock();

		point a0(a, f(a));
		point b0(b, f(b));
		p.push_back(a0);
		p.push_back(b0);

		k = 1;
		itr = itl = p.begin();
		itr++;
		while (itr != p.end())
		{
			double max = fabs((itr->z - itl->z) / (itr->x - itl->x));
			if (mm < max)
				mm = max;
			itr++;
			itl++;
		}
		if (mm > 0.0)
			m = r * mm;
		else
			m = 1.0;

		if (mold != m)
		{
			while (!q.empty())
				q.pop();

			itr = itl = p.begin();
			itr++;
			while (itr != p.end())
			{
				q.push(interval(Rfunc(*itl, *itr, m), &(*itl), &(*itr)));
				itl++;
				itr++;
			}
		}
		do
		{
			switch (k)
			{
			case 1:
			{
				zk1 = q.top();
				q.pop();
				thread p1(agp_potok, ref(q), ref(p), ref(zk1), ref(m), ref(f1), ee);
				p1.join();
				k++;
				break;
			}
			case 2:
			{
				zk1 = q.top();
				q.pop();
				zk2 = q.top();
				q.pop();
				thread p1(agp_potok, ref(q), ref(p), ref(zk1), ref(m), ref(f1), ee);
				thread p2(agp_potok, ref(q), ref(p), ref(zk2), ref(m), ref(f2), ee);
				p1.join();
				p2.join();
				k += 2;
				break;
			}
			case 3:
			{
				zk1 = q.top();
				q.pop();
				zk2 = q.top();
				q.pop();
				zk3 = q.top();
				q.pop();
				thread p1(agp_potok, ref(q), ref(p), ref(zk1), ref(m), ref(f1), ee);
				thread p2(agp_potok, ref(q), ref(p), ref(zk2), ref(m), ref(f2), ee);
				thread p3(agp_potok, ref(q), ref(p), ref(zk3), ref(m), ref(f3), ee);
				p1.join();
				p2.join();
				p3.join();
				k += 3;
				break;
			}
			default:
			{
				zk1 = q.top();
				q.pop();
				zk2 = q.top();
				q.pop();
				zk3 = q.top();
				q.pop();
				zk4 = q.top();
				q.pop();
				thread p1(agp_potok, ref(q), ref(p), ref(zk1), ref(m), ref(f1), ee);
				thread p2(agp_potok, ref(q), ref(p), ref(zk2), ref(m), ref(f2), ee);
				thread p3(agp_potok, ref(q), ref(p), ref(zk3), ref(m), ref(f3), ee);
				thread p4(agp_potok, ref(q), ref(p), ref(zk4), ref(m), ref(f4), ee);
				p1.join();
				p2.join();
				p3.join();
				p4.join();
				k += 4;
				break;
			}
			}
		} while (((f1) || ((f2) && (k > 1)) || ((f3) && (k > 2)) || ((f4) && (k > 3))) && (k < n));

		itl = p.begin();
		minf = itl->z;
		minx = itl->x;
		itl++;

		while (itl != p.end())
		{
			if (minf > itl->z)
			{
				minf = itl->z;
				minx = itl->x;
			}
			itl++;
		}

		en = clock();

		cout << "Time work = " << (double)(en - st) / 1000 << endl;
		cout << "Min f = " << minf << endl;
		cout << "Arg min f = " << minx << endl;
		cout << "Number iterations = " << k << endl;
		cout << "Want continue?" << endl;
		cout << "1 - yes" << endl;
		cout << "2 - no" << endl;
		cin >> v;
	} while (v == 1);
	return 0;
}