/**
 * Delaunay Triangulation for the estimation of unsteady aerodynamics.
 * The program is based on the CGAL, which is a powerful computational geometry library.
 * 
 */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef Delaunay::Point Point;

const double pi = atan(1) * 4;
const size_t cntTask = 4;

size_t cnt = 0;
std::vector< std::pair<Point, unsigned> > samplePoint;
std::vector<double> sampleData;

double am[cntTask] = { 0, -10, 10, 0 };
double aa[cntTask] = { -24, -5, -5, -5 };
double freq[cntTask] = { 0.35, 0.2, 0.2, 0.2 };

void read_sample();
int work(const Delaunay &DT, double alpha_m, double alpha_a, double f, std::ofstream &out, int step);
inline double calcTriangleArea(const Point &p1, const Point &p2, const Point &p3);
inline double calcSegmentDistance(const Point &p1, const Point &p2);

double handle_three(Point p[], int index[]);
double handle_two(Point p[], int index[]);
double handle_one(Point p[], int index[]);

int main()
{
	read_sample();

	Delaunay T;
	T.insert(samplePoint.begin(), samplePoint.end());

	for (int i = 0; i < cntTask; i++)
	{
		std::stringstream ofn;
		ofn << "../results/Case" << i + 1 << ".txt";
		std::ofstream fout(ofn.str());
		if (!fout)
			throw "Failed to open output file!\n";

		work(T, am[i], aa[i], freq[i], fout, 80);

		fout.close();
	}

	return 0;
}

void read_sample()
{
	const double alpha_a = -24.0;

	char tmp[1024] = { 0 };
	double dataset[20] = { 0 };

	std::ifstream fin;
	for (int i = 1; i <= 7; i++)
	{
		double f = i*0.1;
		std::stringstream file_name;
		file_name << "../data/ba" << i << ".dat";
		fin.open(file_name.str());
		if (!fin)
			throw "Open data file failed!\n";

		int n = 0, k = 0;
		fin >> n >> k;

		fin.getline(tmp, sizeof(tmp));
		fin.getline(tmp, sizeof(tmp));

		for (int j = 0; j < n; j++)
		{
			for (int l = 1; l <= k; l++)
				fin >> dataset[l];

			double Cn = dataset[12];
			double psij = dataset[18];

			double psij_d = 2 * pi * f * std::sqrt(std::pow(alpha_a, 2) - std::pow(psij, 2));
			if (j >= n / 2)
				psij_d = -psij_d;

			samplePoint.push_back(std::make_pair(Point(psij, psij_d), cnt));
			sampleData.push_back(Cn);

			++cnt;
		}
		fin.close();
	}

	//for graphic view and test
	std::ofstream f2("../results/raw_2D.dat");
	std::ofstream f3("../results/raw_3D.dat");
	for (int i = 0; i < cnt; i++)
	{
		f2 << samplePoint[i].first.x() << ' ' << samplePoint[i].first.y() << std::endl;
		f3 << samplePoint[i].first.x() << ' ' << samplePoint[i].first.y() << ' ' << sampleData[i] << std::endl;
	}
	f2.close();
	f3.close();
}

int work(const Delaunay &DT, double alpha_m, double alpha_a, double f, std::ofstream &out, int step = 50)
{
	out.setf(std::ios::right, std::ios::adjustfield);
	out << std::setw(4) << "Index"
		<< std::setw(16) << "Psi_j"
		<< std::setw(16) << "d_Psi_j"
		<< std::setw(16) << "Cn"
		<< std::endl;

	for (int i = 0; i < step; i++)
	{
		double t = i*1.0 / (f*step);
		const double psi_j = alpha_m + alpha_a*cos(2 * pi*f*t);
		const double d_psi_j = -2 * pi*f*alpha_a*sin(2 * pi*f*t);

		Point p[4];
		int index[4];
		
		p[0] = Point(psi_j, d_psi_j);
		auto fh = DT.locate(p[0]);

		int validCnt = 3;
		for (int j = 1; j <= 3; j++)
		{
			auto cur_vit = fh->vertex(j - 1);
			p[j] = cur_vit->point();
			index[j] = cur_vit->info();
			if (index[j] < 0 || index[j] >= cnt)
				--validCnt;
		}

		double cur_ans = 0;
		if (validCnt == 3)
			cur_ans = handle_three(p, index);
		else if (validCnt == 2)
			cur_ans = handle_two(p, index);
		else if (validCnt == 1)
			cur_ans = handle_one(p, index);
		else
			continue;

		out.setf(std::ios::right, std::ios::adjustfield);
		out << std::setw(4) << i
			<< std::setw(16) << p[0].x()
			<< std::setw(16) << p[0].y()
			<< std::setw(16) << cur_ans
			<< std::endl;
	}
	return 0;
}

double handle_three(Point p[],int index[])
{
	double s[4];
	s[0] = calcTriangleArea(p[1], p[2], p[3]);
	s[1] = calcTriangleArea(p[0], p[2], p[3]);
	s[2] = calcTriangleArea(p[1], p[0], p[3]);
	s[3] = calcTriangleArea(p[1], p[2], p[0]);

	double esti_cn = 0;
	for (int j = 1; j <= 3; j++)
		esti_cn += s[j] * sampleData[index[j]];

	esti_cn /= s[0];

	return esti_cn;
}

double handle_two(Point p[], int index[])
{
	int cnt = 1, t[3];

	for (int i = 1; i <= 3; i++)
		if (index[i] >= 0)
			t[cnt++] = i;

	double l[3];
	l[0] = calcSegmentDistance(p[t[1]], p[t[2]]);
	l[1] = calcSegmentDistance(p[t[2]], p[0]);
	l[2] = calcSegmentDistance(p[0], p[t[1]]);

	double esti_cn = 0;
	for (int j = 1; j <= 2; j++)
		esti_cn += l[j] * sampleData[index[t[j]]];

	esti_cn /= l[0];

	return esti_cn;
}

double handle_one(Point p[], int index[])
{
	int i = 0;
	for (i = 1; i <= 3; i++)
		if (index[i] >= 0)
			break;

	return sampleData[index[i]];
}

inline double calcTriangleArea(const Point &p1, const Point &p2, const Point &p3)
{
	double x[4], y[4];
	x[1] = p1.x(); y[1] = p1.y();
	x[2] = p2.x(); y[2] = p2.y();
	x[3] = p3.x(); y[3] = p3.y();

	return 0.5*fabs((x[2] * y[3] + x[1] * y[2] + x[3] * y[1]) - (x[2] * y[1] + x[1] * y[3] + x[3] * y[2]));
}

inline double calcSegmentDistance(const Point &p1, const Point &p2)
{
	return std::sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2));
}
