#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <cassert>

using namespace std;

const int Dimension = 2;

const double eps = 1e-8;
const double pi = 4 * atan(1);

int cnt_point = 0;
int cnt_edge = 0;
int cnt_cell = 0;

const double gama = 1.4;
double Ma = 0.8;
double angle_of_attack = 3;//degree
double p_inf = 101325.0;//Pa
double rho_inf = 1.225;//Kg/m^3
double c_inf = sqrt(gama*p_inf / rho_inf);//来流音速，无粘，用等熵关系式
double u_inf = Ma*c_inf*cos(angle_of_attack*pi / 180);
double v_inf = Ma*c_inf*sin(angle_of_attack*pi / 180);

vector<double> w_inf//来流的守恒量
{ 
	rho_inf,
	rho_inf*u_inf,
	rho_inf*v_inf,
	p_inf / (gama - 1) + rho_inf*(pow(u_inf,2) + pow(v_inf,2))*0.5 
};
vector<double> primitives_inf
{
	rho_inf,
	u_inf,
	v_inf,
	p_inf
};

double k2 = 0.9;
double k4 = 0.015;
double CFL = 1;

const int RungeKutta_STEP = 4;
const double RungeKutta_alpha[RungeKutta_STEP + 1] = { 0,1.0 / 4,1.0 / 3,1.0 / 2,1 };

int STEP = 5000;

class Point
{
private:
	void exchange(Point &rhs)
	{
		swap(coordinate, rhs.coordinate);
	}

public:
	vector<double> coordinate;

	Point() = default;
	Point(double _x, double _y) :coordinate(vector<double>{_x, _y}) {}
	Point(const Point &rhs) = default;
	~Point() = default;
	
	Point &operator=(Point rhs)
	{
		exchange(rhs);
		return *this;
	}

	double DistanceTo(const Point &rhs) const
	{
		assert(coordinate.size() == rhs.coordinate.size());

		double len = 0;
		for (int i = 0; i < coordinate.size(); i++)
			len += pow(coordinate[i] - rhs.coordinate[i], 2);
		return sqrt(len);
	}

	vector<double> DeltaTo(const Point &rhs) const
	{
		assert(coordinate.size() == rhs.coordinate.size());

		vector<double> ans(coordinate.size(), 0);
		for (int i = 0; i < coordinate.size(); i++)
			ans[i] = rhs.coordinate[i] - coordinate[i];
		return ans;
	}
};

class Edge
{
private:
	void exchange(Edge &rhs)
	{
		swap(start, rhs.start);
		swap(end, rhs.end);
		swap(len, rhs.len);
		swap(delta, rhs.delta);	
		swap(delta_n, rhs.delta_n);
		swap(c, rhs.c);
	}

public:
	int start, end;
	int leftCell, rightCell;
	double len;
	vector<double> delta;
	vector<double> delta_n;//单位外法向量(右手)
	
	vector<double> average_w;//两边单元的物理量的平均
	vector<double> diff_w;//两边单元物理量之差
	double Z;//用于辅助计算
	double c;//当地声速
	double scale_factor;//谱半径
	double v;//激波探测器


	Edge() = default;
	Edge(const vector<Point> &pts, int _start, int _end, int _lc, int _rc) 
		:start(_start),
		end(_end), 
		leftCell(_lc), 
		rightCell(_rc),
		c(0),
		delta(pts[_start].DeltaTo(pts[_end]))
	{
		double length = 0;
		for (int i = 0; i < delta.size(); i++)
			length += pow(delta[i], 2);
		len = sqrt(length);

		delta_n = vector<double>{ delta[1] / len,-delta[0] / len };
	}
	Edge(const Edge &rhs) = default;
	~Edge() = default;

	Edge &operator=(Edge rhs)
	{
		exchange(rhs);
		return *this;
	}

};

class Cell
{
private:
	void exchange(Cell &rhs)
	{
		swap(point, rhs.point);
		swap(volum, rhs.volum);
		swap(w, rhs.w);
		swap(w_next, rhs.w_next);
		swap(physicalVar, rhs.physicalVar);
		swap(pre_physicalVar, rhs.pre_physicalVar);
	}

public:
	vector<int> point;
	double volum;
	vector<double> w;//rho,rho*u,rho*v,rho*E
	vector<double> w_next;
	vector<double> physicalVar;//rho,u,v,p
	vector<double> pre_physicalVar;

	double localTimeStep;
	vector<double> Q;
	vector<double> D;
	vector<double> ddw;

	Cell() = default;
	Cell(const vector<int> &_pts) :point(_pts), volum(0), w(w_inf), physicalVar(primitives_inf), w_next(vector<double>(4, 0)), pre_physicalVar(vector<double>(4, 0)) {}
	Cell(const Cell &rhs) = default;
	~Cell() = default;

	Cell &operator=(Cell rhs)
	{
		exchange(rhs);
		return *this;
	}

	void updatePrevRecord()
	{
		pre_physicalVar = physicalVar;
	}

	void calcPhysicalVar()
	{
		physicalVar[0] = w[0];
		physicalVar[1] = w[1]/w[0];
		physicalVar[0] = w[2]/w[0];
		physicalVar[0] = (gama - 1)*(w[3] - 0.5*w[0] * (pow(physicalVar[1], 2) + pow(physicalVar[2], 2)));
	}
	

	double getDensity() const
	{
		return physicalVar[0];
	}
	double getDensityDiff() const
	{
		return fabs(physicalVar[0] - pre_physicalVar[0]);
	}
	
	double getVelocity() const 
	{ 
		return sqrt(pow(physicalVar[1], 2) + pow(physicalVar[2], 2)); 
	}
	double getVelocityDiff() const 
	{ 
		return sqrt(pow(physicalVar[1] - pre_physicalVar[1], 2) + pow(physicalVar[2] - pre_physicalVar[2], 2));
	}

	double getPressure() const
	{
		return physicalVar[3];
	}
	double getPressureDiff() const
	{
		return fabs(physicalVar[3] - pre_physicalVar[3]);
	}
};

vector<Point> point;
vector<Edge> edge;
vector<Cell> cell;

double err_rho = 1.0;
double err_vel = 1.0;
double err_p = 1.0;
double cl = 0;
double cd = 0;

ifstream mesh("../data/naca0012.grd");
ofstream conv("../result/convergence_history.dat");
ofstream flow("../result/flow_field.dat");

void Init();
void TimeMarching(int step);
void Output();

int main(int argc, char **argv)
{
	//读入网格数据并初始化
	Init();

	//迭代
	for (int step = 0; step < STEP; step++)
		TimeMarching(step);

	//输出稳态结果
	Output();

	mesh.close();
	conv.close();
	flow.close();
	return 0;
}

void Init()
{
	//输入网格，在构造函数中设置初值
	mesh >> cnt_point >> cnt_edge >> cnt_cell;

	double x, y;
	for (int i = 0; i < cnt_point; i++)
	{
		mesh >> x >> y;
		point.push_back(Point(x, y));
	}

	int a, b, lc, rc;
	for (int i = 0; i < cnt_edge; i++)
	{
		mesh >> a >> b >> lc >> rc;
		edge.push_back(Edge(point, a - 1, b - 1, lc - 1, rc > 0 ? rc - 1 : rc));
	}

	vector<int> pts(3, -1);
	for (int i = 0; i < cnt_cell; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			mesh >> pts[k];
			--pts[k];
		}
		cell.push_back(Cell(pts));
	}

	for (int i = 0; i < cnt_cell; i++)
		mesh >> cell[i].volum;

	//初始化边上的衍生数据，并计算当地时间步长
	//TODO
}

void TimeMarching(int step)
{
	for (int i = 0; i < cnt_cell; i++)
		cell[i].w_next = cell[i].w;

	RungeKutta(cell);

	for (int i = 0; i < cnt_cell; i++)
	{
		cell[i].w = cell[i].w_next;
		cell[i].updatePrevRecord();
		cell[i].calcPhysicalVar();
	}

	calcResidue();
	calcAerodynamics();

	conv << setw(8) << step
		<< setw(16) << err_rho
		<< setw(16) << err_vel
		<< setw(16) << err_p
		<< setw(16) << cl
		<< setw(16) << cd
		<< endl;
}

void Output()
{
	flow << "TITLE = \"Flow Field for NACA0012\"" << endl;
	flow << "VARIABLES = \"X\", \"Y\", \"DENSITY\", \"U\", \"V\", \"P\"" << endl;
	flow << "ZONE  N = " << cnt_point << ", E = " << cnt_cell << ", ZONETYPE = FETRIANGLE" << endl;
	flow << "DATAPACKING = BLOCK" << endl;
	flow << "VARLOCATION = ([1-2] = NODAL, [3-6] = CELLCENTERED)" << endl;

	//节点X,Y坐标
	for (int k = 0; k < 2; k++)
		for (int i = 0; i < cnt_point; i++)
			flow << point[i].coordinate[k] << endl;
	
	//每个单元的DENSITY,U,V,P
	for (int k = 0; k < 4; k++)
		for (int i = 0; i < cnt_cell; i++)
			flow << cell[i].physicalVar[k] << endl;

	//单元与节点的对应关系
	for (int i = 0; i < cnt_cell; i++)
	{
		for (int k = 0; k < 3; k++)
			flow << cell[i].point[k] + 1 << ' ';
		flow << endl;
	}

}

void RungeKutta(vector<Cell> &elem)
{
	for (int i = 1; i <= RungeKutta_STEP; i++)
		RungeKutta_SubStep(elem, i);

	//计算下一时间步每个单元的当地时间步长
	//TODO
}

void calcResidue(vector<Cell> &elem)
{
	double a = 0, b = 0;

	for (int i = 0; i < cnt_cell; i++)
	{
		a += elem[i].getDensityDiff();
		b += elem[i].getDensity();
	}
	err_rho = a / b;

	a = b = 0;
	for (int i = 0; i < cnt_cell; i++)
	{
		a += elem[i].getVelocityDiff();
		b += elem[i].getVelocity();
	}
	err_vel = a / b;

	a = b = 0;
	for (int i = 0; i < cnt_cell; i++)
	{
		a += elem[i].getPressureDiff();
		b += elem[i].getPressure();
	}
	err_p = a / b;
}

void calcAerodynamics(vector<Cell> &elem)
{

}

void RungeKutta_SubStep(vector<Cell> &elem, int step)
{
	//计算每条边的平均物理量、物理量的差、Z、当地音速
	//将1/scale_factor累加到每个cell中，用来计算当地时间步长
	//将每条边上的ddw累加到关联的两个cell，用来计算人工粘性项
	for (int i = 0; i < cnt_edge; i++)
	{
		const int n = cell[edge[i].leftCell].physicalVar.size();
		if (edge[i].rightCell >= 0)
		{
			for (int k = 0; k < n; i++)
			{
				edge[i].average_w[k] = 0.5*(cell[edge[i].leftCell].w[k] + cell[edge[i].rightCell].w[k]);
				edge[i].diff_w[k] = cell[edge[i].leftCell].w[k] - cell[edge[i].rightCell].w[k];
			}
		}
		else if (edge[i].rightCell == -1)
		{
			edge[i].average_w = edge[i].diff_w = cell[edge[i].leftCell].physicalVar;
		}
		else
		{
			edge[i].average_w = edge[i].diff_w = cell[edge[i].leftCell].physicalVar;
		}
	}

	for (int i = 0; i < cnt_edge; i++)
	{
		edge[i].Z = fabs(edge[i].average_w[1] * edge[i].delta[1] - edge[i].average_w[2] * edge[i].delta[0]);
		edge[i].c = sqrt(gama*edge[i].average_w[3] / edge[i].average_w[0]);
		edge[i].scale_factor = edge[i].Z + edge[i].c*edge[i].len;
		edge[i].v = fabs((cell[edge[i].leftCell].getPressure() - cell[edge[i].rightCell].getPressure()) / (cell[edge[i].leftCell].getPressure() + cell[edge[i].rightCell].getPressure()));
	}

	for (int i = 0; i < cnt_edge; i++)
	{
		cell[edge[i].leftCell].localTimeStep += 1 / edge[i].scale_factor;
		for (int k = 0; k < cell[edge[i].leftCell].ddw.size();k++)
			cell[edge[i].leftCell].ddw[k] += -(edge[i].diff_w[k]);
		
		if (edge[i].rightCell >= 0)
		{
			cell[edge[i].rightCell].localTimeStep += 1 / edge[i].scale_factor;
			for (int k = 0; k < cell[edge[i].rightCell].ddw.size(); k++)
				cell[edge[i].rightCell].ddw[k] += edge[i].diff_w[k];
		}
	}



}

