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

ifstream mesh("../data/naca0012.grd");
ofstream conv("../result/convergence_history.dat");
ofstream flow("../result/flow_field.dat");
ofstream pd("../result/pressure_distribution.dat");

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
double s_inf = p_inf / pow(rho_inf, gama);
double ke_inf = 0.5*rho_inf*pow(Ma*c_inf, 2);

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

double err_rho = 1.0;
double err_vel = 1.0;
double err_p = 1.0;
double cl = 0;
double cd = 0;

int STEP = 5000;

class Point
{
	friend ostream &operator<<(ostream &fout, const Point &p);

private:
	//opertaor=的辅助函数
	void exchange(Point &rhs)
	{
		swap(coordinate, rhs.coordinate);
	}

public:
	vector<double> coordinate;

	//构造函数
	Point() = default;
	Point(double _x, double _y) :coordinate(vector<double>{_x, _y}) {}
	Point(double _x, double _y, double _z) :coordinate(vector<double>{_x, _y, _z}) {}
	Point(const vector<double> &_coord) :coordinate(_coord) {}
	Point(const Point &rhs) :coordinate(rhs.coordinate) {}
	
	//析构函数
	~Point() = default;
	
	//重载赋值运算符
	Point &operator=(Point rhs)
	{
		exchange(rhs);
		return *this;
	}

	//到另一点的几何距离
	double DistanceTo(const Point &rhs) const
	{
		assert(coordinate.size() == rhs.coordinate.size());

		double len = 0;
		for (int i = 0; i < coordinate.size(); i++)
			len += pow(coordinate[i] - rhs.coordinate[i], 2);
		return sqrt(len);
	}

	//到另一点的矢量
	vector<double> DeltaTo(const Point &rhs) const
	{
		assert(coordinate.size() == rhs.coordinate.size());

		vector<double> ans(coordinate.size(), 0);
		for (int i = 0; i < coordinate.size(); i++)
			ans[i] = rhs.coordinate[i] - coordinate[i];
		return ans;
	}
};

//Point类的<<运算符
ostream &operator<<(ostream &fout, const Point &p)
{
	for (auto c : p.coordinate)
		fout << setw(16) << c;

	return fout;
}

class Edge
{
	friend void convertToConservativeVar(vector<double> &primitive, vector<double> &conservative);
	friend void convertToPhysicalVar(vector<double> &conservative, vector<double> &primitive);

private:
	//opertaor=的辅助函数
	void exchange(Edge &rhs)
	{
		swap(start, rhs.start);
		swap(end, rhs.end);
		swap(len, rhs.len);
		swap(delta, rhs.delta);	
		swap(delta_n, rhs.delta_n);
		swap(c, rhs.c);
	}

	//根据无反射边界条件计算远场边
	//得到w_av,rho,u,v,p,c
	void calcFarBound()
	{
		double Vn_inf = u_inf*delta_n[0] + v_inf*delta_n[1];
		if (fabs(Vn_inf / c_inf) >= 1)
		{
			if (Vn_inf < 0)//超音速入流,边上的物理量和来流一致
			{
				w_av = w_inf;
				physicalVar = primitives_inf;
				c = c_inf;
			}
			else//超音速出流，边上的物理量和左侧的单元一致
			{
				w_av = cell[leftCell].w_next;
				convertToPhysicalVar(w_av, physicalVar);
				c = sqrt(gama*physicalVar[3] / physicalVar[0]);
			}
		}
		else
		{
			w_av = cell[leftCell].w_next;
			convertToPhysicalVar(w_av, physicalVar);
			c = sqrt(gama*physicalVar[3] / physicalVar[0]);
			double Vn_cell = physicalVar[1] * delta_n[0] + physicalVar[2] * delta_n[1];

			//Riemann Invariants
			double R_plus = Vn_cell + 2 * c / (gama - 1);
			double R_minus = Vn_inf - 2 * c_inf / (gama - 1);

			//由两个Riemann不变量得到边界上的法向速度与音速
			double Vn = 0.5*(R_plus + R_minus);
			c = 0.25*(gama - 1)*(R_plus - R_minus);

			if (Vn < 0)//亚音速入流
			{
				double s = s_inf;
				physicalVar[0] = pow(pow(c,2) / (gama*s), 1.0 / (gama - 1));
				physicalVar[1] = u_inf + (Vn - Vn_inf)*delta_n[0];
				physicalVar[2] = v_inf + (Vn - Vn_inf)*delta_n[1];
				physicalVar[3] = physicalVar[0] * pow(c, 2) / gama;

				convertToConservativeVar(physicalVar, w_av);
			}
			else//亚音速出流
			{
				double s = w_av[3] / pow(w_av[0], gama);
				physicalVar[0] = pow(pow(c, 2) / (gama*s), 1.0 / (gama - 1));
				physicalVar[1] += (Vn - Vn_cell)*delta_n[0];
				physicalVar[2] += (Vn - Vn_cell)*delta_n[1];
				physicalVar[3] = physicalVar[0] * pow(c, 2) / gama;

				convertToConservativeVar(physicalVar, w_av);
			}
		}
	}

public:
	int start, end;
	int leftCell, rightCell;
	double len;
	vector<double> delta;
	vector<double> delta_n;//单位外法向量(右手)
	
	vector<double> w_av;//两边单元的物理量的平均
	vector<double> physicalVar;//边上的物理量：rho,u,v,p
	double c;//当地声速
	double Z;//用于辅助计算
	vector<double> convective_flux;//该边上的对流通量

	vector<double> w_diff;//两边单元物理量之差,右-左，没有右则不算
	double ScalingFactor;//谱半径*边长+fabs(Z)
	double v;//激波探测器
	double eps2, eps4;//流场探测器
	vector<double> d2, d4;//components of the dissipation_flux
	
	//构造函数
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
	
	//析构函数
	~Edge() = default;

	//重载赋值运算
	Edge &operator=(Edge rhs)
	{
		exchange(rhs);
		return *this;
	}

	//取得中点
	Point getMidPoint() const
	{
		assert(point[start].coordinate.size() == point[end].coordinate.size());

		vector<double> mid_point(point[start].coordinate.size());
		for (int i = 0; i < mid_point.size(); i++)
			mid_point[i] = 0.5*(point[start].coordinate[i] + point[end].coordinate[i]);

		return Point(mid_point);
	}

	//计算该边上的w_av，根据不同类型做不同处理
	//进一步得到边上的rho,u,v,p,c,Z
	//注意：由于这只用于Runge-Kutta中的迭代，所以从每个cell中的w_next中取数据
	void calcBasicVariables()
	{
		if (rightCell >= 0)//内部边
		{
			for (int i = 0; i < w_av.size(); i++)
				w_av[i] = 0.5*(cell[leftCell].w_next[i] + cell[rightCell].w_next[i]);

			convertToPhysicalVar(w_av, physicalVar);
			c = sqrt(gama*physicalVar[3] / physicalVar[0]);
			Z = physicalVar[1] * delta[1] - physicalVar[2] * delta[0];
		}
		else if (rightCell == -1)//物面
		{
			for (int i = 0; i < w_av.size(); i++)
				w_av[i] = cell[leftCell].w_next[i];

			convertToPhysicalVar(w_av, physicalVar);
			c = sqrt(gama*physicalVar[3] / physicalVar[0]);
			Z = 0;
		}
		else//远场
		{
			calcFarBound();
			Z = physicalVar[1] * delta[1] - physicalVar[2] * delta[0];
		}
	}

	//计算该边上的对流通量Q
	void calcConvection()
	{
		convective_flux[0] = Z*physicalVar[0];
		convective_flux[1] = convective_flux[0] * physicalVar[1] + physicalVar[3] * delta[1];
		convective_flux[2] = convective_flux[0] * physicalVar[2] - physicalVar[3] * delta[0];
		convective_flux[3] = Z*(physicalVar[3] * (gama / (gama - 1)) + 0.5*physicalVar[0] * (pow(physicalVar[1], 2) + pow(physicalVar[2], 2)));
	}

	//将该边上的对流通量Q累加到相关联的cell
	void addConvectionToAdjcentCell()
	{
		for (int i = 0; i < convective_flux.size(); i++)
			cell[leftCell].Q[i] += convective_flux[i];

		if (rightCell >= 0)
		{
			for (int i = 0; i < convective_flux.size(); i++)
				cell[rightCell].Q[i] -= convective_flux[i];
		}
	}

	//计算与人工耗散相关的量：scaling_factor，v，eps2,eps4,w_diff(右-左),d2
	//在过程中将w_diff累加到相关联的单元
	//注意：耗散相关的量需要边的两侧都有cell,但此处并未显示要求，应在调用过程中体现
	//由于只用于Runge-Kutta中的迭代，所以用cell的w_next而不是w
	void calcDissipationRelatedVariables()
	{
		//人工耗散相关的量
		ScalingFactor = fabs(Z) + c*len;
		v = fabs((cell[rightCell].getPressure() - cell[leftCell].getPressure()) / (cell[rightCell].getPressure() + cell[leftCell].getPressure()));
		eps2 = k2*v;
		eps4 = max(0.0, k4 - eps2);
		for (int i = 0; i < w_diff.size(); i++)
			w_diff[i] = cell[rightCell].w_next[i] - cell[leftCell].w_next[i];
		for (int i = 0; i < d2.size(); i++)
			d2[i] = ScalingFactor*eps2*w_diff[i];

		//累加w_diff
		for (int i = 0; i < w_diff.size(); i++)
		{
			cell[leftCell].ddw[i] += w_diff[i];
			cell[rightCell].ddw[i] -= w_diff[i];
		}
	}

	//计算d4
	//注意：耗散相关的量需要边的两侧都有cell,但此处并未显示要求，应在调用过程中体现
	void calc_d4()
	{
		for (int i = 0; i < d4.size(); i++)
			d4[i] = -ScalingFactor*eps4*(cell[rightCell].ddw[i] - cell[leftCell].ddw[i]);
	}

	//将该边上的耗散通量D=d2+d4累加到相关联的cell,注意正负
	//注意：耗散相关的量需要边的两侧都有cell,但此处并未显示要求，应在调用过程中体现
	void addDissipationToAdjcentCell()
	{
		assert(d2.size() == d4.size());

		for (int i = 0; i < d2.size(); i++)
		{
			cell[leftCell].D[i] += (d2[i] + d4[i]);
			cell[rightCell].D[i] -= (d2[i] + d4[i]);
		}
	}
};

class Cell
{
	friend void convertToConservativeVar(vector<double> &primitive, vector<double> &conservative);
	friend void convertToPhysicalVar(vector<double> &conservative, vector<double> &primitive);

private:
	//opertaor=的辅助函数
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

void convertToPhysicalVar(vector<double> &conservative, vector<double> &primitive);
void convertToConservativeVar(vector<double> &primitive, vector<double> &conservative);
void Init();
void TimeMarching(int step);
void Output_FlowField();
void Output_PressureDistribution();
void RungeKutta();
void calcResidue();
void calcAerodynamics();
void RungeKutta_SubStep(int step);

int main(int argc, char **argv)
{
	//读入网格数据并初始化
	Init();

	//迭代
	for (int step = 0; step < STEP; step++)
		TimeMarching(step);

	//输出稳态结果
	Output_FlowField();
	Output_PressureDistribution();

	mesh.close();
	conv.close();
	flow.close();
	pd.close();
	return 0;
}

void convertToPhysicalVar(vector<double> &conservative, vector<double> &primitive)
{
	primitive[0] = conservative[0];
	primitive[1] = conservative[1] / conservative[0];
	primitive[0] = conservative[2] / conservative[0];
	primitive[0] = (gama - 1)*(conservative[3] - 0.5*conservative[0] * (pow(primitive[1], 2) + pow(primitive[2], 2)));
}

void convertToConservativeVar(vector<double> &primitive, vector<double> &conservative)
{
	conservative[0] = primitive[0];
	conservative[1] = primitive[0] * primitive[1];
	conservative[2] = primitive[0] * primitive[2];
	conservative[3] = primitive[3] / (gama - 1) + 0.5*primitive[0] * (pow(primitive[1], 2) + pow(primitive[2], 2));
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

	//初始化边上的衍生数据
	//TODO
}

void TimeMarching(int step)
{
	for (int i = 0; i < cnt_cell; i++)
		cell[i].w_next = cell[i].w;

	RungeKutta();

	for (int i = 0; i < cnt_cell; i++)
	{
		cell[i].w = cell[i].w_next;
		cell[i].updatePrevRecord();
		convertToPhysicalVar(cell[i].w, cell[i].physicalVar);
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

void Output_FlowField()
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

void Output_PressureDistribution()
{
	vector<pair<Point,double>> up, down;

	for (auto e : edge)
	{
		if (e.rightCell != -1)
			continue;

		if (fabs(e.delta_n[1]) > eps)//上翼面
			up.push_back(make_pair(e.getMidPoint(), cell[e.leftCell].getPressure()));
		else//下翼面
			down.push_back(make_pair(e.getMidPoint(), cell[e.leftCell].getPressure()));
	}

	for (auto elem : up)
		pd << elem.first << setw(16) << elem.second << endl;

	pd << endl;

	for (auto elem : down)
		pd << elem.first << setw(16) << elem.second << endl;
}

void RungeKutta()
{
	//计算每个cell的当地时间步长
	for (auto c : cell)
		c.localTimeStep = 0;

	for (auto e : edge)
	{
		cell[e.leftCell].localTimeStep += e.ScalingFactor;
		if (e.rightCell>=0)
			cell[e.rightCell].localTimeStep += e.ScalingFactor;
	}
	for (auto c : cell)
		c.localTimeStep = CFL*c.volum / c.localTimeStep;

	//4步Runge-Kutta迭代
	for (int i = 1; i <= RungeKutta_STEP; i++)
		RungeKutta_SubStep(i);
}

void calcResidue()
{
	double a = 0, b = 0;

	for (int i = 0; i < cnt_cell; i++)
	{
		a += cell[i].getDensityDiff();
		b += cell[i].getDensity();
	}
	err_rho = a / b;

	a = b = 0;
	for (int i = 0; i < cnt_cell; i++)
	{
		a += cell[i].getVelocityDiff();
		b += cell[i].getVelocity();
	}
	err_vel = a / b;

	a = b = 0;
	for (int i = 0; i < cnt_cell; i++)
	{
		a += cell[i].getPressureDiff();
		b += cell[i].getPressure();
	}
	err_p = a / b;
}

void calcAerodynamics()
{
	double Fx = 0, Fy = 0;
	double p = 0;

	for (auto e : edge)
	{
		if (e.rightCell != -1)
			continue;

		p = cell[e.leftCell].getPressure();

		Fx += p*e.delta_n[0];
		Fy += p*e.delta_n[1];
	}

	cd = Fx / ke_inf;
	cl = Fy / ke_inf;
}

void RungeKutta_SubStep(int step)
{
	//清零每个单元中的Q,D,ddw
	for (auto c : cell)
	{
		fill(c.Q.begin(), c.Q.end(), 0);
		fill(c.D.begin(), c.D.end(), 0);
		fill(c.ddw.begin(), c.ddw.end(), 0);
	}

	//求每个单元中的Q
	for (auto e : edge)
	{
		//得到每条边上的w,rho,u,v,p,c,Z
		e.calcBasicVariables();
		
		//得到每条边上的Q
		e.calcConvection();

		//每条边上的Q累加到关联的cell
		e.addConvectionToAdjcentCell();
	}

	//求每个单元中的D
	for (auto e : edge)
	{
		if (e.rightCell < 0)
			continue;

		//计算每条边上与人工耗散相关的量：scaling_factor，v，eps2, eps4, w_diff,d2
		//并在过程中将w_diff累加到相关联的单元
		e.calcDissipationRelatedVariables();
	}

	//每个单元的ddw都聚齐了，继续算每条边上的d4，然后累加到关联的cell
	for (auto e : edge)
	{
		if (e.rightCell < 0)
			continue;
		
		//计算d4
		e.calc_d4();

		//每条边上的D=d2+d4累加到关联的cell,注意正负
		e.addDissipationToAdjcentCell();
	}

	//Runge-Kutta更新w_next
	for (auto c : cell)
		for (int i = 0; i < c.w_next.size(); i++)
			c.w_next[i] = c.w[i] + RungeKutta_alpha[step] * c.localTimeStep*(-1.0*(c.Q[i] - c.D[i]) / c.volum);
}
