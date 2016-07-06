#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <limits>
#include <cassert>
#include <Eigen/Dense>

using namespace std;

typedef void(*xferMethod)(const int &i, const int &j, double &x, double &y);

const double pi = 4 * atan(1);

//D2Q9参数
const int D = 2;
const int Q = 9;
const int cx[Q] = { 0,1,0,-1,0,1,-1,-1,1 };
const int cy[Q] = { 0,0,1,0,-1,1,1,-1,-1 };
const double w[Q] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };
const int reverse_dir[Q] = { 0,3,4,1,2,7,8,5,6 };
const double cs_square = 1.0 / 3;

//流场尺寸与网格定义
const double radius = 0.5;
const double ratio = 50;
const int MeshType = 1;
const double min_gap = radius / 100;
const double min_angle = 1.5;
const int N = 120;
const int M = 360/min_angle;
const double cha = 1.2;
const double dt = 0.8*(radius*min_angle*pi / 180.0);

//流场初始变量
const double Ma = 0.1;
const int STEPS = 3000;
const double DENSITY = 1.0;
const double U = Ma*sqrt(cs_square);
const double V = 0.0;
const double Re = 200;
const double dynamic_vis = U*radius * 2 / Re;
const double vis = DENSITY*dynamic_vis;
const double tau = 3 * dynamic_vis / dt + 0.5;
const double I = 0.01;//紊流度

//残差
double err_vel = 1.0;
double err_rho = 1.0;

//气动力系数
double CL = 0.0;
double CD = 0.0;

class LBM_Node;
void calcResidue(const vector<vector<LBM_Node>> &mesh);
void calcAerodynamics(const vector<vector<LBM_Node>> &mesh);
void Output_Velocity(ofstream &fout, const vector<vector<LBM_Node>> &mesh);
void Output_Mesh(ofstream &fout, const vector<vector<LBM_Node>> &mesh);

void comp2phy_linear(const int &i, const int &j, double &x, double &y)
{
	double r = radius + min_gap*(i - 1);
	double theta = (min_angle*pi / 180)*j;
	x = r*cos(theta);
	y = r*sin(theta);
}

void comp2phy_nonuniform(const int &i, const int &j, double &x, double &y)
{
	double part = 1.0*(i - 1) / (N - 2);

	double r = radius + (ratio - 1)*radius*(1.0 - 1.0 / cha*atan((1.0 - part)*tan(cha)));
	double theta = (min_angle*pi / 180)*j;
	x = r*cos(theta);
	y = r*sin(theta);
}

vector<xferMethod> xfer = { comp2phy_linear, comp2phy_nonuniform };

class LBM_Node
{
private:
	//Helper function for performing operator=.
	//Adopt the copy-on-swap idiom.
	//Not exchange one's position (x,y)!
	void exchange(LBM_Node &rhs)
	{
		swap(density, rhs.density);
		swap(u, rhs.u);
		swap(v, rhs.v);

		swap(density_pre, rhs.density_pre);
		swap(u_pre, rhs.u_pre);
		swap(v_pre, rhs.v_pre);

		for (int i = 0; i < Q; i++)
		{
			swap(f[i], rhs.f[i]);
			swap(f_eq[i], rhs.f_eq[i]);
			swap(f_next[i], rhs.f_next[i]);
		}
	}

public:
	double x, y;

	double density;
	double u, v;
	
	double density_pre;
	double u_pre, v_pre;

	double f[Q], f_eq[Q], f_next[Q];

	//Constructors
	LBM_Node() = default;
	LBM_Node(double _den, double _u, double _v) :density(_den), u(_u), v(_v)
	{
		density_pre = 0;
		u_pre = 0;
		v_pre = 0;
	}
	LBM_Node(const LBM_Node &rhs) :x(rhs.x), y(rhs.y), density(rhs.density), u(rhs.u), v(rhs.v), density_pre(rhs.density_pre), u_pre(rhs.u_pre), v_pre(rhs.v_pre)
	{
		memcpy(f, rhs.f, Q * sizeof(double));
		memcpy(f_eq, rhs.f_eq, Q * sizeof(double));
		memcpy(f_next, rhs.f_next, Q * sizeof(double));
	}

	//Destructor
	~LBM_Node() = default;

	//operator= for two LBM_Node
	LBM_Node & operator=(LBM_Node rhs)
	{
		exchange(rhs);
		return *this;
	}

	//Calculate the equilibrim part of current density distribution.
	//All the directions of f_eq[] will be updated according to EDF.
	void calcDistribution_EQ()
	{
		for (int i = 0; i < Q; i++)
			calcDistribution_EQ_Dir(i);
	}

	//Get the equilibrim part of the density distribution from current local physical variables.
	//This is the kernel of the Lattice Boltzmann Method!
	double calcDistribution_EQ_Dir(int k)
	{
		double cu = cx[k] * u + cy[k] * v;
		double U2 = u*u + v*v;
		return f_eq[k] = w[k] * density*(1 + 3 * cu + 4.5*cu*cu - 1.5*U2);
	}

	//Update the previous records
	//Should be done first before calculating newly updated physical variables.
	void updatePrevRecord()
	{
		density_pre = density;
		u_pre = u;
		v_pre = v;
	}

	//Calculate the physical variables from current density distribution.
	//The conversion also reveals the physical characteristics of the density distribution.
	void calcPhysicalVar()
	{
		density = u = v = 0;
		for (int i = 0; i < Q; i++)
		{
			density += f[i];
			u += cx[i] * f[i];
			v += cy[i] * f[i];
		}
		u /= density;
		v /= density;
	}

	//Set the position of this point in 2D.
	//It should be noted that the operator= doesn't change one's position!
	void setPos(double _x, double _y)
	{
		x = _x;
		y = _y;
	}

	//Set new velocity components of this point
	void setVelocity(double _u, double _v)
	{
		u = _u;
		v = _v;
	}

	//Get the value of certain physical variable.
	//The results are all in SI.
	inline double getVelocity() const { return sqrt(pow(u, 2.0) + pow(v, 2.0)); }
	inline double getVelocityDiff() const { return sqrt(pow(u - u_pre, 2.0) + pow(v - v_pre, 2.0)); }
	inline double getDensityDiff() const { return fabs(density - density_pre); }
};

//所有格点
vector<vector<LBM_Node>> node(N, vector<LBM_Node>(M, LBM_Node(DENSITY, U, V)));

//每点1-8方向的A矩阵的第一行
vector<vector<vector<vector<double>>>> A(N, vector<vector<vector<double>>>(M, vector<vector<double>>(Q - 1, vector<double>(Q, 0.0))));

int main(int argc, char **argv)
{
	ofstream vel("velocity.dat");
	ofstream resi("residue.dat"); 
	ofstream ad("aerodynamics.dat"); 
	ofstream msh("mesh.dat");

	srand(time(NULL));

	//设置一定的紊流度
	for (int i = 2; i < node.size(); i++)
	{
		for (int j = 0; j < node[i].size(); j++)
		{
			//node[i][j].u += (I*U)*(-1 + 2 * (rand() % 100) *0.01);
			node[i][j].v += (I*U)*(-1 + 2 * (rand() % 100) *0.01);
		}
	}

	/*
	for (int i = 1; i < node.size(); i++)
	{
		cerr << "Layer " << i << ":" << endl;
		for (int j = 0; j < node[i].size(); j++)
			cerr << node[i][j].u << " " << node[i][j].v << endl;
	}
	*/

	//壁面速度
	for (int j = 0; j < M; j++)
		node[1][j].setVelocity(0, 0);

	//所有点的位置
	for (int i = 1; i < N; i++)
		for (int j = 0; j < M; j++)
			(*xfer[MeshType])(i, j, node[i][j].x, node[i][j].y);

	for (int j = 0; j < M; j++)
		node[0][j].setPos(2 * node[1][j].x - node[2][j].x, 2 * node[1][j].y - node[2][j].y);

	Output_Mesh(msh, node);

	//初始化f=f_eq
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			node[i][j].calcDistribution_EQ();
			memcpy(node[i][j].f, node[i][j].f_eq, Q*sizeof(double));
		}
	}

	//处理最外层的pre记录，迭代过程中始终保持最外层不变
	for (int j = 0; j < M; j++)
		node[N - 1][j].updatePrevRecord();

	//预先算出每点每个方向上的A矩阵的第一行
	int posX[Q], posY[Q];
	double deltaX[Q], deltaY[Q];
	Eigen::Matrix<double, Q, 6> S;

	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//周围点的标号
			for (int k = 0; k < Q; k++)
			{
				posX[k] = i + cx[k];
				posY[k] = (j + cy[k] + M) % M;
			}

			//计算1-8方向上的位形矩阵，由于0方向没有速度，LBM得出的f就是下一时间步的f
			for (int dir = 1; dir < Q; dir++)
			{
				//周围点做LBM后得到的中间点与(i,j)点之间的delta信息
				for (int k = 0; k < Q; k++)
				{
					deltaX[k] = node[posX[k]][posY[k]].x + cx[dir] * dt - node[i][j].x;
					deltaY[k] = node[posX[k]][posY[k]].y + cy[dir] * dt - node[i][j].y;
				}

				//构造矩阵S
				for (int k = 0; k < Q; k++)
				{
					S(k, 0) = 1;
					S(k, 1) = deltaX[k];
					S(k, 2) = deltaY[k];
					S(k, 3) = pow(deltaX[k], 2) / 2;
					S(k, 4) = pow(deltaY[k], 2) / 2;
					S(k, 5) = deltaX[k] * deltaY[k];
				}

				//当前点的A矩阵
				auto curA = ((S.transpose()*S).inverse()*S.transpose()).eval();

				//记下第一行
				for (int k = 0; k < Q; k++)
					A[i][j][dir - 1][k] = curA(0, k);
			}
		}
	}

	//Iterate
	int step = 0;

	while (step < STEPS)
	{
		cout << step << endl;

		//Update the virtual mesh inside the cylinder (noted as i = 0) for next computation,
		//and we don't change its x,y position.
		for (int j = 0; j < M; j++)
			node[0][j] = node[2][j];

		//Calculate f(P,t+dt) for each point P.
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 0; j < M; j++)
			{	
				//0方向，由于没有速度，直接LBM
				node[i][j].f_next[0] = node[i][j].f[0] + (1.0 / tau)*(node[i][j].f_eq[0] - node[i][j].f[0]);

				//1-8方向，LBM周围点后再做最小二乘
				for (int dir = 1; dir < Q; dir++)
				{
					//周围点的标号
					for (int k = 0; k < Q; k++)
					{
						posX[k] = i + cx[k];
						posY[k] = (j + cy[k] + M) % M;
					}

					double g[Q];
					for (int k = 0; k < Q; k++)
						g[k] = node[posX[k]][posY[k]].f[dir] + (1.0 / tau)*(node[posX[k]][posY[k]].f_eq[dir] - node[posX[k]][posY[k]].f[dir]);

					node[i][j].f_next[dir] = 0;
					for (int k = 0; k < Q; k++)
						node[i][j].f_next[dir] += A[i][j][dir - 1][k] * g[k];
				}
			}
		}		

		//When all the nodes get their density distribution in next time-step,update them all together.
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 0; j < M; j++)
			{
				memcpy(node[i][j].f, node[i][j].f_next, Q*sizeof(double));
				node[i][j].updatePrevRecord();
				node[i][j].calcPhysicalVar();
			}
		}

		//Keep the boundary condition.
		for (int j = 0; j < M; j++)
			node[1][j].setVelocity(0, 0);

		//Calculate the f_eq for newly updated physical variables.
		for (int i = 1; i < N - 1; i++)
			for (int j = 0; j < M; j++)
				node[i][j].calcDistribution_EQ();
		
		//The residues and aerodynamic coefficents
		calcResidue(node);
		calcAerodynamics(node);
		
		//Output some intermediate results
		resi << setw(6) << step << setw(16) << err_vel << setw(16) << err_rho << endl;
		ad << setw(6) << step << setw(16) << CD << setw(16) << CL << endl;

		++step;
	}

	//Output the streamlines and vorticity contours
	Output_Velocity(vel, node);

	vel.close();
	resi.close();
	ad.close();
	msh.close();
	return 0;
}

void calcResidue(const vector<vector<LBM_Node>> &mesh)
{
	double e1 = 0, e2 = 0;
	double r1=0, r2 = 0;
	
	for (int i = 1; i < mesh.size(); i++)
	{
		for (int j = 0; j < mesh[i].size(); j++)
		{
			e1 += mesh[i][j].getVelocityDiff();
			e2 += mesh[i][j].getVelocity();

			r1 += mesh[i][j].getDensityDiff();
			r2 += mesh[i][j].density;
		}
	}

	err_vel = e1 / e2;
	err_rho = r1 / r2;
}

void calcAerodynamics(const vector<vector<LBM_Node>> &mesh)
{
	double Fx = 0.0, Fy = 0.0;

	//P引起的法向力
	for (int i = 0; i < mesh[1].size(); i++)
	{
		int curIndex = i;
		int nextIndex = (i + 1) % mesh[1].size();

		double p1 = mesh[1][curIndex].density*cs_square;
		double p2 = mesh[1][nextIndex].density*cs_square;
		double p_av = (p1 + p2) / 2;
		double len = sqrt(pow(mesh[1][curIndex].x - mesh[1][nextIndex].x, 2) + pow(mesh[1][curIndex].y - mesh[1][nextIndex].y, 2));
		double mid_x = (mesh[1][curIndex].x + mesh[1][nextIndex].x) / 2;
		double mid_y = (mesh[1][curIndex].y + mesh[1][nextIndex].y) / 2;

		double f = -p_av*len;

		Fx += f*mid_x / sqrt(pow(mid_x, 2) + pow(mid_y, 2));
		Fy += f*mid_y / sqrt(pow(mid_x, 2) + pow(mid_y, 2));
	}

	//切向的粘性力
	for (int i = 0; i < mesh[1].size(); i++)
	{
		double sin_theta = mesh[1][i].y / sqrt(pow(mesh[1][i].x, 2) + pow(mesh[1][i].y, 2));
		double cos_theta = mesh[1][i].x / sqrt(pow(mesh[1][i].x, 2) + pow(mesh[1][i].y, 2));

		double vel1_t = -mesh[2][i].u*sin_theta + mesh[2][i].v*cos_theta;
		double vel2_t = -mesh[3][i].u*sin_theta + mesh[3][i].v*cos_theta;

		double d1 = sqrt(pow(mesh[2][i].x - mesh[1][i].x, 2) + pow(mesh[2][i].y - mesh[1][i].y, 2));
		double d2 = sqrt(pow(mesh[3][i].x - mesh[1][i].x, 2) + pow(mesh[3][i].y - mesh[1][i].y, 2));

		double a1 = (d2*d2*vel1_t - d1*d1*vel2_t) / ((d2 - d1)*d1*d2);

		double f = vis*a1*(radius*min_angle*pi / 180);

		Fx += f*(-sin_theta);
		Fy += f*cos_theta;
	}

	CD = Fx / (0.5*DENSITY*U*U);
	CL = Fy / (0.5*DENSITY*U*U);
}

void Output_Mesh(ofstream &fout, const vector<vector<LBM_Node>> &mesh)
{
	//按照有限元中的四边形单元输出
	fout << "TITLE = \"O-Grid for the cylinder\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "ZONE N = " << (mesh.size()-1)*mesh[0].size() << ", E = " << (mesh.size()-2)*mesh[0].size() << ", F = FEPOINT, ET = QUADRILATERAL" << endl;

	//先输出节点,虽然这里是以0为base，但是Tecplot中是以1为base！
	for (int i = 1; i < mesh.size(); i++)
		for (int j = 0; j < mesh[0].size(); j++)
			fout << mesh[i][j].x << " " << mesh[i][j].y << endl;

	//再输出每个单元中的节点的连接关系，注意顺序！
	for (int i = 1; i < mesh.size()-1; i++)
	{
		int minIndex = mesh[0].size() * (i - 1) + 1;
		int maxIndex = mesh[0].size() * i;

		int p[4] = { 0 };

		for (int j = minIndex; j <= maxIndex; j++)
		{
			p[0] = j;
			p[1] = j + 1 > maxIndex ? minIndex : j + 1;
			p[2] = p[1] + mesh[0].size();
			p[3] = p[0] + mesh[0].size();

			for (int k = 0; k < 4; k++)
				fout << p[k] << ' ';
			fout << endl;
		}
	}
}

void Output_Velocity(ofstream &fout, const vector<vector<LBM_Node>> &mesh)
{
	//按照有限元中的四边形单元输出
	fout << "TITLE = \"flow field velocity info\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"rho\", \"U\", \"V\"" << endl;
	fout << "ZONE N = " << (mesh.size() - 1)*mesh[0].size() << ", E = " << (mesh.size() - 2)*mesh[0].size() << ", F = FEPOINT, ET = QUADRILATERAL" << endl;

	//先输出节点,虽然这里是以0为base，但是Tecplot中是以1为base！
	for (int i = 1; i < mesh.size(); i++)
		for (int j = 0; j < mesh[0].size(); j++)
			fout << mesh[i][j].x << " " 
				 << mesh[i][j].y << " "
				 << mesh[i][j].density << " "
				 << mesh[i][j].u << " " 
				 << mesh[i][j].v 
				 << endl;

	//再输出每个单元中的节点的连接关系，注意顺序！
	for (int i = 1; i < mesh.size()-1; i++)
	{
		int minIndex = mesh[1].size()*(i-1) + 1;
		int maxIndex = mesh[1].size()*i;

		int p[4] = { 0 };

		for (int j = minIndex; j <= maxIndex; j++)
		{
			p[0] = j;
			p[1] = j + 1 > maxIndex ? minIndex : j + 1;
			p[2] = p[1] + mesh[0].size();
			p[3] = p[0] + mesh[0].size();

			for (int k = 0; k < 4; k++)
				fout << p[k] << ' ';
			fout << endl;
		}
	}
}
