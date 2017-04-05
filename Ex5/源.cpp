#include "CImg.h"
#include "canny_class.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm> 
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_linalg.h>

using namespace std;
using namespace cimg_library;

int n = 8;
double pi = 3.14159265358979323846;
vector<vector<double>> ver;
vector<vector<double>>ver_order;

//进行霍夫变换找到四个顶点
CImg<int> hough(CImg<float> img, CImg<float> src, double val)
{
	//设置阈值参数
	const int angle_num = 180;
	int range_thre = 5;
	double fun_thre = 0.8;
	int point_thre = 3;
	double ang = pi / angle_num;
	int width = img._width, height = img._height;
	int X = width / 2, Y = height / 2;

	const int d = int(round(sqrt(pow(width, 2) + pow(height, 2))));

	//初始化hough空间存储矩阵
	vector<vector<int>> s;
	s.resize(angle_num);
	for (int i = 0; i < angle_num; i++)
	{
		s[i].resize(2 * d);
	}

	//计算hough空间矩阵
	int max = 0;
	cimg_forXY(img, x, y)
	{
		if (img(x, y) != 0)
		{
			for (int i = 0; i < angle_num; i++)
			{
				int p = int(round((x - X)*cos(i*ang) + (y - Y)*sin(i*ang)));
				p += d;
				if (p < 0 || (p >= 2 * d))
					continue;
				else
				{
					s[i][p]++;
					if (s[i][p] > max)
						max = s[i][p];
				}
			}
		}
	}

	//寻找hough空间中线的交点，且其交点次数大于阈值，并是一定区域内的最大值，以免受到光栅化的影响而有相似的值
	vector<pair<int, int>> lines;
	int thre = int(max*val);
	for (int i = 0; i < angle_num; i++)
	{
		for (int j = 0; j < 2 * d; j++)
		{
			bool flag = 1;
			if (s[i][j] > thre)
			{
				for (int k = 0; k < int(lines.size()); k++)
				{
					if ((abs(lines[k].first - i) < range_thre || abs((angle_num - lines[k].first) + i) < range_thre) && abs(lines[k].second - j) < range_thre)
					{
						if (s[i][j] > s[lines[k].first][lines[k].second])
						{
							lines[k].first = i;
							lines[k].second = j;
						}
						flag = 0;
					}
				}
				if (flag)
				{
					lines.push_back(make_pair(i, j));
				}
			}
		}
	}

	//由交点得到原空间的线性方程
	vector<vector<double>> fun;
	cout << "lines : " << endl;
	for (int i = 0; i < int(lines.size()); i++)
	{
		int row = lines[i].first;
		int col = lines[i].second;
		double dy = sin(row*ang), dx = cos(row*ang);
		cout << dx << "x + " << dy << "y = " << col - d + Y*dy + X*dx << endl;
		vector<double> t;
		t.push_back(dx);
		t.push_back(dy);
		t.push_back(col - d + Y*dy + X*dx);
		fun.push_back(t);
	}

	//计算交点，且在一定区域内最准确，以免重复
	CImg<int> dst(src);
	unsigned char red[] = { 255,0,0 };
	unsigned char green[] = { 0,255,0 };
	//vector<vector<double>> ver;
	cimg_forXY(dst, x, y)
	{
		int sum = 0;
		double error = 0;
		for (int i = 0; i<int(fun.size()); i++)
		{
			double re = x*fun[i][0] + y*fun[i][1];
			if (re<fun[i][2] + fun_thre&&re>fun[i][2] - fun_thre)
			{
				sum++;
				error += abs(fun[i][2] - re);
				dst(x, y, 0) = 0; dst(x, y, 1) = 255; dst(x, y, 0) = 0;
			}
		}
		if (sum >= 2)
		{		
			bool flag = 1;
			for (int i = 0; i<int(ver.size()); i++)
			{
				if (abs(x - ver[i][0]) < point_thre&&abs(y - ver[i][1]) < point_thre)
				{
					if (error < ver[i][2])
					{
						ver[i][0] = x;
						ver[i][1] = y;
						ver[i][2] = error;
					}
					flag = 0;
				}
			}
			if (flag)
			{
				vector<double> t;
				t.push_back(x);
				t.push_back(y);
				t.push_back(error);
				ver.push_back(t);
			}
		}
	}

	//画点
	cout << "vertex : " << endl;
	for (int i = 0; i<int(ver.size()); i++)
	{
		cout << "( " << int(ver[i][0]) << " , " << int(ver[i][1]) << " )" << endl;
		dst.draw_circle(int(ver[i][0]), int(ver[i][1]), 6, red, 1);
	}
	cout << endl;
	return dst;
}
//将确定的4个点按照左上，左下，右下，右上的顺序排列
void order(CImg<float>img) {
	int width = img.width();
	vector<double>d00;//到（0,0）点的距离
	vector<double>dw0;//到（width,0）点的距离
	for (int i = 0; i<int(ver.size()); i++)
	{
		double d1 = 0, d2 = 0;
		d1 = sqrt(ver[i][0] * ver[i][0] + ver[i][1] * ver[i][1]);
		d2 = sqrt((ver[i][0] - width)*(ver[i][0] - width) + ver[i][1] * ver[i][1]);
		d00.push_back(d1);
		dw0.push_back(d2);
	}
	int left_up, left_down, right_up, right_down;
	vector<double>::iterator d00_max, d00_min, dw0_max, dw0_min;
	d00_max = std::max_element(begin(d00), end(d00));
	d00_min = std::min_element(begin(d00), end(d00));
	dw0_max = std::max_element(begin(dw0), end(dw0));
	dw0_min = std::min_element(begin(dw0), end(dw0));
	left_up = distance(begin(d00), d00_min);
	left_down = distance(begin(dw0), dw0_max);
	right_up = distance(begin(dw0), dw0_min);
	right_down = distance(begin(d00), d00_max);
	cout << left_up << " " << left_down << " " << right_down << " " << right_up << endl;
	//vector<vector<double>>ver_order;
	vector<double>temp;
	temp.push_back(ver[left_up][0]);
	temp.push_back(ver[left_up][1]);
	ver_order.push_back(temp);
	temp.clear();
	temp.push_back(ver[left_down][0]);
	temp.push_back(ver[left_down][1]);
	ver_order.push_back(temp);
	temp.clear();
	temp.push_back(ver[right_down][0]);
	temp.push_back(ver[right_down][1]);
	ver_order.push_back(temp);
	temp.clear();
	temp.push_back(ver[right_up][0]);
	temp.push_back(ver[right_up][1]);
	ver_order.push_back(temp);
	temp.clear();
	for (int i = 0; i < 4; i++) {
		cout << ver_order[i][0] << "  " << ver_order[i][1] << endl;
	}
}
//系数矩阵初始化
void Init(vector<vector<double> > & A)
{
	double u1, u2, u3, u4, v1, v2, v3, v4;
	double x1 = 0, x2 = 0, x3 = 500, x4 = 500, y1 = 0, y2 = 707, y3 = 707, y4 = 0;
	//u1 = 75; u2 = 69; u3 = 263; u4 = 271;
	//v1 = 59; v2 = 355; v3 = 351; v4 = 77;
	u1 = ver_order[0][0]; u2 = ver_order[1][0]; u3 = ver_order[2][0]; u4 = ver_order[3][0];
	v1 = ver_order[0][1]; v2 = ver_order[1][1]; v3 = ver_order[2][1]; v4 = ver_order[3][1];
	double a[8][9] = {
		x1,y1,1,0,0,0,-x1*u1,-y1*u1,u1,
		x2,y2,1,0,0,0,-x2*u2,-y2*u2,u2,
		x3,y3,1,0,0,0,-x3*u3,-y3*u3,u3,
		x4,y4,1,0,0,0,-x4*u4,-y4*u4,u4,
		0,0,0,x1,y1,1,-x1*v1,-y1*v1,v1,
		0,0,0,x2,y2,1,-x2*v2,-y2*v2,v2,
		0,0,0,x3,y3,1,-x3*v3,-y3*v3,v3,
		0,0,0,x4,y4,1,-x4*v4,-y4*v4,v4 };
	A.clear();
	for (int i = 0; i < n; i++)
	{
		vector<double> temp(n + 1);
		for (int j = 0; j < n + 1; j++)
		{
			temp[j] = a[i][j];
		}
		A.push_back(temp);
		temp.clear();
	}
}
//输出解向量
void DisplaySolver(const vector<double> & solve)
{
	int n = solve.size();
	for (int i = 0; i < n; i++)
	{
		cout << "x" << i << " = " << solve[i] << endl;
	}
}
//用gsl库函数解方程，A为系数矩阵，x为解向量
void GSLSolve(const vector<vector<double> > & A, vector<double> & x)
{
	x.clear();
	int n = A.size();
	int m = A[0].size();
	x.resize(n);
	double * a_data;
	double * b_data;
	try
	{
		a_data = new double[n*n];
		b_data = new double[n];
	}
	catch (bad_alloc)
	{
		cout << "内存分配失败！" << endl;
		return ;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a_data[i*n + j] = A[i][j];
		}
		b_data[i] = A[i][m - 1];
	}
	gsl_matrix_view gm = gsl_matrix_view_array(a_data, n, n);
	gsl_vector_view gb = gsl_vector_view_array(b_data, n);
	gsl_vector *sx = gsl_vector_alloc(n);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(&gm.matrix, p, &s);
	gsl_linalg_LU_solve(&gm.matrix, p, &gb.vector, sx);
	for (int i = 0; i < n; i++)
	{
		x[i] = gsl_vector_get(sx, i);
	}
	gsl_permutation_free(p);
	gsl_vector_free(sx);
	delete[] a_data;
	delete[] b_data;
}
int main()
{
	//设置参数
	string filename[] = { "dataset (1).bmp", "dataset (2).bmp", "dataset (3).bmp", "dataset (4).bmp", "dataset (5).bmp", "dataset (6).bmp" };
	double thre_val[] = { 0.5,0.6,0.6,0.5,0.5,0.5 };
	float sigma[] = { 6.0f, 6.0f, 9.0f, 6.0f, 6.0f, 6.0f };
	float threshold[] = { 3.5f, 3.5f, 5.0f, 3.5f, 3.5f, 3.5f };
	//先进行canny处理，再进行hough变换，并存在result文件夹中
	for (int ni = 0; ni < sizeof(filename)/sizeof(string); ni++)
	{
		cout << endl << filename[ni] << " : " << endl;
		canny_img img("Dataset\\" + filename[ni], sigma[ni], threshold[ni]);
		CImg<float> c_img=img.CannyDiscrete();
		CImg<float> origin(("Dataset\\" + filename[ni]).c_str());
		hough(c_img, origin, thre_val[ni]);
		CImg<float> goal(500, 707, 1, 3);
		order(origin);
		vector<vector<double>> A;
		//cout.width(10);
		Init(A);
		vector<double> x(A.size()); //解向量
		GSLSolve(A, x);
		DisplaySolver(x);
		int _x, _y;
		cimg_forXY(goal, a, b) {
			_x = int((x[0] * a + x[1] * b + x[2]) / (x[6] * a + x[7] * b + 1));
			_y = int((x[3] * a + x[4] * b + x[5]) / (x[6] * a + x[7] * b + 1));
			goal(a, b, 0) = origin(_x, _y, 0);
			goal(a, b, 1) = origin(_x, _y, 1);
			goal(a, b, 2) = origin(_x, _y, 2);
		}
		//(origin, goal).display();
		goal.save(("result\\" + filename[ni]).c_str());
		ver.clear();
		ver_order.clear();
	}
	system("pause");
}
