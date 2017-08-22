#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <climits>
using namespace std;
//-----------------函数的定义------------------
void inputdata();
void clearvector();
//数据输入
double chazhi(vector<double> v1, vector<double>v2, double value);
double O_Z(double o);
double Z_V(double z);
double V_Z(double v);
double H_Nmax(double h);
double func_e(double v1, double v2, double qrk, double r, int i);
void func_Z(double dv, double Zend);
//插值函数
//变量的定义
double zstar = 0;
double zend = 0;
double Qck = 0;
double dt = 24 * 3600;
double qss = 1;
double Z = 0;
double Qfd = 0;
double dh = 1.1;
double P = 0;//出力
double k = 8.5;
double pb = 4990000; //kw
					 //-----------------变量数组的定义------------------
vector<double> Qrk;  //入库流量
vector<double> Zw;  //上游水位
vector<double> Wz;  //库容
vector<double> Zxy;  //上游对应水位
vector<double> Qxy;  //上游下泄量
vector<double> Zmin;  //水位限制下限
vector<double> Zmax;  //水位限制上限
vector<double> hyx;  //预想水头
vector<double> pyx;  //预想出力
vector<double> Qmin;  //生态流量限制
vector<double> Zo;  //上游水位
vector<double> Oz;  //泄流流量限制
vector< vector<double> > V;    // 库容离散点

							   //----------------------------------------------
							   // 计算结果
vector< vector<int> > max_id;      // 最优条件下各指标数值
vector< vector<double> > max_e;
vector< vector<double> > z;
vector< vector<double> > qck;
vector< vector<double> > qfd;
vector< vector<double> > p;
vector< vector<double> > e;

//---------------数据输入------------------------
stringstream stream;     //声明一个stringstream变量
void inputdata(string filename, int year)
{
	string yearstr = filename;
	filename = filename + ".txt";
	cout << string(filename) << endl;
	//读取入库流量
	double temp1 = 0;
	double temp2 = 0;
	ifstream Qrkput(filename, ios::in);
	if (!Qrkput)
	{
		cerr << "读取日流量文件错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (Qrkput.good())
		{
			Qrkput >> temp1;
			Qrk.push_back(temp1);
		}
	}
	cout << "读取日流量数据:" << Qrk.size() << endl;
	Qrkput.close();

	//读取水位库容曲线(亿立方米)
	ifstream ZVinput("ZV.txt", ios::in);
	if (!ZVinput)
	{
		cerr << "读取水位库容文件错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (ZVinput.good())
		{
			ZVinput >> temp1 >> temp2;
			Zw.push_back(temp1);
			Wz.push_back(temp2);
		}
	}
	cout << "读取水位库容数据:" << Wz.size() << endl;
	ZVinput.close();

	//读取泄流下游水位曲线
	ifstream qZinput("qz.txt", ios::in);
	if (!qZinput)
	{
		cerr << "读取泄流下游水位文件错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (qZinput.good())
		{
			qZinput >> temp1 >> temp2;
			Qxy.push_back(temp2);
			Zxy.push_back(temp1);
		}
	}
	cout << "读取泄流下游水位数据:" << Zxy.size() << endl;
	qZinput.close();

	//读取水位限制数据
	if (year % 4 == 0)
	{
		ifstream  Zinput("Zrun.txt", ios::in);  //闰年
		if (!Zinput)
		{
			cerr << "读取水位限制错误" << endl;
			system("pause");
			exit(1);
		}
		else
		{
			while (Zinput.good())
			{
				Zinput >> temp1 >> temp2;
				Zmin.push_back(temp1);
				Zmax.push_back(temp2);

			}
		}
		cout << "读取水位限制数据:" << Zmin.size() << endl;
		Zinput.close();

	}
	else
	{
		ifstream   Zinput("Z.txt", ios::in);   //非闰年
		if (!Zinput)
		{
			cerr << "读取水位限制错误" << endl;
			system("pause");
			exit(1);
		}
		else
		{
			while (Zinput.good())
			{
				Zinput >> temp1 >> temp2;
				Zmin.push_back(temp1);
				Zmax.push_back(temp2);

			}
		}
		cout << "读取水位限制数据:" << Zmin.size() << endl;
		Zinput.close();

	}



	//读取预想出力数据
	string yxname;
	yxname = "yx2.txt";
	//if (year > 2009)
	//{
	//	yxname = "yx2.txt";
	//}
	ifstream Ninput(yxname, ios::in);
	if (!Ninput)
	{
		cerr << "读取预想出力数据错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (Ninput.good())
		{
			Ninput >> temp1 >> temp2;
			hyx.push_back(temp1);
			pyx.push_back(temp2);

		}
	}
	cout << "读取预想出力数据:" << hyx.size() << endl;
	Ninput.close();

	//读取生态流量限制数据
	//if (year % 4 == 0)
	//{
	//	ifstream Qmininput("qmin无约束run.txt", ios::in);  //闰年
	//	if (!Qmininput)
	//	{
	//		cerr << "读取生态流量限制数据错误" << endl;
	//		system("pause");
	//		exit(1);
	//	}
	//	else
	//	{
	//		while (Qmininput.good())
	//		{
	//			Qmininput >> temp1;
	//			Qmin.push_back(temp1);

	//		}
	//	}
	//	cout << "读取生态流量限制数据:" << Qmin.size() << endl;
	//	Qmininput.close();
	//}
	//else
	//{
	//	ifstream Qmininput("qmin无约束.txt", ios::in);     //非闰年
	//	if (!Qmininput)
	//	{
	//		cerr << "读取生态流量限制数据错误" << endl;
	//		system("pause");
	//		exit(1);
	//	}
	//	else
	//	{
	//		while (Qmininput.good())
	//		{
	//			Qmininput >> temp1;
	//			Qmin.push_back(temp1);

	//		}
	//	}
	//	cout << "读取生态流量限制数据:" << Qmin.size() << endl;
	//	Qmininput.close();
	//}
	string name = "Qmin" + yearstr + ".txt";
	/*name = "Qmin2004.txt";*/
	ifstream Qmininput(name, ios::in);
	if (!Qmininput)
	{
		cerr << "读取生态流量限制数据错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (Qmininput.good())
		{
			Qmininput >> temp1;
			Qmin.push_back(temp1);

		}
	}
	cout << "读取生态流量限制数据:" << Qmin.size() << endl;
	Qmininput.close();


	//读取泄流能力限制数据
	ifstream ZOinput("ZOmax.txt", ios::in);
	if (!ZOinput)
	{
		cerr << "读取泄流能力限制数据错误" << endl;
		system("pause");
		exit(1);
	}
	else
	{
		while (ZOinput.good())
		{
			ZOinput >> temp1 >> temp2;
			Zo.push_back(temp1);
			Oz.push_back(temp2);

		}
	}
	cout << "读取泄流能力限制数据:" << Zo.size() << endl;
	ZOinput.close();

	////读取起调水位数据
	//ifstream Ztiao("Ztiao"+ yearstr +".txt", ios::in);
	//if (!Ztiao)
	//{
	//	cerr << "读取起调水位数据错误" << endl;
	//	system("pause");
	//	exit(1);
	//}
	//else
	//{
	//	while (Ztiao.good())
	//	{
	//		Ztiao >> temp1 >> temp2;
	//		zstar=temp1;
	//		zend=temp2;
	//	}
	//}
	//cout << "读取泄流能力限制数据:" << Zo.size() << endl;
	//ZOinput.close();


}
//-----------------------清空函数-------------------
void clearvector()
{
	Qrk.clear();  //入库流量
	Zw.clear();  //上游水位
	Wz.clear();  //库容
	Zxy.clear();  //上游对应水位
	Qxy.clear();  //上游下泄量
	Zmin.clear();  //水位限制下限
	Zmax.clear();  //水位限制上限
	hyx.clear();  //预想水头
	pyx.clear();  //预想出力
	Qmin.clear();  //生态流量限制
	Zo.clear();  //上游水位
	Oz.clear();  //泄流流量限制
	V.clear();    // 库容离散点
	max_id.clear();      // 最优条件下
	max_e.clear();
	z.clear();
	qck.clear();
	qfd.clear();
	p.clear();
	e.clear();
}

//-----------------------插值函数---------------------------------------
double chazhi(vector<double> v1, vector<double>v2, double value)
{
	double k;
	int size = v1.size();

	if (value <= v1[0])
	{
		k = (v2[1] - v2[0]) / (v1[1] - v1[0]);
		if (v2[0] - k * (v1[0] - value) <= 0)
		{
			return  0;
		}
		return v2[0] - k * (v1[0] - value);

	}
	else if (value >= v1[size - 1])
	{
		k = (v2[size - 1] - v2[size - 2]) / (v1[size - 1] - v1[size - 2]);
		if (v2[size - 1] + k * (value - v1[size - 1]) <= 0)
		{
			return  0;
		}
		return v2[size - 1] + k * (value - v1[size - 1]);
	}
	else
	{
		for (int i = 0; i<size - 1; i++)
		{
			if (value == v1[i])
				return v2[i];
			if (value < v1[i + 1])
			{
				k = (v2[i + 1] - v2[i]) / (v1[i + 1] - v1[i]);
				return k * (value - v1[i]) + v2[i];
			}
		}
		return 0;
	}
}

double O_Z(double o)   // 由出库流量得到下游水位
{
	return chazhi(Qxy, Zxy, o);
}

double Z_V(double z)    // 由水位得到库容
{
	return chazhi(Zw, Wz, z);
}

double V_Z(double v)    // 由库容得到水位
{
	return chazhi(Wz, Zw, v);
}

double H_Nmax(double h)    // 由水头得到预想出力
{
	return chazhi(hyx, pyx, h);
}
double Z_Omax(double z)    // 由水位得到泄流能力
{
	return chazhi(Zo, Oz, z);
}

double func_e(double v1, double v2, double qrk, double r, int i, int year,int a ,int b) // 根据初始库容和水库基本参数计算目标函数
{
	double E;
	double zxy, zsy1, zsy2;
	Qck = qrk - (v2 - v1) * 1e8 / dt - qss;

	/*if (Qck < Qmin[i])*/
	if (Qck < 0)
	{
		Qck = -1;
		return  -pow(10, 20);
	}
	else
	{
		//if (Qck < 5000) //下游约束
		if (Qck < Qmin[i])
		{
			if (a == b)
			{
				goto loop;
			}
			return  -pow(10, 20);
		}

	}

loop:

	zxy = O_Z(Qck);
	zsy1 = V_Z(v1);
	zsy2 = V_Z(v2);
	Z = zsy2;
	double qmax;
	qmax = Z_Omax((zsy1 + zsy2) / 2.0);
	Qfd = (qmax<Qck) ? qmax : Qck;

	double h, pyx;
	h = (zsy1 + zsy2) / 2 - zxy;
	h = h - dh;
	if (h < 0)
	{
		return  -pow(10, 20);
	}
	pyx = H_Nmax(h) * 1000;
	P = k * Qfd * h;

	int n;
	n = int(P / pyx);
	if (n > 32)
	{
		P = 32 * pyx;
		Qfd = P / (k*h);
	}
		
	E = P * dt / 3600;

    if(P>pb)
		return E;
	else
		return (P + (P - pb)*r) * dt / 3600 ;
}
//-----------------------------------------------------------
void func_Z(double dv = 1, double Zstar = 170, double Zend = 170)   // 根据水位上限、下限和初末水位计算每个时段的库容离散点(二维)
{
	double Vmax, Vmin, v,vend;
	int i = 0;
	double v_length = Zmax.size(); //整个序列长度
	V.resize(v_length + 1);

	v = Z_V(Zstar);
	V[0].push_back(v);

	for (i = 1; i < v_length; i++)
	{
		Vmax = Z_V(Zmax[i]);
		Vmin = Z_V(Zmin[i]);
		while (Vmax > Vmin)
		{
			V[i].push_back(Vmax);
			Vmax -= dv;
		}
		V[i].push_back(Vmin);
	}
	vend = Z_V(Zend);
	V[i].push_back(vend);
	//输出检测
	ofstream outfile("库容离散输出.txt");
	for (int i = 0; i < V.size(); i++)
	{
		outfile << "第" << i + 1 << " 阶段： ";
		for (int j = 0; j < V[i].size(); j++)
		{
			outfile << V[i][j] << " ";
		}
		outfile << endl;
	}
	cout << "离散完成" << endl;
}
void main()
{
	string classset = "173.5HY-";
	ofstream zongjie(classset + " Result.txt");
	ofstream IHA(classset + " IHA.txt");
	
	IHA << "年份 " << "阶段 " << " 出流" << " 入流" << " 水位" <<" 生态"<<" 出力"<< endl;
	zongjie << "年份 " << "发电量 " << "发电保证率 " <<" 生态保证率"<<" 缺水量"<< endl;
	//---------------------------------------
	for (int year =1981; year <2013; year++)
	//---------------------------------------
	{
		clearvector();
		string filename;
		char t[256];
		sprintf(t, "%d", year);
		filename = t;
		//stream << year;
		//stream >> filename;
		//---------------------------------------
		inputdata(filename, year);
		//------------------------------------
		if (year == 1992 || year == 2007 || year == 2009||year==1997||year==2002)
		{
			func_Z(1, 170, 170);//初始化完成（1：222）
		}
		else
		{
			func_Z(1, 173.5, 173.5);//初始化完成（1：222）
		}
		//--------------------------------------

	   //------------二维数组------------------
		double r = 500;//罚函数
		double E = 0;
		double length = Qrk.size();
		max_id.resize(length + 1);
		max_e.resize(length + 1);
		e.resize(length + 1);
		z.resize(length);
		qck.resize(length);
		qfd.resize(length);
		p.resize(length);

		for (int i = 0; i < length + 1; i++)
		{
			max_id[i].resize(V[i].size());
			max_e[i].resize(V[i].size());
			e[i].resize(V[i].size());
		}
		max_e[length][0] = 0;
		e[length][0] = 1;
		for (int i = 0; i < length; i++)
		{
			z[i].resize(V[i].size());
			qck[i].resize(V[i].size());
			qfd[i].resize(V[i].size());
			p[i].resize(V[i].size());
		}
		//---------------顺序递推----------------
		ofstream outf("输出TEST"+filename +".txt");
		for (int i = length - 1; i >= 0; i--)
		{
			cout << "-----------第" << i << " 阶段------------- " << endl;
			outf << "-----------第" << i << " 阶段------------- " << endl;
			for (int j = 0; j < V[i].size(); j++)
			{
				z[i][j] = 10;
				max_e[i][j] = 0;
				e[i][j] = 0;
				for (int kk = 0; kk < V[i + 1].size(); kk++)
				{
					E = func_e(V[i][j], V[i + 1][kk], Qrk[i], r, i, year,kk, V[i + 1].size()-1);
					E = E + max_e[i + 1][kk];
					if (E >= max_e[i][j] || max_e[i][j] == 0)
					{

						max_e[i][j] = E;
						max_id[i][j] = kk;
						z[i][j] = Z;
						e[i][j] = E;
						qck[i][j] = Qck;
						qfd[i][j] = Qfd;
						p[i][j] = P / 1000; //MW
					}
					/*outf <<j<<" "<< V_Z(V[i][j])<<" "<<z[i][j]<<" " <<max_id[i][j]<<Qrk[i]<<" "<<qck[i][j]<<" "<<qfd[i][j]<<" "<<p[i][j]<<"  "<<e[i][j]<<endl;*/
				}
				outf << j << " " << V_Z(V[i][j]) << " " << z[i][j] << " " << max_id[i][j] <<" "<< Qrk[i] << " " << qck[i][j] << " " << qfd[i][j] << " " << p[i][j] << "  " << e[i][j] << endl;
			}
		}
		cout << "-------逆序计算完成------" << endl;
		outf.close();
		//-------------顺序回代(最优序列)--------------
		ofstream output(classset + " 流量限制" + filename + ".txt");
		//output << "阶段 " << "水位(m)" << " " << "最佳ID" << " " << "入库流量(m3/s)" << " 出库流量(m3/s)" << " " << "发电流量(m3/s)" << " " <<"出力(MW)" << "  " <<"指标函数"<< endl;
		output << "阶段 " << "水位(m)" << " " << "入库流量(m3/s)" << " 出库流量(m3/s)" << " " << "发电流量(m3/s)" << " " << "出力(MW)"  << endl;
		double sume = 0;//总发电量
		double sumque = 0;//总缺水
		int id = 0;
		int count = 0;
		int count2 = 0;
		double bzl = 0;
		double bzl2 = 0;
		double queshui = 0;
		for (int i = 0; i < length; i++)
		{
			if (p[i][id] > pb / 1000)
			{
				count++;
			}
			if (qck[i][id]>Qmin[i])
			{
				count2++;
				queshui = 0;
			}
			else
			{
				queshui = (qck[i][id] - Qmin[i])*24*3600/1e8;
			}

			//output<<"第"<<i+1<<"阶段 "<< z[i][id] << " " << max_id[i][id] << " "<<Qrk[i]<<" " << qck[i][id] << " " << qfd[i][id] << " " << p[i][id] << "  " << e[i][id] << endl;
			output << "第" << i + 1 << "阶段 " << z[i][id] << " " << Qrk[i] << " " << qck[i][id] << " " << qfd[i][id] << " " << p[i][id] <<" "<< queshui << endl;
			IHA << year << " " << i + 1 << " " << qck[i][id] << " " << Qrk[i]  << " " << z[i][id] <<" "<<Qmin[i]<<" " << p[i][id] << endl;
			sume += p[i][id] * 24;//MW
			sumque += queshui;//亿m³
			id = max_id[i][id];
		}
		sume = sume * 1000 / 1e8;
		bzl = count*1.0 / length;
		bzl2 = count2*1.0 / length;
		cout << filename << "发电量 " << sume << " 亿kW·h 发电保证率 " << bzl * 100 <<  endl;
		output << filename << "发电量 " << sume << " 亿kW·h 发电保证率 " << bzl * 100 <<  endl;
		output << "罚系数r " << r << endl;
		//output <<filename << "发电量 " << sume << " 亿kW·h 保证率 " << bzl * 100 << "%" << endl;
	/*	zongjie << "起始 " << zstar << " 末水位 " << zend << endl;*/
		zongjie << year << " " << sume << " " << bzl * 100 <<" "<< bzl2 *100<<" "<< sumque <<  endl;
		//-------------------------------------------
		output.close();
	}
	zongjie.close();
	IHA.close();
	system("pause");
}
