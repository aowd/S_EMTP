#include "SynchGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
SynchGenerator::SynchGenerator(int a,int b,int c)
{
	isUserDef = 1;
	need_NEC=1;
	//构造函数，a、b、c为节点编号
	nPort=3;
	nodeNumber=new int[3];
	type=6;

	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	w = 2*PI*60;//电角速度，有名值,f=60Hz
	Ef = 1.3255;//励磁电压，标幺值
	//发电机的有功和无功，有名值
	P0 = 126391.44e6;//MW
	Q0 = 238242.25e6;//MVar
	Sm = 0;//视在功率，标幺值，可计算得到,先赋为0，初始化时再重新赋值
	Vm = 1.0;//相电压幅值，标幺值
	Im = 0;//相电流幅值，标幺值，可计算得到,先赋为0，初始化时再重新赋值
	ph = 1.642;//A相电压的初始相位，弧度表示
	Sangle = 0;//功率因素角，可计算得到,先赋为0，初始化时再重新赋值
	Vangle = 0;//A相电压初始相角，可计算得到,先赋为0，初始化时再重新赋值
	Iangle = 0;//A相电流初始相角，可计算得到,先赋为0，初始化时再重新赋值
	Angle = 0;//转子角的初始值，可计算得到，这个变量没有用上

	//dq0坐标系下的电阻参数,标幺值
	Ra = 0.00621;
	Rf = 0.00025765;
	RD = 0.003033;
	RQ = 0.0030469;
	//dq0坐标系下的电抗参数,标幺值
	Xad = 0.34288;
	Xaq = 0.15453;
	Xd = 0.40075;
	Xq = 0.2124;
	Xf = 0.37918;
	XD = 0.35411;
	XQ = 0.162953;
	XfD = Xad;
	Xaf = Xad;
	XaD = Xad;
	XaQ = Xaq;

	//定子侧电压、电流和阻抗的基值
	Vbase = 277e3*sqrt(2.0);//V
	Ibase = 366e3*sqrt(2.0);//A
	Zbase = Vbase/Ibase;

	//机械侧参数
	H = 1.0;
	D = 0;
	Tm = 0.4221;//输入转矩

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/
	//计算coff和T
	coff = 99.0/101.0;
	T = 1.0/deltaT/w*0.5*(1+coff);

	//励磁电压的系数
	Vf=Rf/Xad;

	//派克变换矩阵及其逆矩阵
	P.ResizeTo(1,2,1,3);
	invP.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,5,1,5);
	BB.ResizeTo(1,5,1,5);
	////
	Rss.ResizeTo(1,2,1,2);
	Rsr.ResizeTo(1,2,1,3);
	Rrs.ResizeTo(1,3,1,2);
	Rrr.ResizeTo(1,3,1,3);

	////
	Bss.ResizeTo(1,2,1,2);
	Bsr.ResizeTo(1,2,1,3);
	Brs.ResizeTo(1,3,1,2);
	Brr.ResizeTo(1,3,1,3);

	////处理系数矩阵，使得转换到abc坐标系时能解耦
	Rs.ResizeTo(1,2,1,2);
	R_ave.ResizeTo(1,2,1,2);
	R_res.ResizeTo(1,2,1,2);
	inv_Rs.ResizeTo(1,2,1,2);
	inv_Rr.ResizeTo(1,3,1,3);

	////
	Gs1.ResizeTo(1,2,1,2);
	Gs2.ResizeTo(1,2,1,3);
	Gs3.ResizeTo(1,2,1,2);
	Gs4.ResizeTo(1,2,1,3);
	Gs5.ResizeTo(1,2,1,3);
	Gs6.ResizeTo(1,2,1,2);
	Gs13.ResizeTo(1,2,1,2);
	Gs25.ResizeTo(1,2,1,3);

	////
	Gr1.ResizeTo(1,3,1,3);
	Gr2.ResizeTo(1,3,1,2);
	Gr3.ResizeTo(1,3,1,2);
	Gr4.ResizeTo(1,3,1,3);
	Gr5.ResizeTo(1,3,1,3);
	Gr15.ResizeTo(1,3,1,3);

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值，先令其等于0，初始化时将赋初值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vabc[i-1] = 0;
		Vabc_1[i-1] = 0;
		Vabc_2[i-1] = 0;
		Iabc[i-1] = 0;
		Iabc_1[i-1] = 0;
		Iabc_2[i-1] = 0;
	}
	//dq0坐标系下定子线圈的电压电流
	for (int i=1;i<=2;i++)
	{
		Vdq0[i-1] = 0;
		Idq0[i-1] = 0;
		Idq0_his[i-1] = 0;
		Idq0_his2[i-1] = 0;
		Idq0_forcast[i-1] = 0;
	}
	//dq0坐标系下转子侧的电压电流
	for (int i=1;i<=3;i++)
	{
		VfDQ[i-1] = 0;
		IfDQ[i-1] = 0;
		IfDQ_his[i-1] = 0; 
	}

	/*  机械侧物理量  */
	wr = 1;//转子转速
	wr_his = 1;//wr的历史值
	curAngle = 0;//转子角,先赋为0，初始化时再重新赋值
	curAngle_his = 0;//curAngle的历史值,没用上

	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//等效导纳的历史值，用于更新节点导纳矩阵

	for(int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0坐标下的诺顿等效电流，标幺值
	}

	/************************************************************************/
	/*              与EMTP接口所需的变量                                    */
	/************************************************************************/
	nodeNumber[0] = a;
	nodeNumber[1] = b;
	nodeNumber[2] = c;

	/************************************************************************/
	/*              计算系数矩阵和诺等等值导纳                              */
	/************************************************************************/
	calculateCoefficientMatrix();
}
SynchGenerator::~SynchGenerator()
{//析构函数
	
}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/
void SynchGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time)
{//初始化支路电压电流
	//初始化wr和curAngle
	Sm = sqrt(P0*P0+Q0*Q0);
	Sangle = atan(Q0/P0);
	Sm = Sm/(1.5*Vbase*Ibase);
	Im = Sm/Vm;
	Vangle = ph/180*PI;
	Iangle = Vangle - Sangle;
	curAngle = atan((Vm*sin(Vangle)+Xq*Im*cos(Iangle)+Ra*Im*sin(Iangle))/
		(Vm*cos(Vangle)-Xq*Im*sin(Iangle)+Ra*Im*cos(Iangle)))
		-0.5*PI;

	Iabc[0] = Im*cos(Iangle);
	Iabc[1] = Im*cos(Iangle-(2.0/3.0)*PI);
	Iabc[2] = Im*cos(Iangle+(2.0/3.0)*PI);

	Vabc[0] = Vm*cos(Vangle);
	Vabc[1] = Vm*cos(Vangle-(2.0/3.0)*PI);
	Vabc[2] = Vm*cos(Vangle+(2.0/3.0)*PI);

	parkTransMatrix(P,curAngle);

	for (int i=0;i<2;i++)
	{
		Vdq0[i] = P(i+1,1)*Vabc[0]+P(i+1,2)*Vabc[1]+P(i+1,3)*Vabc[2];
		Vdq0_his[i] = Vdq0[i];
		Vdq0_his2[i] = Vdq0[i];
		Vdq0_emend[i] = Vdq0[i];
		Idq0[i] = P(i+1,1)*Iabc[0]+P(i+1,2)*Iabc[1]+P(i+1,3)*Iabc[2];
		Idq0_his[i] = Idq0[i];
		Idq0_his2[i] = Idq0[i];
		Idq0_forcast[i] = Idq0[i];
		Idq0_emend[i] = Idq0[i];
	}

	IfDQ[0] = (Vdq0[1]+Ra*Idq0[1]+Xd*Idq0[0])/Xad;
	IfDQ[1] = 0;
	IfDQ[2] = 0;

	VfDQ[0] = IfDQ[0]*Rf;
	VfDQ[1] = 0;
	VfDQ[2] = 0;

	//解机械方程，初始化wr和curAngle
	calculateOmega();

	//诺顿等值电流初始化
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			- Gs1(i,1)*Idq0_forcast[0]-Gs1(i,2)*Idq0_forcast[1]
		- Gs3(i,1)*Idq0[0]-Gs3(i,2)*Idq0[1]
		+ Gs4(i,1)*IfDQ[0]+Gs4(i,2)*IfDQ[1]+Gs4(i,3)*IfDQ[2]
		+ Gs6(i,1)*Vdq0[0]+Gs6(i,2)*Vdq0[1]
		+ Gs25(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invP,curAngle);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invP(i,1)*nortonEquivalentCurrent_dq0[0]
		+ invP(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}

	//用PSCAD数据初始化定子电压电流
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<3;i++)
	{
		Iabc[i] = initialCurrentArray[ptr+i];
	}
	ptr+=3;

	calculateCoefficientMatrix();
}
void SynchGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//从节点电压数组中读支路节点电压

	for (int i=0;i<3;i++)
	{
		Vabc_1[i] =Vabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vabc[i]=0;
		}
		else
		{
			Vabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}
}
void SynchGenerator::calculateBranchVoltage()
{//计算支路电压
}

void SynchGenerator::calculateBranchCurrent()
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Iabc_1[i] = Iabc[i];
		Iabc[i] = -Vabc[i]*nortonEquivalentConductance-nortonEquivalentCurrent[i];
	}
}

void SynchGenerator::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项
	calculateCoefficientMatrix();
	calculateDq0Results();//计算dq0坐标系下的结果，预报定子电流
	calculateOmega();
	
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			- Gs1(i,1)*Idq0_forcast[0]-Gs1(i,2)*Idq0_forcast[1]
			- Gs3(i,1)*Idq0[0]-Gs3(i,2)*Idq0[1]
			+ Gs4(i,1)*IfDQ[0]+Gs4(i,2)*IfDQ[1]+Gs4(i,3)*IfDQ[2]
			+ Gs6(i,1)*Vdq0[0]+Gs6(i,2)*Vdq0[1]
			+ Gs25(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invP,curAngle);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invP(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invP(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}
}

void SynchGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//形成节点诺顿等效电流向量
	int N;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}
}
void SynchGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//形成节点导纳阵
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance;
	}
}

void SynchGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix(counter,ptr) = Iabc[0];
	branchCurrentMatrix(counter,ptr+1) = Iabc[1];
	branchCurrentMatrix(counter,ptr+2) = Iabc[2];
	ptr+=3;
}

void SynchGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix_1[counter-1][ptr] = Iabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Iabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Iabc[2];
	ptr+=3;
}

void SynchGenerator::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值

	//电压电流变量插值
	for (int k=0;k<3;k++)
	{
		Vabc[k] = (1-ratio)*Vabc_1[k] + ratio*Vabc[k];
		Iabc[k] = (1-ratio)*Iabc_1[k] + ratio*Iabc[k];
	}

	//诺顿等值电流插值
	double temp1,temp2;
	for (int k=0;k<3;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void SynchGenerator::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			Vabc_2[i] = Vabc_1[i];
			Iabc_2[i] = Iabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<3;i++)
		{
			Vabc_1[i] = Vabc_2[i];
			Iabc_1[i] = Iabc_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

void SynchGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
{//转子转速变化时更新导纳矩阵
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his+this->nortonEquivalentConductance;
	}
}

/************************************************************************/
/*                发电机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
//同时重算nortonEquivalentConductance
void SynchGenerator::calculateCoefficientMatrix()
{
	//计算AA矩阵，即伴随矩阵
	AA(1,1) = Ra+2.0*T*Xd;
	AA(1,2) = -wr*Xq;
	AA(1,3) = 2.0*T*Xaf;
	AA(1,4) = 2.0*T*XaD;
	AA(1,5) = -wr*XaQ;

	AA(2,1) = wr*Xd;
	AA(2,2) = Ra+2.0*T*Xq;
	AA(2,3) = wr*Xaf;
	AA(2,4) = wr*XaD;
	AA(2,5) = 2.0*T*XaQ;

	AA(3,1) = 2.0*T*Xaf;
	AA(3,2) = 0;
	AA(3,3) = Rf+2.0*T*Xf;
	AA(3,4) = 2.0*T*XfD;
	AA(3,5) = 0;

	AA(4,1) = 2.0*T*XaD;
	AA(4,2) = 0;
	AA(4,3) = 2.0*T*XfD;
	AA(4,4) = RD+2.0*T*XD;
	AA(4,5) = 0;

	AA(5,1) = 0;
	AA(5,2) = 2.0*T*XaQ;
	AA(5,3) = 0;
	AA(5,4) = 0;
	AA(5,5) = RQ+2.0*T*XQ;

	//计算BB矩阵
	BB(1,1) = coff*Ra-2.0*T*Xd;
	BB(1,2) = -coff*wr*Xq;
	BB(1,3) = -2.0*T*Xaf;
	BB(1,4) = -2.0*T*XaD;
	BB(1,5) = -coff*wr*XaQ;

	BB(2,1) = coff*wr*Xd;
	BB(2,2) = coff*Ra-2.0*T*Xq;
	BB(2,3) = coff*wr*Xaf;
	BB(2,4) = coff*wr*XaD;
	BB(2,5) = -2.0*T*XaQ;

	BB(3,1) = -2.0*T*Xaf;
	BB(3,2) = 0;
	BB(3,3) = coff*Rf-2.0*T*Xf;
	BB(3,4) = -2.0*T*XfD;
	BB(3,5) = 0;

	BB(4,1) = -2.0*T*XaD;
	BB(4,2) = 0;
	BB(4,3) = -2.0*T*XfD;
	BB(4,4) = coff*RD-2.0*T*XD;
	BB(4,5) = 0;

	BB(5,1) = 0;
	BB(5,2) = -2.0*T*XaQ;
	BB(5,3) = 0;
	BB(5,4) = 0;
	BB(5,5) = coff*RQ-2.0*T*XQ;

	//计算矩阵Rss,Rsr,Rrs,Rrr,Bss,Bsr,Brs,Brr
	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rsr(i,j) = AA(i,j+2);
			Bsr(i,j) = BB(i,j+2);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrs(i,j) = AA(i+2,j);
			Brs(i,j) = BB(i+2,j);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rrr(i,j) = AA(i+2,j+2);
			Brr(i,j) = BB(i+2,j+2);
		}
	}

	//计算Rs,R_ave,R_res,inv_Rs,inv_Rr
	inv_Rr = Rrr;
	inv_Rr.Invert();

	Rs = Rss-Rsr*inv_Rr*Rrs;

	double resistance = (Rs(1,1)+Rs(2,2))/2.0;
	R_ave(1,1) = resistance;
	R_ave(1,2) = 0;
	R_ave(2,1) = 0;
	R_ave(2,2) = resistance;
	
	R_res = Rs - R_ave;

	//计算诺顿等效电导
	nortonEquivalentConductance = 1/resistance;
	inv_Rs = R_ave;
	inv_Rs.Invert();
	nortonEquivalentConductance /= Zbase;

	//计算Gs1...Gs6
	Gs1 = inv_Rs*R_res;
	Gs1 = (-1.0)*Gs1;
	Gs2 = inv_Rs*Rsr*inv_Rr;
	Gs2 = (-1.0)*Gs2;
	Gs3 = Gs2*Brs + inv_Rs*Bss;
	Gs3 = (-1.0)*Gs3;
	Gs4 = Gs2*Brr + inv_Rs*Bsr;
	Gs4 = (-1.0)*Gs4;
	Gs5 = coff*Gs2;
	Gs6 = coff*inv_Rs;
	Gs13 = Gs1+Gs3;
	Gs25 = Gs2+Gs5;

	//计算Gr1...Gr5
	Gr1 = inv_Rr;
	Gr2 = inv_Rr*Rrs;
	Gr2 = (-1.0)*Gr2;
	Gr3 = inv_Rr*Brs;
	Gr3 = (-1.0)*Gr3;
	Gr4 = inv_Rr*Brr;
	Gr4 = (-1.0)*Gr4;
	Gr5 = coff*inv_Rr;
	Gr15 = Gr1+Gr5;
}

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
//为下一时步计算诺顿等效电流做准备
void SynchGenerator::calculateDq0Results()
{
	double Vabc_pu[3];
	double Iabc_pu[3];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		Vabc_pu[i] = Vabc[i]/Vbase;
		Iabc_pu[i] = Iabc[i]/Ibase;
	}

	//计算park变换矩阵
	parkTransMatrix(P,curAngle);

	//计算定子电压
	for (int i=0;i<2;i++)
	{
		Vdq0_his2[i] = Vdq0_his[i];
		Vdq0_his[i] = Vdq0[i];
		Vdq0[i] = P(i+1,1)*Vabc_pu[0]+P(i+1,2)*Vabc_pu[1]+P(i+1,3)*Vabc_pu[2];
	}

	//计算定子电流
	for (int i=0;i<2;i++)
	{
		Idq0_his2[i] = Idq0_his[i];
		Idq0_his[i] = Idq0[i];
		Idq0[i] = P(i+1,1)*Iabc_pu[0]+P(i+1,2)*Iabc_pu[1]+P(i+1,3)*Iabc_pu[2];
	}
	
	//预报下一时刻定子电流
	for (int i=0;i<2;i++)
	{
		Idq0_forcast[i] = 1.3333333*Idq0[i]+0.3333333*Idq0_his[i]-0.6666666*Idq0_his2[i];
	}

	//计算转子电压
	Vf = Rf/Xad;
	VfDQ[0] = Ef*Vf;
	VfDQ[1] = 0;
	VfDQ[2] = 0;

	//计算转子电流
	for (int i=0;i<3;i++)
	{
		IfDQ_his[i] = IfDQ[i];
	}
	for (int i=1;i<=3;i++)
	{
		IfDQ[i-1] =- Gr2(i,1)*Idq0[0] - Gr2(i,2)*Idq0[1]
		- Gr3(i,1)*Idq0_his[0] - Gr3(i,2)*Idq0_his[1]
		+ Gr4(i,1)*IfDQ_his[0] + Gr4(i,2)*IfDQ_his[1] + Gr4(i,3)*IfDQ_his[2]
		+ Gr15(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}
}

//求解机械方程，计算转子转速
//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
void SynchGenerator::calculateOmega()
{
	double tt1,tt2,tt3;
	tt1 = (Xd-Xq)*Idq0[0]*Idq0[1];
	tt2 = -Xad*(IfDQ[0]*Idq0[1]+IfDQ[1]*Idq0[1]);
	tt3 = -Xaq*IfDQ[2]*Idq0[0];

	wr_his = wr;
	wr = (Tm+(tt1+tt2-tt3))*deltaT/2/H+wr_his;

	curAngle = curAngle + deltaT*wr*w/2 + deltaT*wr_his*w/2;
	if (curAngle>2*PI)
	{
		curAngle -= 2*PI;
	}
}

//计算派克变换与反变换矩阵
void SynchGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void SynchGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}
