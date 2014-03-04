#include "WoundInducGenerator2.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                构造函数与析构函数                                    */
/************************************************************************/
WoundInducGenerator2::WoundInducGenerator2(int id, int firstNode,int lastNode, int control,double Vw)
{
	this->id = id;
	isUserDef = 0;
	need_NEC=1;
	//构造函数，a、b、c为节点编号
	nPort=6;
	nodeNumber=new int[nPort];
	type=29;

	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	Wbase = 2*PI*50;//电角速度，有名值,f=60Hz
	w=1;

	//dq0坐标系下的电阻参数,标幺值
	Rs = 0.048;
	Rr = 0.018;
	//dq0坐标系下的电抗参数,标幺值
	Xm = 3.8;
	Xls = 0.075;
	Xlr = 0.12;

	//定子侧电压、电流和阻抗的基值
	Vbase=690;
	Sbase=2e6;
	Vbase = Vbase/sqrt(3.0);//V
	Ibase = Sbase/3/Vbase;//A
	Zbase = Vbase/Ibase;

	//定转子绕组匝数比
	turnratio=0.3;

	this->control = control;

	//机械侧参数
	H = 0.1;
	D = 0;
	TL =-0.7;//输入转矩

	//风机参数
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=Vw-0.1*(id-1);

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/

	//派克变换矩阵及其逆矩阵
	Ps.ResizeTo(1,3,1,3);
	invPs.ResizeTo(1,3,1,3);
	Pr.ResizeTo(1,3,1,3);
	invPr.ResizeTo(1,3,1,3);
	P.ResizeTo(1,6,1,6);
	invP.ResizeTo(1,6,1,6);

	//
	AA.ResizeTo(1,6,1,6);
	BB.ResizeTo(1,6,1,6);
	AA_inv.ResizeTo(1,6,1,6);
	//
	nortonEquivalentConductanceMatrix.ResizeTo(1,6,1,6);

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值，先令其等于0，初始化时将赋初值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vsabc[i-1] = 0;
		Vsabc_1[i-1] = 0;
		Vsabc_2[i-1] = 0;
		Isabc[i-1] = 0;
		Isabc_1[i-1] = 0;
		Isabc_2[i-1] = 0;
	}
	//dq0坐标系下定子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vsdq0[i-1] = 0;
		Isdq0[i-1] = 0;
	}

	//abc坐标系下转子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vrabc[i-1] = 0;
		Vrabc_1[i-1] = 0;
		Vrabc_2[i-1] = 0;
		Irabc[i-1] = 0;
		Irabc_1[i-1] = 0;
		Irabc_2[i-1] = 0;
	}
	//dq0坐标系下转子线圈的电压电流
	for (int i=1;i<=3;i++)
	{
		Vrdq0[i-1] = 0;
		Irdq0[i-1] = 0;
	}

	/*  机械侧物理量  */
	if (control==1)
	{
		wr = 1.15;//转子转速
		wr_his = 1.15;//wr的历史值
	} 
	else
	{
		wr = 0;//转子转速
		wr_his = 0;//wr的历史值
	}
	curAngle_s= 0;//转子角,先赋为0，初始化时再重新赋值
	curAngle_r= 0;//转子角,先赋为0，初始化时再重新赋值
	Te=0;
	Te_his=0;


	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	for (int i=0;i<36;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_pu[i]=0;
	}

	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0坐标下的诺顿等效电流，标幺值
		nortonEquivalentCurrent_dq0_his[i-1] = 0;
	}

	/************************************************************************/
	/*              与EMTP接口所需的变量                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}
#ifdef WIND_VELOCITY_DATA_INPUT
	WindVelocityCounter=0;
	WindVelocityCounter_bak=0;
#endif
}
WoundInducGenerator2::~WoundInducGenerator2()
{//析构函数

}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/

void WoundInducGenerator2::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time){ptr+=6;}

void WoundInducGenerator2::readNodeVoltage(TVectorD &nodeVoltageArray)
{//从节点电压数组中读支路节点电压

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
		Vsabc_his2[i] = Vsabc_his[i];
		Vsabc_his[i] = Vsabc[i];
	}

	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vsabc[i]=0;
		}
		else
		{
			Vsabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}

	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
		Vrabc_his2[i] = Vrabc_his[i];
		Vrabc_his[i] = Vrabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i+3]==0)
		{
			Vrabc[i]=0;
		}
		else
		{
			Vrabc[i]= nodeVoltageArray(nodeNumber[i+3]);
		}	
	}
}
void WoundInducGenerator2::calculateBranchVoltage()
{//计算支路电压
}

void WoundInducGenerator2::calculateBranchCurrent()
{//计算支路电流
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i]=0;
		for (int j=0;j<3;j++)
		{
			Isabc[i]+=Vsabc[j]*nortonEquivalentConductance[6*i+j];
			Isabc[i]+=Vrabc[j]*nortonEquivalentConductance[6*i+j+3];
		}
		Isabc[i]+=nortonEquivalentCurrent[i];

		Irabc_1[i] = Irabc[i];
		Irabc[i]=0;
		for (int j=0;j<3;j++)
		{
			Irabc[i]+=Vsabc[j]*nortonEquivalentConductance[6*3+6*i+j];
			Irabc[i]+=Vrabc[j]*nortonEquivalentConductance[6*3+6*i+j+3];
		}
		Irabc[i]+=nortonEquivalentCurrent[i+3];
	}

	calculateDq0Results();
}

void WoundInducGenerator2::calculateNortonEquivalentCurrent(double time)
{//计算支路的诺顿等效电路中的电流项

	//calculateDq0Results();//计算dq0坐标系下的结果
	//calculateOmega();
	//calculateCoefficientMatrix();

	TMatrixD CC=AA_inv*BB;
	//计算dq0坐标下的历史电流项
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];

		nortonEquivalentCurrent_dq0[i-1] =-(CC(i,1)*Isdq0[0]+CC(i,2)*Isdq0[1]+CC(i,3)*Isdq0[2]+CC(i,4)*Irdq0[0]+CC(i,5)*Irdq0[1]+CC(i,6)*Irdq0[2])
			+ AA_inv(i,1)*Vsdq0[0]+AA_inv(i,2)*Vsdq0[1]+AA_inv(i,3)*Vsdq0[2]+AA_inv(i,4)*Vrdq0[0]+AA_inv(i,5)*Vrdq0[1]+AA_inv(i,6)*Vrdq0[2];


		nortonEquivalentCurrent_dq0_his[i+2] = nortonEquivalentCurrent_dq0[i+2];

		nortonEquivalentCurrent_dq0[i+2] =-(CC(i+3,1)*Isdq0[0]+CC(i+3,2)*Isdq0[1]+CC(i+3,3)*Isdq0[2]+CC(i+3,4)*Irdq0[0]+CC(i+3,5)*Irdq0[1]+CC(i+3,6)*Irdq0[2])
			+ AA_inv(i+3,1)*Vsdq0[0]+AA_inv(i+3,2)*Vsdq0[1]+AA_inv(i+3,3)*Vsdq0[2]+AA_inv(i+3,4)*Vrdq0[0]+AA_inv(i+3,5)*Vrdq0[1]+AA_inv(i+3,6)*Vrdq0[2];
	}

	//计算park逆变换矩阵
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//计算abc坐标下的历史电流项
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]+ invPs(i,2)*nortonEquivalentCurrent_dq0[1]+invPs(i,3)*nortonEquivalentCurrent_dq0[2];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[3]+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[4]+invPr(i-3,3)*nortonEquivalentCurrent_dq0[5];			
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}		
	}
}

void WoundInducGenerator2::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//形成节点诺顿等效电流向量
	int N;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) = nodeNortonEquivalentCurrentArray(N)- nortonEquivalentCurrent[i];
	}
}
void WoundInducGenerator2::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//形成节点导纳阵
	int N,M;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumber[i];
		for (int j=0;j<nPort;j++)
		{
			M=nodeNumber[j];
			conductanceMatrix(N,M)+=this->nortonEquivalentConductance[nPort*i+j];
		}
	}
}

void WoundInducGenerator2::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix(counter,ptr) = Isabc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2];
	ptr+=3;
	branchCurrentMatrix(counter,ptr) = Irabc[0];
	branchCurrentMatrix(counter,ptr+1) = Irabc[1];
	branchCurrentMatrix(counter,ptr+2) = Irabc[2];
	ptr+=3;
}

void WoundInducGenerator2::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//保存支路电流
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2];
	ptr+=3;
	branchCurrentMatrix_1[counter-1][ptr] = Irabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Irabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Irabc[2];
	ptr+=3;
}
void WoundInducGenerator2::interpolate(double ratio)
{//给定比值，对支路的电压电流进行插值

	//电压电流变量插值
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
	}

	for (int k=0;k<3;k++)
	{
		Vsdq0[k] = (1-ratio)*Vsdq0_his[k] + ratio*Vsdq0[k];
		Isdq0[k] = (1-ratio)*Isdq0_his[k] + ratio*Isdq0[k];
		Vrdq0[k] = (1-ratio)*Vrdq0_his[k] + ratio*Vrdq0[k];
		Irdq0[k] = (1-ratio)*Irdq0_his[k] + ratio*Irdq0[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_dq0[k] = (1-ratio)*nortonEquivalentCurrent_dq0_his[k] + ratio*nortonEquivalentCurrent_dq0[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent[k] = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	}

	curAngle_r = (1-ratio)*curAngle_r_1 + ratio*curAngle_r;
	curAngle_s = (1-ratio)*curAngle_s_1 + ratio*curAngle_s;
	wr = (1-ratio)*wr_his + ratio*wr;

	//诺顿等值电流插值
	//double temp1,temp2;
	//for (int k=0;k<6;k++)
	//{
	//	temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	//	temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

	//	nortonEquivalentCurrent_1[k] = temp1;
	//	nortonEquivalentCurrent[k] = temp2;
	//}
}

void WoundInducGenerator2::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_2[i] = Vsabc_1[i];
			Isabc_2[i] = Isabc_1[i];
			Vrabc_2[i] = Vrabc_1[i];
			Irabc_2[i] = Irabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_2[i+3] = nortonEquivalentCurrent_1[i+3];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			Vrabc_1[i] = Vrabc_2[i];
			Irabc_1[i] = Irabc_2[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
			nortonEquivalentCurrent_1[i+3] = nortonEquivalentCurrent_2[i+3];
		}
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

/************************************************************************/
/*                发电机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB
//同时重算nortonEquivalentConductance
void WoundInducGenerator2::calculateCoefficientMatrix()
{
	//计算AA矩阵，即伴随矩阵
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w*(Xm+Xls);
	AA(1,4) = T*Xm;
	AA(1,5) = -w*Xm;

	AA(2,1) = w*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,4) = w*Xm;
	AA(2,5) = T*Xm;

	AA(3,3)=Rs+T*(Xm+Xls);
	AA(3,6)=T*Xm;

	AA(4,1) = T*Xm;
	AA(4,2) = (wr-w)*Xm;
	AA(4,4) = Rr+T*(Xm+Xlr);
	AA(4,5) = (wr-w)*(Xm+Xlr);

	AA(5,1) = (w-wr)*Xm;
	AA(5,2) = T*Xm;
	AA(5,4) = (w-wr)*(Xm+Xlr);
	AA(5,5) = Rr+T*(Xm+Xlr);

	AA(6,3)=T*Xm;
	AA(6,6)=Rr+T*(Xm+Xlr);

	//计算BB矩阵
	BB(1,1) = Rs-T*(Xm+Xls);
	BB(1,2) = -w*(Xm+Xls);
	BB(1,4) = -T*Xm;
	BB(1,5) = -w*Xm;

	BB(2,1) = w*(Xm+Xls);
	BB(2,2) = Rs-T*(Xm+Xls);
	BB(2,4) = w*Xm;
	BB(2,5) = -T*Xm;

	BB(3,3)=Rs-T*(Xm+Xls);
	BB(3,6)=-T*Xm;

	BB(4,1) = -T*Xm;
	BB(4,2) = (wr-w)*Xm;
	BB(4,4) = Rr-T*(Xm+Xlr);
	BB(4,5) = (wr-w)*(Xm+Xlr);

	BB(5,1) = (w-wr)*Xm;
	BB(5,2) = -T*Xm;
	BB(5,4) = (w-wr)*(Xm+Xlr);
	BB(5,5) = Rr-T*(Xm+Xlr);

	BB(6,3)=-T*Xm;
	BB(6,6)=Rr-T*(Xm+Xlr);

	AA_inv = AA;
	AA_inv.Invert();
}

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
//为下一时步计算诺顿等效电流做准备
void WoundInducGenerator2::calculateDq0Results()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];
	double Vrabc_pu[3];
	double Irabc_pu[3];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
		Vrabc_pu[i] = Vrabc[i]/Vbase*turnratio;
		Irabc_pu[i] = Irabc[i]/Ibase/turnratio;
	}

	//计算park变换矩阵
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	//计算电压
	for (int i=0;i<3;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//计算电流
	for (int i=0;i<3;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}
}
//求解机械方程，计算转子转速
//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
void WoundInducGenerator2::calculateOmega()
{
	//采用风力机拖动方式
	if (control==2)
	{
		//计算风力机转矩
		calculateTwind();
		//风力机转矩为拖动转矩
		TL=-Twind;
	}
	wr_his = wr;
	curAngle_r_1= curAngle_r;
	curAngle_s_1=curAngle_s;

	if (control==0||control==2)
	{
		Te_his=Te;
		Te=0.5*Xm*(Irdq0[0]*Isdq0[1]-Irdq0[1]*Isdq0[0]);
		wr = ((Te+Te_his)/2-TL)*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*Wbase*(w-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w-wr)*Wbase;
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
}

//计算派克变换与反变换矩阵
void WoundInducGenerator2::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
	T(3,1)= 2.0/3.0*1.0/2.0;
	T(3,2)= 2.0/3.0*1.0/2.0;
	T(3,3)= 2.0/3.0*1.0/2.0;
}
void WoundInducGenerator2::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(1,3)= 1.0/2.0;
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(2,3)= 1.0/2.0;
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
	T_(3,3)= 1.0/2.0;
}


// 平均模型相关，xuyin，20121224
// 存储内部变量，以便在迭代校正时恢复
void WoundInducGenerator2::storeInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc_bak[k] = Vsabc[k];
		Vsabc_his_bak[k] = Vsabc_his[k];
		Vsabc_his2_bak[k] = Vsabc_his2[k];
		Isabc_bak[k] = Isabc[k];
		Vrabc_bak[k] = Vrabc[k];		
		Vrabc_his_bak[k] = Vrabc_his[k];
		Vrabc_his2_bak[k] = Vrabc_his2[k];
		Irabc_bak[k] = Irabc[k];
	}

	for (int k=0; k<3; k++) {
		Vsdq0_bak[k] = Vsdq0[k];
		Vsdq0_his_bak[k] = Vsdq0_his[k];
		Vsdq0_his2_bak[k] = Vsdq0_his2[k];
		Isdq0_bak[k] = Isdq0[k];
		Isdq0_his_bak[k] = Isdq0_his[k];
		Isdq0_his2_bak[k] = Isdq0_his2[k];
		Vrdq0_bak[k] = Vrdq0[k];
		Vrdq0_his_bak[k] = Vrdq0_his[k];
		Vrdq0_his2_bak[k] = Vrdq0_his2[k];
		Irdq0_bak[k] = Irdq0[k];
		Irdq0_his_bak[k] = Irdq0_his[k];
		Irdq0_his2_bak[k] = Irdq0_his2[k];
	}

	wr_bak = wr;
	wr_his_bak = wr_his;
	Te_bak = Te;
	Te_his_bak=Te_his;
	curAngle_s_bak = curAngle_s;
	curAngle_r_bak = curAngle_r;
	WindVelocityCounter_bak = WindVelocityCounter;
}

// 迭代校正时恢复内部变量
void WoundInducGenerator2::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Vsabc_his[k] = Vsabc_his_bak[k];
		Vsabc_his2[k] = Vsabc_his2_bak[k];
		Isabc[k] = Isabc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Vrabc_his[k] = Vrabc_his_bak[k];
		Vrabc_his2[k] = Vrabc_his2_bak[k];
		Irabc[k] = Irabc_bak[k];
	}

	for (int k=0; k<3; k++) {
		Vsdq0[k] = Vsdq0_bak[k];
		Vsdq0_his[k] = Vsdq0_his_bak[k];
		Vsdq0_his2[k] = Vsdq0_his2_bak[k];
		Isdq0[k] = Isdq0_bak[k];
		Isdq0_his[k] = Isdq0_his_bak[k];
		Isdq0_his2[k] = Isdq0_his2_bak[k];
		Vrdq0[k] = Vrdq0_bak[k];
		Vrdq0_his[k] = Vrdq0_his_bak[k];
		Vrdq0_his2[k] = Vrdq0_his2_bak[k];
		Irdq0[k] = Irdq0_bak[k];
		Irdq0_his[k] = Irdq0_his_bak[k];
		Irdq0_his2[k] = Irdq0_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te = Te_bak;
	Te_his = Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
	WindVelocityCounter = WindVelocityCounter_bak;
}

void WoundInducGenerator2::calculateNortonEquivalentResistance(double time)
{
	T = 2.0/deltaT/Wbase;
	calculateOmega();
	calculateCoefficientMatrix();
	
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<4;i++)
	{
		for (int j=1;j<4;j++)
		{
			P(i,j)=Ps(i,j);
			P(i+3,j+3)=Pr(i,j);
			invP(i,j)=invPs(i,j);
			invP(i+3,j+3)=invPr(i,j);
		}
	}

	nortonEquivalentConductanceMatrix=invP*AA_inv*P;

	for (int i=0;i<6;i++)
	{
		for (int j=0;j<6;j++)
		{
			nortonEquivalentConductance[6*i+j]=nortonEquivalentConductanceMatrix(i+1,j+1);
		}
	}

	//将等效从有名值换算到标幺值
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			nortonEquivalentConductance[6*i+j]=nortonEquivalentConductance[6*i+j]/Zbase;
			nortonEquivalentConductance[6*(i+3)+j]=nortonEquivalentConductance[6*(i+3)+j]/Zbase*turnratio;
			nortonEquivalentConductance[6*i+(j+3)]=nortonEquivalentConductance[6*i+(j+3)]/Zbase*turnratio;
			nortonEquivalentConductance[6*(i+3)+(j+3)]=nortonEquivalentConductance[6*(i+3)+(j+3)]/Zbase*turnratio*turnratio;
		}
	}
}

// 电机初始化专用函数, xuyin, 20121226
void WoundInducGenerator2::initializeGen(double** GenInitialMatrix, int& ptr)
{
	// 除了第一列为时间t外，每个电机对应GenInitialMatrix的15列
	// 1-3列：Vsabc，4-6列：Vrabc，7-9列：Isabc，10-12列：Irabc，13-15列：Tm，Theta，wr
	// GenInitialMatrix共三行，分别对应当前时步以及前两个时步的值（按照时间顺序）

	// 电气量数值*1000
	for (int i=0; i<3; i++)
		for (int j=0; j<12; j++)
			GenInitialMatrix[i][ptr+j] *= 1000;

	// 读取时间
	double t = GenInitialMatrix[2][0]; //当前时间的坐标为第3行第1列

	// 初始化当前时刻的abc相电压电流
	for (int k=0; k<3; k++) {
		Vsabc[k] = GenInitialMatrix[2][ptr+k];
		Vsabc_his[k] = GenInitialMatrix[1][ptr+k];
		Vsabc_his2[k] = GenInitialMatrix[0][ptr+k];
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Vrabc_his[k] = GenInitialMatrix[1][ptr+3+k];
		Vrabc_his2[k] = GenInitialMatrix[0][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// 机械转矩、转速和转子角
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	curAngle_s = w*Wbase*t;
	//curAngle_s = 0; // curAngle_s的初值可以任意选取
	double Theta = GenInitialMatrix[2][ptr+13];
	curAngle_r = curAngle_s - Theta;

	// 初始化电压电流的dq分量
	double dqMatrix[3][12];

	//将abc坐标下的定子电压电流转化为标幺值
	for (int i=0;i<3;i++)
	{
		for (int j=ptr;j<ptr+3;j++)
			GenInitialMatrix[i][j] /= Vbase;

		for (int j=ptr+3;j<ptr+6;j++)
			GenInitialMatrix[i][j] /= (Vbase/turnratio);

		for (int j=ptr+6;j<ptr+9;j++)
			GenInitialMatrix[i][j] /= Ibase;

		for (int j=ptr+9;j<ptr+12;j++)
			GenInitialMatrix[i][j] /= (Ibase*turnratio);
	}

	// 计算dq分量
	for (int i=0; i<3; i++)
	{
		// 计算park变换角度
		double theta_s = w*Wbase*GenInitialMatrix[i][0];
		double theta_r = theta_s - GenInitialMatrix[i][ptr+13];

		// 计算park变换矩阵
		parkTransMatrix(Ps,theta_s);
		parkTransMatrix(Pr,theta_r);

		// park变换
		// 定子电压
		for (int j=0; j<3; j++)
		{
			dqMatrix[i][j] = Ps(j+1,1)*GenInitialMatrix[i][ptr] + Ps(j+1,2)*GenInitialMatrix[i][ptr+1]
			+ Ps(j+1,3)*GenInitialMatrix[i][ptr+2];
		}
		// 转子电压
		for (int j=3; j<6; j++)
		{
			dqMatrix[i][j] = Pr(j-2,1)*GenInitialMatrix[i][ptr+3] + Pr(j-2,2)*GenInitialMatrix[i][ptr+4]
			+ Pr(j-2,3)*GenInitialMatrix[i][ptr+5];
		}
		// 定子电流
		for (int j=6; j<9; j++)
		{
			dqMatrix[i][j] = Ps(j-5,1)*GenInitialMatrix[i][ptr+6] + Ps(j-5,2)*GenInitialMatrix[i][ptr+7]
			+ Ps(j-5,3)*GenInitialMatrix[i][ptr+8];
		}
		// 转子电压
		for (int j=9; j<12; j++)
		{
			dqMatrix[i][j] = Pr(j-8,1)*GenInitialMatrix[i][ptr+9] + Pr(j-8,2)*GenInitialMatrix[i][ptr+10]
			+ Pr(j-8,3)*GenInitialMatrix[i][ptr+11];
		}
	}

	// dq分量赋值
	for (int k=0; k<3; k++) {
		Vsdq0[k] = dqMatrix[2][k];
		Vsdq0_his[k] = dqMatrix[1][k];
		Vsdq0_his2[k] = dqMatrix[0][k];
		Isdq0[k] = dqMatrix[2][k+6];
		Isdq0_his[k] = dqMatrix[1][k+6];
		Isdq0_his2[k] = dqMatrix[0][k+6];
		Vrdq0[k] = dqMatrix[2][k+3];
		Vrdq0_his[k] = dqMatrix[1][k+3];
		Vrdq0_his2[k] = dqMatrix[0][k+3];
		Irdq0[k] = dqMatrix[2][k+9];
		Irdq0_his[k] = dqMatrix[1][k+9];
		Irdq0_his2[k] = dqMatrix[0][k+9];
	}

	// 节点导纳矩阵计算
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq0_his2[0]*Isdq0_his2[1]-Irdq0_his2[1]*Isdq0_his2[0]);
	Te=0.5*Xm*(Irdq0_his[0]*Isdq0_his[1]-Irdq0_his[1]*Isdq0_his[0]);

	ptr += 15;
}

// 在风力机输入模式下，计算风力机输出转矩
void WoundInducGenerator2::calculateTwind()
{
#ifdef WIND_VELOCITY_DATA_INPUT
	setWindVelocity();
#endif
	TSR=2.237*vWind/(wr*Wbase/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}

//保存电机转速, By Gao Haixiang
void WoundInducGenerator2::saveMachineWr(double** machineWrMatrix, int& ptr, int counter)
{
	machineWrMatrix[counter-1][ptr] = wr;
	ptr++;
}

//	将每个时步的风速数据保存在每台风机中 ,By Gao Haixiang
void WoundInducGenerator2::getWindVelocityData(double** VwMatrix, int rows, int& ptr)
{
	WindVelocityVector = new double[rows];
	for (int i=0;i<rows;i++)
	{
		WindVelocityVector[i] = VwMatrix[i][ptr];
	}
	ptr++;
}

// 设定风速, By Gao Haixiang
void WoundInducGenerator2::setWindVelocity()
{
	vWind = WindVelocityVector[WindVelocityCounter];
	WindVelocityCounter++;
}