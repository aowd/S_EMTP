#include "Asych5.h"

#include <cmath>
#include <iostream>
using namespace std;

Asych5::Asych5(int firstNode,int lastNode)
{
	//自定义元件
	isUserDef = 1;
	need_NEC=1;
	type=555;
	nPort=10;
	nodeNumber=new int[nPort];
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}

	//dq0坐标系下的电阻参数
	Rp = 0.0086917;
	Rr1 = 0.0047;
	Rrt = 0.161;
	//dq0坐标系下的电抗参数
	Xl1 = 0.202625;
	Xlr1 = 0.077525;
	Xm1 = 1.7819;
	Xlt1 = 0.01;
	Xlrt = 0.1492;
	Xmt = 8.117;
	Xs0 = 0.15045;
	//中间变量
	Xs1 = Xl1 + Xm1; 
	Xr1 = Xlr1 + Xm1;
	Xst = Xlt1 + Xmt;
	Xrt = Xlrt + Xmt;

	//定子侧电压、电流和角频率的基值
	Vbase = 2400;//V,peak value
	Ibase = 2052;//A,peak value
	wbase = 2*PI*20;//fbase=20Hz;

	//机械侧参数
	Je = 25000;//转动惯量


	//派克变换矩阵及其逆矩阵
	P.ResizeTo(1,3,1,5);
	invP.ResizeTo(1,5,1,3);

	/*     梯形积分法需要使用的系数矩阵      */
	coff = 39.0/41.0;//alpha
	T=1.0/(deltaT*wbase);
	//
	A.ResizeTo(1,5,1,5);//(5,5)
	B.ResizeTo(1,5,1,5);//(5,5)
	AA.ResizeTo(1,5,1,5);//(5,5)
	BB.ResizeTo(1,5,1,5);//(5,5)
	//
	Rss.ResizeTo(1,3,1,3);//AA(1:3,1:3)
	Rsr.ResizeTo(1,3,1,2);//AA(1:3,4:5)
	Rrs.ResizeTo(1,2,1,3);//AA(4:5,1:3)
	Rrr.ResizeTo(1,2,1,2);//AA(4:5,4:5)
	//
	Bss.ResizeTo(1,3,1,3);//BB(1:3,1:3)
	Bsr.ResizeTo(1,3,1,2);//BB(1:3,4:5)
	Brs.ResizeTo(1,2,1,3);//BB(4:5,1:3)
	Brr.ResizeTo(1,2,1,2);//BB(4:5,4:5)

	//处理系数矩阵，使得转换到abcde坐标系时能解耦
	Rs.ResizeTo(1,3,1,3);//
	Rs1.ResizeTo(1,3,1,3);//
	Rs2.ResizeTo(1,3,1,3);//Rs2=Rs-Rs1;

	//
	Gs1.ResizeTo(1,3,1,3);//
	Gs2.ResizeTo(1,3,1,3);//
	Gs3.ResizeTo(1,3,1,3);//
	Gs4.ResizeTo(1,3,1,2);//

	//
	Gr1.ResizeTo(1,2,1,2);//
	Gr2.ResizeTo(1,2,1,3);//
	Gr3.ResizeTo(1,2,1,3);//
	Gr4.ResizeTo(1,2,1,2);//

	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	for (int i=0;i<10;i++)
	{
		Vnode[i] = 0;//节点电压
	}
	for (int i=0;i<5;i++)
	{
		Vabcde[i] = 0;//支路电压
		Iabcde[i] = 0;//支路电流
		Vabcde_1[i] = 0;
		Iabcde_1[i] = 0;
		Vabcde_2[i] = 0;
		Iabcde_2[i] = 0;
	}

	//dq0坐标系下定子线圈的电压电流
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		Idqs[i] = 0;
		Idqs_his[i] = 0;//预报用
		Idqs_his2[i] = 0;//预报用
		Idqs_forcast[i] = 0;//预报值
	}

	//dq0坐标系下定子侧的磁链
	psi.ResizeTo(1,5);
	Isr.ResizeTo(1,5);//计算磁链时用的电流

	//dq0坐标系下转子侧的电压电流
	Vr[0] = 0;
	Ir[0] = 0;
	Ir_his[0] = 0;
	Vr[1] = 0;
	Ir[1] = 0;
	Ir_his[1] = 0;
	/*  机械侧物理量  */
	speed0 = 0;;//转子转速初始值（标幺值）
	wr0 = speed0*wbase;//转子角频率初始值	
	wr = 0;//转子转速
	wr_his = 0;//wr的历史值
	theta = 0;//同步坐标系d轴领先a轴的角度（Park变换时使用）
	Te = 0;//电磁转矩
	Te_his = 0;
	Tm = 0;//负载机械转矩
	Tm_his = 0;

	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//等效导纳的历史值，用于更新节点导纳矩阵
	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;//abcde坐标下的诺顿等效电流，有名值
		nortonEquivalentCurrent_1[i] = 0;
		nortonEquivalentCurrent_2[i] = 0;
	}
	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = 0;//dq0坐标下的诺顿等效电流，有名值
	}

	//计算Gs2,Gs3,Gs4
	calculateCoefficientMatrix();
}

//析构函数
Asych5::~Asych5()
{
}

/************************************************************************/
/*                与EMTP接口的函数                                      */
/************************************************************************/
//初始化支路电压电流
void Asych5::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	for (int k=0;k<5;k++)
	{
		Vabcde[k] = 0;
		Iabcde[k] = 0;
	}

	Iabcde[0]=1000.0;
	Iabcde[1]=1000.0;
	Iabcde[2]=0.0;
	Iabcde[3]=-1000.0;
	Iabcde[4]=-1000.0;

	//计算Park变换矩阵
	double Angle;
	Angle = theta;
	parkTransMatrix(P,Angle);

	//计算dq0坐标系下的定子电压
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Vdqs[i] += P(i+1,j+1)*Vabcde[j];
		}
	}

	//计算dq0坐标系下的定子电流
	for (int i=0;i<3;i++)
	{
		Idqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Idqs[i] += P(i+1,j+1)*Iabcde[j];
		}
	}

	for (int k=0;k<3;k++)
	{
		Idqs_his2[k] = Idqs[k];
		Idqs_his[k] = Idqs[k];
	}

	//预报下一时刻定子电流
	for (int i=0;i<3;i++)
	{
		Idqs_forcast[i] = (202.0/120.0)*Idqs[i]-(44.0/120.0)*Idqs_his[i]-(38.0/120.0)*Idqs_his2[i];
	}

	//计算转子电流
	Ir[0] = 0;
	Ir[1] = 0;

	Ir_his[0] = Ir[0];
	Ir_his[1] = Ir[1];

	for (int i=0;i<2;i++)
	{
		Ir[i] = 0;
		for (int j=0;j<3;j++)
		{
			Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	}

	//计算Inew
	TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	tem1.ResizeTo(1,3);
	Idqs_forcast_V.ResizeTo(1,3);
	Idqs_V.ResizeTo(1,3);
	Vdqs_V.ResizeTo(1,3);
	Ir_V.ResizeTo(1,2);
	Inew.ResizeTo(1,3);

	for (int i=1;i<=3;i++)
	{
		Idqs_forcast_V(i) = Idqs_forcast[i-1];
		Idqs_V(i) = Idqs[i-1];
		Vdqs_V(i) = Vdqs[i-1];
	}
	Ir_V(1) = Ir[0];
	Ir_V(2) = Ir[1];

	tem1 = Gs2*Idqs_forcast_V + Gs3*Idqs_V + Gs4*Ir_V + coff*Vdqs_V;
	Inew = Gs1*tem1;

	//计算abc坐标下的诺顿等值电流
	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	Angle = wr*deltaT + theta;
	invParkTransMatrix(invP,Angle);

	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;
		for (int j=0;j<3;j++)
		{
			nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
		}
	}

	//重新计算Idqs，Idqs_his，Idqs_his2
	for (int k=0;k<3;k++)
	{
		Idqs[k] = Inew(k+1);
		Idqs_his[k] = Idqs[k];
		Idqs_his2[k] = Idqs[k];
	}

	//求解机械方程
	//计算Isr
	for (int i=1;i<=3;i++)
	{
		Isr(i) = Idqs[i-1];
	}
	for (int i=1;i<=2;i++)
	{
		Isr(i+3) = Ir[i-1];
	}

	//计算磁链psi
	psi = B*Isr;

	for (int i=1;i<=5;i++)
	{
		psi(i) /= wbase;
	}

	//计算电磁转矩Te
	Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0])*15.0/1e6;

	//计算转子角频率wr
	wr_his = wr;
	wr = wr_his + (Te-Tm)*deltaT/Je;

	//计算转子角theta
	theta = theta+(wr-wr0)*deltaT;

	//用PSCAD数据初始化定子电压电流
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<5;i++)
	{
		Iabcde[i] = initialCurrentArray[ptr+i];
	}
	ptr+=5;
}

//从节点电压数组中读支路两节点电压
void Asych5::readNodeVoltage(TVectorD& nodeVoltageArray)
{//从节点电压数组中读电动机端节点电压
	for (int i=0;i<10;i++)
	{
		if (nodeNumber[i] == 0)
		{
			Vnode[i] = 0;
		}else
		{
			Vnode[i] = nodeVoltageArray(nodeNumber[i]);
		}
	}
}

//计算支路电压
void Asych5::calculateBranchVoltage()
{
	for (int i=0;i<5;i++)
	{
		Vabcde_1[i] = Vabcde[i];
		Vabcde[i] = Vnode[2*i] - Vnode[2*i+1];
	}
}

//计算支路电流
void Asych5::calculateBranchCurrent()
{
	for (int i=0;i<5;i++)
	{
		Iabcde_1[i] = Iabcde[i];
		Iabcde[i] = nortonEquivalentConductance*Vabcde[i] + nortonEquivalentCurrent[i];
	}

	//cout<<endl;
	//for (int i=0;i<5;i++)
	//{
	//	cout<<Vabcde[i]<<endl;
	//}
}

//计算支路的诺顿等效电路中的电流项
void Asych5::calculateNortonEquivalentCurrent(double time)
{
	calculateCoefficientMatrix();
	calculateDq0Results(time);

	TVectorD tem1,Idqs_forcast_V,Idqs_V,Vdqs_V,Ir_V,Inew;
	tem1.ResizeTo(1,3);
	Idqs_forcast_V.ResizeTo(1,3);
	Idqs_V.ResizeTo(1,3);
	Vdqs_V.ResizeTo(1,3);
	Ir_V.ResizeTo(1,2);
	Inew.ResizeTo(1,3);

	for (int i=1;i<=3;i++)
	{
		Idqs_forcast_V(i) = Idqs_forcast[i-1];
		Idqs_V(i) = Idqs[i-1];
		Vdqs_V(i) = Vdqs[i-1];
	}
	Ir_V(1) = Ir[0];
	Ir_V(2) = Ir[1];

	tem1 = Gs2*Idqs_forcast_V + Gs3*Idqs_V + Gs4*Ir_V + coff*Vdqs_V;
	Inew = Gs1*tem1;

	for (int i=0;i<3;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	//cout<<endl<<"Idqs_forcast_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Idqs_forcast_V(i)<<endl;
	//}

	//cout<<endl<<"Idqs_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Idqs_V(i)<<endl;
	//}

	//cout<<endl<<"Ir_V:"<<endl;
	//for (int i=1;i<=2;i++)
	//{
	//	cout<<Ir_V(i)<<endl;
	//}

	//cout<<endl<<"Vdqs_V:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Vdqs_V(i)<<endl;
	//}

	//cout<<endl<<"Inew:"<<endl;
	//for (int i=1;i<=3;i++)
	//{
	//	cout<<Inew(i)<<endl;
	//}

	double Angle;
	Angle = wr0*time + wr*deltaT + theta;
	invParkTransMatrix(invP,Angle);

	for (int k=0;k<5;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	for (int i=0;i<5;i++)
	{
		nortonEquivalentCurrent[i] = 0;
		for (int j=0;j<3;j++)
		{
			nortonEquivalentCurrent[i] += invP(i+1,j+1)*nortonEquivalentCurrent_dq0[j];
		}
	}

	//cout<<endl<<"nortonEquivalentCurrent:"<<endl;
	//for (int i=0;i<5;i++)
	//{
	//	cout<<nortonEquivalentCurrent[i]<<endl;
	//}

	calculateOmega();
}

//形成节点诺顿等效电流向量
void Asych5::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{
	int k,m;
	for (int i=0;i<5;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		nodeNortonEquivalentCurrentArray(k) -= nortonEquivalentCurrent[i];
		nodeNortonEquivalentCurrentArray(m) += nortonEquivalentCurrent[i];
	}
}
void Asych5::formConductanceMatrix(TMatrixD &conductanceMatrix)
{
	int k,m;
	for (int i=0;i<5;i++)
	{
		k = nodeNumber[2*i];
		m = nodeNumber[2*i+1];

		conductanceMatrix(k,k) += nortonEquivalentConductance;
		conductanceMatrix(m,m) += nortonEquivalentConductance;
		conductanceMatrix(k,m) -= nortonEquivalentConductance;
		conductanceMatrix(m,k) -= nortonEquivalentConductance;
	}
}

//保存支路电流
void Asych5::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{
	for (int i=0;i<5;i++)
	{
		branchCurrentMatrix(counter,ptr+i) = Iabcde[i];
	}
	ptr+=5;
}

//保存支路电流
void Asych5::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int i=0;i<5;i++)
	{
		branchCurrentMatrix_1[counter-1][ptr+i] = Iabcde[i];
	}
	ptr+=5;
}

//给定比值，对支路的电压电流进行插值
void Asych5::interpolate(double ratio)
{
	//电压电流变量插值
	for (int k=0;k<5;k++)
	{
		Vabcde[k] = (1-ratio)*Vabcde_1[k] + ratio*Vabcde[k];
		Iabcde[k] = (1-ratio)*Iabcde_1[k] + ratio*Iabcde[k];
	}

	//诺顿等值电流插值
	double temp1,temp2;
	for (int k=0;k<5;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void Asych5::updateResult(int updateTypeNum)
{
	double V_temp[5],I_temp[5];//临时变量，交换数据时使用

	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
		for (int i=0;i<5;i++)
		{
			Vabcde_2[i] = Vabcde_1[i];
			Iabcde_2[i] = Iabcde_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://将_2变量的数值存入_1变量中
		for (int i=0;i<5;i++)
		{
			Vabcde_1[i] = Vabcde_2[i];
			Iabcde_1[i] = Iabcde_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://诺顿等值电流存储：将_1变量的数值存入_2变量中
		for (int i=0;i<5;i++)
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
/*                电动机内部调用的函数                                  */
/************************************************************************/
//计算系数矩阵AA,BB,Rss...Brr,Reqs...inv_Rr,Gs1...Gs4,Gr1...Gr4,...
void Asych5::calculateCoefficientMatrix()
{
	double w;
	w=wr/wbase;

	//计算系数矩阵A和B
	A(1,1)=Rp;
	A(1,2)=-w*Xs1;
	A(1,5)=-w*Xm1;

	A(2,1)=w*Xs1;
	A(2,2)=Rp;
	A(2,4)=w*Xm1;

	A(3,3)=Rp;

	A(4,4)=Rr1;
	A(5,5)=Rr1;

	B(1,1)=Xs1;
	B(1,4)=Xm1;

	B(2,2)=Xs1;
	B(2,5)=Xm1;

	B(3,3)=Xs0;

	B(4,1)=Xm1;
	B(4,4)=Xr1;

	B(5,2)=Xm1;
	B(5,5)=Xr1;

	//计算系数矩阵AA和BB
	AA=A+(1.0+coff)*T*B;
	BB=coff*A-(1.0+coff)*T*B;

	//计算系数矩阵Rss...Brr
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rsr(i,j) = AA(i,j+3);
			Bsr(i,j) = BB(i,j+3);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rrs(i,j) = AA(i+3,j);
			Brs(i,j) = BB(i+3,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrr(i,j) = AA(i+3,j+3);
			Brr(i,j) = BB(i+3,j+3);
		}
	}

	//计算系数矩阵Gr1...Gr4
	Gr1(1,1)=-1.0/Rrr(1,1);
	Gr1(2,2)=Gr1(1,1);

	Gr2 = Gr1*Rrs;
	Gr3 = Gr1*Brs;
	Gr4 = Gr1*Brr;

	//计算系数矩阵Reqs...inv_Rr,Gs1...Gs4
	Rs = Rss+Rsr*Gr2;
	Tp = Rs(1,1);

	Rs1(1,1)=Tp;
	Rs1(2,2)=Tp;
	Rs1(3,3)=Tp;

	Rs2 = Rs-Rs1;
 
	nortonEquivalentConductance = 1/Tp;
	//cout<<nortonEquivalentConductance<<endl;

	Gs1(1,1)=1.0/Rs(1,1);
	Gs1(2,2)=Gs1(1,1);
	Gs1(3,3)=Gs1(1,1);

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Gs2(i,j) = -1*Rs2(i,j);
		}
	}

	Gs3 = Bss+Rsr*Gr3;
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Gs3(i,j) = -1*Gs3(i,j);
		}
	}

	Gs4 = Bsr+Rsr*Gr4;
	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Gs4(i,j) = -1*Gs4(i,j);
		}
	}
}

//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
void Asych5::calculateDq0Results(double time)
{
	//计算Park变换矩阵
	double Angle;
	Angle = wr0*time+theta;
	parkTransMatrix(P,Angle);

	//计算dq0坐标系下的定子电压
	for (int i=0;i<3;i++)
	{
		Vdqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Vdqs[i] += P(i+1,j+1)*Vabcde[j];
		}
	}

	//计算dq0坐标系下的定子电流
	for (int i=0;i<3;i++)
	{
		Idqs_his2[i] = Idqs_his[i];
		Idqs_his[i] = Idqs[i];
		Idqs[i] = 0;
		for (int j=0;j<5;j++)
		{
			Idqs[i] += P(i+1,j+1)*Iabcde[j];
		}
	}

	//预报下一时刻定子电流
	for (int i=0;i<3;i++)
	{
		Idqs_forcast[i] = (202.0/120.0)*Idqs[i]-(44.0/120.0)*Idqs_his[i]-(38.0/120.0)*Idqs_his2[i];
	}

	//计算转子电流
	Ir_his[0] = Ir[0];
	Ir_his[1] = Ir[1];

	for (int i=0;i<2;i++)
	{
		Ir[i] = 0;
		for (int j=0;j<3;j++)
		{
			Ir[i] = Ir[i] + Gr2(i+1,j+1)*Idqs[j] + Gr3(i+1,j+1)*Idqs_his[j];
		}
		Ir[i] = Ir[i] + Gr4(i+1,1)*Ir_his[0] + Gr4(i+1,2)*Ir_his[1];
	}

	//for (int i=0;i<2;i++)
	//{
	//	cout<<Ir[i]<<endl;
	//}
}

//求解机械方程，计算转子转速
void Asych5::calculateOmega()
{
	////计算Isr
	for (int i=1;i<=3;i++)
	{
		Isr(i) = Idqs[i-1];
	}
	for (int i=1;i<=2;i++)
	{
		Isr(i+3) = Ir[i-1];
	}

	//计算磁链psi
	psi = B*Isr;

	for (int i=1;i<=5;i++)
	{
		psi(i) /= wbase;
	}

	//计算电磁转矩Te
	Te = (psi(1)*Idqs[1]-psi(2)*Idqs[0])*15.0/1e6;

	//计算转子角频率wr
	wr_his = wr;
	wr = wr_his + (Te-Tm)*deltaT/Je;

	//计算转子角theta
	theta = theta+(wr-wr0)*deltaT;

	//cout<<wr<<endl;
	//cout<<theta<<endl;
}

//计算派克变换矩阵
void Asych5::parkTransMatrix(TMatrixD &P,double Angle)
{
	P(1,1)=2.0/5.0*cos(Angle);
	P(1,2)=2.0/5.0*cos(Angle-2.0/5.0*PI);
	P(1,3)=2.0/5.0*cos(Angle-4.0/5.0*PI);
	P(1,4)=2.0/5.0*cos(Angle+4.0/5.0*PI);
	P(1,5)=2.0/5.0*cos(Angle+2.0/5.0*PI);

	P(2,1)=-2.0/5.0*sin(Angle);
	P(2,2)=-2.0/5.0*sin(Angle-2.0/5.0*PI);
	P(2,3)=-2.0/5.0*sin(Angle-4.0/5.0*PI);
	P(2,4)=-2.0/5.0*sin(Angle+4.0/5.0*PI);
	P(2,5)=-2.0/5.0*sin(Angle+2.0/5.0*PI);

	P(3,1)=1.0/5.0;
	P(3,2)=1.0/5.0;
	P(3,3)=1.0/5.0;
	P(3,4)=1.0/5.0;
	P(3,5)=1.0/5.0;
}

//计算派克反变换矩阵
void Asych5::invParkTransMatrix(TMatrixD &INRP,double Angle)
{
	INRP(1,1)=cos(Angle);
	INRP(2,1)=cos(Angle-2.0/5.0*PI);
	INRP(3,1)=cos(Angle-4.0/5.0*PI);
	INRP(4,1)=cos(Angle+4.0/5.0*PI);
	INRP(5,1)=cos(Angle+2.0/5.0*PI);

	INRP(1,2)=-sin(Angle);
	INRP(2,2)=-sin(Angle-2.0/5.0*PI);
	INRP(3,2)=-sin(Angle-4.0/5.0*PI);
	INRP(4,2)=-sin(Angle+4.0/5.0*PI);
	INRP(5,2)=-sin(Angle+2.0/5.0*PI);

	INRP(1,3)=1.0;
	INRP(2,3)=1.0;
	INRP(3,3)=1.0;
	INRP(4,3)=1.0;
	INRP(5,3)=1.0;
}


//判断是否有节点接地
bool Asych5::isSeries(int from,int to)
{
	int series = 1;

	if(to==0) {to=from;series=0;}
	if(from==0) {from=to;series=0;}
	if(to==0) cerr<<"Both Nodes Zero!!"<<endl;

	return series;
}