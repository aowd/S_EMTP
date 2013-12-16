#ifndef ASYCH5_H
#define ASYCH5_H

#include "Component.h"

class Asych5:public Component {
public:
	/************************************************************************/
	/*                构造函数与析构函数                                    */
	/************************************************************************/
	Asych5(int firstNode,int lastNode);//构造函数,输入为5相电机的节点编号
	~Asych5();

	/************************************************************************/
	/*                与EMTP接口的函数                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	virtual void calculateBranchVoltage();//计算支路电压
	virtual void calculateBranchCurrent();//计算支路电流
	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流
	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量

private:
	/************************************************************************/
	/*                电动机内部调用的函数                                  */
	/************************************************************************/
	//计算系数矩阵AA,BB,Rss...Brr,Reqs...inv_Rr,Gs1...Gs4,Gr1...Gr4,...
	virtual void calculateCoefficientMatrix();
	//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
	//为下一时步计算诺顿等效电流做准备
	virtual void calculateDq0Results(double time);
	//求解机械方程，计算转子转速
	//为下一时步求解系数矩阵、计算诺顿等值导纳做准备
	virtual void calculateOmega();
	//计算派克变换与反变换矩阵
	virtual void parkTransMatrix(TMatrixD &P,double Angle);
	virtual void invParkTransMatrix(TMatrixD &INRP,double Angle);
	//判断是否有节点接地
	virtual bool isSeries(int from,int to);

private:
	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//dq0坐标系下的电阻参数
	double Rp,Rr1,Rrt;
	//dq0坐标系下的电抗参数
	double Xl1,Xlr1,Xm1,Xlt1,Xlrt,Xmt,Xs0;
	//中间变量
	double Xs1,Xr1,Xst,Xrt;

	//定子侧电压、电流和角频率的基值
	double Vbase,Ibase,wbase;

	//机械侧参数
	double Je;//转动惯量

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/

	//派克变换矩阵及其逆矩阵
	TMatrixD P;//(3,5)
	TMatrixD invP;//(5,3)

	/*     梯形积分法需要使用的系数矩阵      */
	double coff;//alpha
	double T;//T=1.0/(deltaT*wbase)
	//
	TMatrixD A;//(5,5)
	TMatrixD B;//(5,5)
	TMatrixD AA;//(5,5)
	TMatrixD BB;//(5,5)

	//double B_1[6][6];

	//
	TMatrixD Rss;//AA(1:3,1:3)
	TMatrixD Rsr;//AA(1:3,4:5)
	TMatrixD Rrs;//AA(4:5,1:3)
	TMatrixD Rrr;//AA(4:5,4:5)
	//
	TMatrixD Bss;//BB(1:3,1:3)
	TMatrixD Bsr;//BB(1:3,4:5)
	TMatrixD Brs;//BB(4:5,1:3)
	TMatrixD Brr;//BB(4:5,4:5)

	//处理系数矩阵，使得转换到abcde坐标系时能解耦
	double Tp;
	TMatrixD Rs;//
	TMatrixD Rs1;//
	TMatrixD Rs2;//Rs2=Rs-Rs1;

	//
	TMatrixD Gs1;//
	TMatrixD Gs2;//
	TMatrixD Gs3;//
	TMatrixD Gs4;//

	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	TMatrixD Gr4;//

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	double Vnode[10];//节点电压
	double Vabcde[5];//支路电压
	double Vabcde_1[5];
	double Vabcde_2[5];
	double Iabcde[5];//支路电流
	double Iabcde_1[5];
	double Iabcde_2[5];
	//dq0坐标系下定子线圈的电压电流
	double Vdqs[3];
	double Idqs[3];
	double Idqs_his[3];//预报用
	double Idqs_his2[3];//预报用
	double Idqs_forcast[3];//预报值
	//dq0坐标系下定子侧的磁链
	TVectorD psi;
	TVectorD Isr;//计算磁链时用的电流

	//double psi_1[6];
	//dq0坐标系下转子侧的电压电流
	double Vr[2];
	double Ir[2];
	double Ir_his[2];
	/*  机械侧物理量  */
	double speed0;//转子转速初始值（标幺值）
	double wr0;//转子角频率初始值	
	double wr;//转子转速
	double wr_his;//wr的历史值
	double theta;//同步坐标系d轴领先a轴的角度（Park变换时使用）
	double Te;//电磁转矩
	double Te_his;
	double Tm;//负载机械转矩
	double Tm_his;

	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//等效导纳的历史值，用于更新节点导纳矩阵
	double nortonEquivalentCurrent[5];//abcde坐标下的诺顿等效电流，有名值
	double nortonEquivalentCurrent_1[5];
	double nortonEquivalentCurrent_2[5];
	double nortonEquivalentCurrent_dq0[3];//dq0坐标下的诺顿等效电流，有名值

	//临时变量
	//double temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18;
};

#endif