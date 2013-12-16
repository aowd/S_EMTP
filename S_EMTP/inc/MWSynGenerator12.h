#ifndef MWSYNGENERATOR12_H
#define MWSYNGENERATOR12_H

#include "Component.h"

class MWSynGenerator12:public Component {
public:
	/************************************************************************/
	/*                构造函数与析构函数                                    */
	/************************************************************************/
	MWSynGenerator12(int firstNode,int lastNode);
	~MWSynGenerator12();

	/************************************************************************/
	/*                与EMTP接口的函数                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//初始化支路电压电流
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//从节点电压数组中读支路两节点电压
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);//计算支路的诺顿等效电路中的电流项
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//形成节点诺顿等效电流向量
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//形成节点导纳阵
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//保存支路电流
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//保存支路电流

	virtual void interpolate(double ratio);//给定比值，对支路的电压电流进行插值
	virtual void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量

private:
	/************************************************************************/
	/*                发电机内部调用的函数                                  */
	/************************************************************************/
	//计算系数矩阵AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
	virtual void calculateCoefficientMatrix();
	//计算dq0坐标系下定子和转子的电压电流，并预报下一时刻定子电流
	//为下一时步计算诺顿等效电流做准备
	virtual void calculateDq0Results();
	//计算派克变换与反变换矩阵
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
	
	void parkTransMatrix(double P[8][3],double angle);
	void invParkTransMatrix(double invP[12][2],double angle);

private:
	/************************************************************************/
	/*                需要从外部输入的参数                                  */
	/************************************************************************/
	//用于初始化的参数
	double w;//电角速度，有名值
	double Ef;//励磁电压，标幺值
	double P0,Q0;//发电机的有功和无功，有名值
	double Sm;//视在功率，标幺值，可计算得到
	double Vm;//相电压幅值，标幺值
	double Im;//相电流幅值，标幺值，可计算得到
	double ph;//A1相电压的初始相位，弧度表示
	double Sangle;//功率因素角，可计算得到
	double Vangle;//A1相电压初始相角，可计算得到
	double Iangle;//A1相电流初始相角，可计算得到
	double Angle;//d轴领先a1轴的角度的初始值，可计算得到//added by chenlj 启用该变量，含义为初始夹角

	//dq0坐标系下的电阻参数
	double Ra,Rf,RD,RQ;
	//dq0坐标系下的电抗参数
	double Xad,Xaq,Xd,Xq,Xf,XD,XQ;
	double Xdm1,Xdm2,Xqm1,Xqm2,XfD,Xaf,XaD,XaQ;

	//定子侧电压、电流和阻抗的基值
	double Vbase,Ibase,Zbase;//其中Zbase可以计算得到

	//机械侧参数
	double H;
	double D;
	double Tm;//输入转矩

	/************************************************************************/
	/*				系数及系数矩阵,各矩阵的意义见相关文档                   */
	/************************************************************************/
	//
	double coff;
	double T;//T=1.0/deltaT/w*0.5*(1+coff)

	//励磁电压的系数
	double Vf;//Vf=Rf/Xad

	//派克变换矩阵及其逆矩阵
	TMatrixD P;
	TMatrixD invP;
	double P_1[8][3];
	double invP_1[12][2];

	//
	TMatrixD AA;//(11,11)
	TMatrixD BB;//(11,11)
	//
	TMatrixD Rss;//AA(1:8,1:8)
	TMatrixD Rsr;//AA(1:8,9:11)
	TMatrixD Rrs;//AA(9:11,1:8)
	TMatrixD Rrr;//AA(9:11,9:11)
	//
	TMatrixD Bss;//BB(1:8,1:8)
	TMatrixD Bsr;//BB(1:8,9:11)
	TMatrixD Brs;//BB(9:11,1:8)
	TMatrixD Brr;//BB(9:11,9:11)

	//处理系数矩阵，使得转换到abc坐标系时能解耦
	TMatrixD Rs;//
	TMatrixD R_ave;//
	TMatrixD R_res;//R_res=Rs-Rave;
	TMatrixD inv_Rs;//inv(R_ave)
	TMatrixD inv_Rr;//inv(Rrr)

	//
	TMatrixD Gs1;//
	TMatrixD Gs2;//
	TMatrixD Gs3;//
	TMatrixD Gs4;//
	TMatrixD Gs5;//
	TMatrixD Gs6;//
	TMatrixD Gs13;//Gs13=Gs1+Gs3
	TMatrixD Gs25;//Gs25=Gs2+Gs5

	double Gs_1[8][8];
	double Gs_3[8][8];
	double Gs_4[8][3];
	double Gs_6[8][8];
	double Gs_25[8][3];

	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	TMatrixD Gr4;//
	TMatrixD Gr5;//
	TMatrixD Gr15;//Gr15=Gr1+Gr5
	
	double Gr_2[3][8];
	double Gr_3[3][8];
	double Gr_4[3][3];
	double Gr_15[3][3];

	/************************************************************************/
	/*              各时步需要计算的物理量                                  */
	/************************************************************************/
	//以下各个物理量均为标幺值
	/*  电气侧物理量  */
	//abc坐标系下定子线圈的电压电流
	double Vabc[12];
	double Vabc_pu[12];
	double Vabc_1[12];
	double Vabc_2[12];
	double Iabc[12];
	double Iabc_pu[12];
	double Iabc_1[12];
	double Iabc_2[12];
	//dq0坐标系下定子线圈的电压电流
	double Vdq0[8];
	double Vdq0_his[8];
	double Vdq0_his2[8];
	double Idq0[8];
	double Idq0_his[8];//预报用
	double Idq0_his2[8];//预报用
	double Idq0_forcast[8];//预报值
	//dq0左边系下转子侧的电压电流
	double VfDQ[3];
	double IfDQ[3];
	double IfDQ_his[3];
	/*  机械侧物理量  */
	double wr;//转子转速
	double wr_his;//wr的历史值
	double curAngle;//d轴领先a轴的角度

	/************************************************************************/
	/*              等效电路中的参数                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//等效导纳的历史值，用于更新节点导纳矩阵
	double nortonEquivalentCurrent[12];//abc坐标下的诺顿等效电流，有名值
	double nortonEquivalentCurrent_1[12];
	double nortonEquivalentCurrent_2[12];
	double nortonEquivalentCurrent_dq0[8];//dq0坐标下的诺顿等效电流，标幺值
};

#endif