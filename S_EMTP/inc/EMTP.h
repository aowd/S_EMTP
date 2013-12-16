#ifndef EMTP_H
#define EMTP_H

#include "Component.h"
#include "SimpleBranch.h"
#include "Resistance.h"
#include "Inductance.h"
#include "Capacitance.h"
#include "AcVoltageSource.h"
#include "SinglePhaseTransformer.h"
#include "Diode.h"
#include "TimeSwitch.h"
#include "SynchGenerator.h"
#include "IGBT.h"
#include "MWSynGenerator12.h"
#include "Motor15.h"
#include "ThreePhaseRectifierBridge.h"
#include "Inverter.h"
#include "DiodeWithSC.h"
#include "Asych5.h"
#include "TimedAcVoltageSource.h"
#include "Impedance.h"
#include "InducGenerator.h"
#include "WoundInducGenerator.h"
#include "WindInducGenerator.h"
#include "newIGBT.h"
#include "ControlledVoltageSource.h"
#include "PWMConverter.h"
#include "ThreePhaseTransformer.h"
#include "PhotoVoltaic.h"
#include "TimedResistance.h"
#include "Battery.h"

#include "MeasureComponent.h"
#include "NodeVoltageMsrComp.h"
#include "BranchVoltageMsrComp.h"
#include "BranchCurrentMsrComp.h"
#include "InducGenWrMsr.h"
#include "InducGenAngleMsr.h"

#include "CtrlComponent.h"
#include "PICtrlComp.h"
#include "SigmaCtrlComp.h"
#include "ConstantCtrlComp.h"
#include "PulseGenComp.h"
#include "TriangleGenComp.h"
#include "SinGenComp.h"
#include "PCtrlComp.h"
#include "T2DTrans.h"
#include "D2TTrans.h"
#include "S2RTrans.h"
#include "R2STrans.h"
#include "PRCoordinate.h"
#include "Comparator.h"
#include "Limiter.h"
#include "TimeConstant.h"
#include "ImpulseGenComp.h"
#include "Sampler.h"
#include "Delay.h"
#include "Selector.h"
#include "SampleHold.h"

#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <cmath>

//是否开启VI文件初始化模式
//#define INITIAL_VI

class EMTP
{
public:
	EMTP();
	~EMTP();
	void initializeSystem();//初始化函数，option=0--启动；option=1--从PSCAD初值启动
	void formConductanceMatrix();
	void calculateBranchNortonEquivalentCurrent(double time);//计算支路诺顿等效电流
	void formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
	void calculateBranchVoltage();
	void calculateBranchCurrent();
	void solveNodeVoltageEquation();//求解节点电压方程
	void saveBranchCurrent();
	void saveMahineWr();
	void advanceTime();//时间前进一个步长
	void saveNodeVoltage();//保存节点电压
	void getCaseDifinitionMatrix(int caseNumber);
	void specifySystem();
	void saveResult(void);
	void checkResult(void);
	void display(void);
	void getPSCADResult(void);
	void getMATLABResult(int option);
	void getCPPResult(int option);

	/* 开关处理相关函数 */
	void switchTreatment();//开关处理主函数
	void checkSwitch(double& ratio,int& switchID,double maxRatio,int& switchCounter);//判断需要动作的开关，以及首先动作的开关对应的ratio
	void doSwitch(double ratio,int switchCounter, int& switchMode);//动作需要动作的开关，重新形成导纳矩阵，并调用插值主函数
	void linearInt(double ratio);//一次线性插值
	void calculateBranchNortonEquivalentCurrent_forBasicComp(double time);//计算支路诺顿等效电流
	void networkSolution();//一次网络求解
	void updateResult(int updateTypeNum);//更新开关处理过程中用于存储结果的变量
	void saveNewResult(int counter);//将插值得到的结果存入结果矩阵中
	void removeChatter(double timeToSwitch);//消除Chatter
	void inverterTreatment();//检查是否有Inverter状态发生变化，若有变化则重新形成导纳矩阵

	/*测量控制系统相关*/
	void getCaseMsrDfnMatrix(int caseNumber);//得到测量系统定义矩阵
	void specifyMsrSystem();//测量系统输入函数
	void getCaseCtrlDfnMatrix(int caseNumber);//得到控制系统定义矩阵
	void specifyCtrlSystem();//控制系统输入函数
	void transferMeasurands();//将电气系统测量值传输至控制系统
	void solveCtrlSystem();//求解控制系统方程
	void transferControlVariables();//将控制信号值传输至电气系统
	void transferControlVariablesForSwitch();//将控制信号值传输至电气系统
	void solveInitCtrlSystem();//求解控制系统初始化方程

	/*PWM变流器平均模型相关*/
	void predictSwitchingInstants(); // 预测开关动作时刻
	void solveCtrlSystemforPrediction(int Nc, double** PWMPredictMatrix); // 预测时求解控制系统
	void correctSwitchingInstants(double tol); // 校正开关动作时刻

	// 求解系统在一个分段子区间上的解,xuyin, 20130224
	void solveG(double * tau_new);
	int check_tau(double * tau_new, double tol);
	void predict_tau();
	void correct_tau(double * tau_new);
	void store_state();
	void restore_state();

	/* 系统计时相关*/
	void my_clock(LARGE_INTEGER &t);//计时
	double cal_timeCost(LARGE_INTEGER &t_start,LARGE_INTEGER &t_end);//计算耗时，返回值单位为us


public:
	/* 仿真设置 */
	double deltaT;//仿真步长
	double saveDeltaT;//存储所用时间步长
	double initDeltaT;//初始化时采用的步长
	double len_subinterval;//平均模型校正的子区间长度
	double startTime;//仿真开始时间
	double finishTime;//仿真结束时间
	int caseNumber;//算例编号
	int counter;//计数器
	int counter2;//控制系统存储所用的计数器
	double curTime;//当前时间
	int Np; // 分段子区间总数
	int Ns; // 一个分段子区间中包含的时步数
	int Nsave;//一个存储子区间中包含的时步数
	double nSolveG;//solveG函数调用总次数
	double ** tau_s;//占空比相关向量
	int cnt_p;//系统当前指针
	/* 系统设置 */
	int nNodes;//节点总数
	int nBranch;//支路总数
	int nColumns;//支路电流的列数
	int rows;//时间，电压电流等全局数组的长度(与结果保存步长saveDeltaT相关)
	int WindRows;//读入风速数组的长度(与仿真步长deltaT相关)
	int caseExistFlag;//算例是否存在的标志。=0，不存在，=1存在
	double costTime;//计算总耗时
	int initializeOption;//初始化选项，0为零值初始化，1为PSCAD初始化
	std::vector<Component*> * branches;//支路数组
	/* 仿真结果存储 */
	TVectorD timeArray;//时间数组
	TVectorD nodeVoltageVec;//节点电压向量
	TVectorD nodeVoltageVec_1;//上一时步节点电压向量_插值时使用
	TVectorD nodeVoltageVec_2;//插值时用于保存nodeVoltageVec_1的信息
	//TMatrixD nodeVoltageMatrix;//节点电压矩阵
	//TMatrixD branchCurrentMatrix;//支路电流矩阵
	double** nodeVoltageMatrix_1;//节点电压矩阵
	double** branchCurrentMatrix_1;//支路电流矩阵
	double** machineWrMatrix;//电机转速保存矩阵
	/* EMTP相关 */
	TMatrixD caseDfnMatrix;//电气系统定义矩阵
	TMatrixD conductanceMatrix;//节点导纳矩阵
	TMatrixD resistanceMatrix;//节点阻抗矩阵
	TDecompLU *lu_conductanceMatrix;//导纳阵的LU分解结果
	TVectorD nodeNortonEquivalentCurrentArray;//节点诺顿等值电流向量

	TMatrixD resultMatrix;//C++程序的计算结果。结果不包含时间列
	TMatrixD resultMatrix0;//其他程序的计算结果，结果不包含时间列

	double testNumber;//for test 09.03.24
	int interpolat_sign;
	/*仿真计时*/
	LARGE_INTEGER tStart,tEnd,tf,tStartAll,tEndAll;
	double t1,t2,t3,t4,t5,t6;

	/************************************************************************/
	/*                        元件分组处理                                  */
	/************************************************************************/
	/* 开关处理相关 */
	int nSwitch;//开关（二极管）总数
	int* diodeNumArray;//所有二极管编号，用于checkSwitch
	int* switchNumArray;//需要动作的开关（二极管）的编号
	double* switchRatioArray;//各个需要动作的开关（二极管）对应的ratio
	int* switchModeArray;
	
	int nInverter;//inverter总数
	int* inverterNumArray;//所有逆变桥编号，用于inverterTreatment

	int nTimeSwitch;//inverter总数
	int* timeSwitchNumArray;//所有逆变桥编号，用于inverterTreatment

	/* 诺顿等值电流计算相关 */
	int nBranch_NEC;//需要计算诺顿等值电流的元件数量
	int* branchNumArray_NEC;//需要计算诺顿等值电流的元件编号

	/* 自定义元件相关 */
	int nBranch_userDef;//自定义元件个数
	int* branchNumArray_userDef;//自定义元件编号

	/*受控源相关*/
	int nControlledBranches;//受控支路个数
	int* controlledBranches;//受控支路编号数组

	/*时变电阻相关*/
	int nTimedResistance;//时变电阻元件个数
	int* timedResistanceNumArray;

	/*控制系统相关*/
	int nCtrlNodes;//控制系统节点总数
	int* nodeCalMark;//控制节点是否计算标志数组
	int* branchCalMark;//控制元件是否计算标志数组
	double* ctrlNodeValue;//控制节点值存储数组
	int nCtrlBranches;//控制支路总数
	std::vector<CtrlComponent*> * ctrlBranches;//控制支路数组
	TMatrixD caseCtrlDfnMatrix;//控制元件定义矩阵
	int nLoopNodes;//闭环节点个数
	int* loopNodeNumber;//闭环节点编号集合
	double* loopNodeValue;//闭环节点值存储集合
	double* loopNodeValue_bak;//闭环节点值存储集合备份

	/* PI控制器相关 */
	int nPIController; // PI控制器数量
	int* PIControllerBranches; // PI控制器支路编号
	double* PIControllerInitialValue; // PI控制器初值

	/* 延迟模块相关 */
	int nDelay; // 延迟模块数量
	int* DelayBranches; // 延迟模块支路编号
	double** DelayInitialValue; // 延迟模块初值

	// xuyin, 20121013
	// double** ctrlResultMatrix; // 存储控制系统的求解结果
	// xuyin, 20130224
	double** ctrlResultMatrix_2; // 仅存储需要存储的控制结果，即迭代校正需要用的信息

	/*测量系统相关*/
	std::vector<MeasureComponent*>* msrComps;//测量系统元件集合
	int nMsrComps;//测量系统元件个数
	TMatrixD caseMsrDfnMatrix;//测量元件定义矩阵

	/*PWM变流器平均模型相关*/
	int nPWMConverter; // PWM变流器数量
	int* branchNumArray_PWMConverter; // PWM变流器编号数组
	int* ctrlNumArray_PWMConverter; // PWM变流器对应的控制节点编号

	/*异步电机初始化相关*/
	int nIndGen;//异步电机台数
	int* branchNumArray_IndGen;//异步电机支路编号数组
};

#endif