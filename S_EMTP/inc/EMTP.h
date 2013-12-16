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

//�Ƿ���VI�ļ���ʼ��ģʽ
//#define INITIAL_VI

class EMTP
{
public:
	EMTP();
	~EMTP();
	void initializeSystem();//��ʼ��������option=0--������option=1--��PSCAD��ֵ����
	void formConductanceMatrix();
	void calculateBranchNortonEquivalentCurrent(double time);//����֧·ŵ�ٵ�Ч����
	void formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
	void calculateBranchVoltage();
	void calculateBranchCurrent();
	void solveNodeVoltageEquation();//���ڵ��ѹ����
	void saveBranchCurrent();
	void saveMahineWr();
	void advanceTime();//ʱ��ǰ��һ������
	void saveNodeVoltage();//����ڵ��ѹ
	void getCaseDifinitionMatrix(int caseNumber);
	void specifySystem();
	void saveResult(void);
	void checkResult(void);
	void display(void);
	void getPSCADResult(void);
	void getMATLABResult(int option);
	void getCPPResult(int option);

	/* ���ش�����غ��� */
	void switchTreatment();//���ش���������
	void checkSwitch(double& ratio,int& switchID,double maxRatio,int& switchCounter);//�ж���Ҫ�����Ŀ��أ��Լ����ȶ����Ŀ��ض�Ӧ��ratio
	void doSwitch(double ratio,int switchCounter, int& switchMode);//������Ҫ�����Ŀ��أ������γɵ��ɾ��󣬲����ò�ֵ������
	void linearInt(double ratio);//һ�����Բ�ֵ
	void calculateBranchNortonEquivalentCurrent_forBasicComp(double time);//����֧·ŵ�ٵ�Ч����
	void networkSolution();//һ���������
	void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���
	void saveNewResult(int counter);//����ֵ�õ��Ľ��������������
	void removeChatter(double timeToSwitch);//����Chatter
	void inverterTreatment();//����Ƿ���Inverter״̬�����仯�����б仯�������γɵ��ɾ���

	/*��������ϵͳ���*/
	void getCaseMsrDfnMatrix(int caseNumber);//�õ�����ϵͳ�������
	void specifyMsrSystem();//����ϵͳ���뺯��
	void getCaseCtrlDfnMatrix(int caseNumber);//�õ�����ϵͳ�������
	void specifyCtrlSystem();//����ϵͳ���뺯��
	void transferMeasurands();//������ϵͳ����ֵ����������ϵͳ
	void solveCtrlSystem();//������ϵͳ����
	void transferControlVariables();//�������ź�ֵ����������ϵͳ
	void transferControlVariablesForSwitch();//�������ź�ֵ����������ϵͳ
	void solveInitCtrlSystem();//������ϵͳ��ʼ������

	/*PWM������ƽ��ģ�����*/
	void predictSwitchingInstants(); // Ԥ�⿪�ض���ʱ��
	void solveCtrlSystemforPrediction(int Nc, double** PWMPredictMatrix); // Ԥ��ʱ������ϵͳ
	void correctSwitchingInstants(double tol); // У�����ض���ʱ��

	// ���ϵͳ��һ���ֶ��������ϵĽ�,xuyin, 20130224
	void solveG(double * tau_new);
	int check_tau(double * tau_new, double tol);
	void predict_tau();
	void correct_tau(double * tau_new);
	void store_state();
	void restore_state();

	/* ϵͳ��ʱ���*/
	void my_clock(LARGE_INTEGER &t);//��ʱ
	double cal_timeCost(LARGE_INTEGER &t_start,LARGE_INTEGER &t_end);//�����ʱ������ֵ��λΪus


public:
	/* �������� */
	double deltaT;//���沽��
	double saveDeltaT;//�洢����ʱ�䲽��
	double initDeltaT;//��ʼ��ʱ���õĲ���
	double len_subinterval;//ƽ��ģ��У���������䳤��
	double startTime;//���濪ʼʱ��
	double finishTime;//�������ʱ��
	int caseNumber;//�������
	int counter;//������
	int counter2;//����ϵͳ�洢���õļ�����
	double curTime;//��ǰʱ��
	int Np; // �ֶ�����������
	int Ns; // һ���ֶ��������а�����ʱ����
	int Nsave;//һ���洢�������а�����ʱ����
	double nSolveG;//solveG���������ܴ���
	double ** tau_s;//ռ�ձ��������
	int cnt_p;//ϵͳ��ǰָ��
	/* ϵͳ���� */
	int nNodes;//�ڵ�����
	int nBranch;//֧·����
	int nColumns;//֧·����������
	int rows;//ʱ�䣬��ѹ������ȫ������ĳ���(�������沽��saveDeltaT���)
	int WindRows;//�����������ĳ���(����沽��deltaT���)
	int caseExistFlag;//�����Ƿ���ڵı�־��=0�������ڣ�=1����
	double costTime;//�����ܺ�ʱ
	int initializeOption;//��ʼ��ѡ�0Ϊ��ֵ��ʼ����1ΪPSCAD��ʼ��
	std::vector<Component*> * branches;//֧·����
	/* �������洢 */
	TVectorD timeArray;//ʱ������
	TVectorD nodeVoltageVec;//�ڵ��ѹ����
	TVectorD nodeVoltageVec_1;//��һʱ���ڵ��ѹ����_��ֵʱʹ��
	TVectorD nodeVoltageVec_2;//��ֵʱ���ڱ���nodeVoltageVec_1����Ϣ
	//TMatrixD nodeVoltageMatrix;//�ڵ��ѹ����
	//TMatrixD branchCurrentMatrix;//֧·��������
	double** nodeVoltageMatrix_1;//�ڵ��ѹ����
	double** branchCurrentMatrix_1;//֧·��������
	double** machineWrMatrix;//���ת�ٱ������
	/* EMTP��� */
	TMatrixD caseDfnMatrix;//����ϵͳ�������
	TMatrixD conductanceMatrix;//�ڵ㵼�ɾ���
	TMatrixD resistanceMatrix;//�ڵ��迹����
	TDecompLU *lu_conductanceMatrix;//�������LU�ֽ���
	TVectorD nodeNortonEquivalentCurrentArray;//�ڵ�ŵ�ٵ�ֵ��������

	TMatrixD resultMatrix;//C++����ļ����������������ʱ����
	TMatrixD resultMatrix0;//��������ļ����������������ʱ����

	double testNumber;//for test 09.03.24
	int interpolat_sign;
	/*�����ʱ*/
	LARGE_INTEGER tStart,tEnd,tf,tStartAll,tEndAll;
	double t1,t2,t3,t4,t5,t6;

	/************************************************************************/
	/*                        Ԫ�����鴦��                                  */
	/************************************************************************/
	/* ���ش������ */
	int nSwitch;//���أ������ܣ�����
	int* diodeNumArray;//���ж����ܱ�ţ�����checkSwitch
	int* switchNumArray;//��Ҫ�����Ŀ��أ������ܣ��ı��
	double* switchRatioArray;//������Ҫ�����Ŀ��أ������ܣ���Ӧ��ratio
	int* switchModeArray;
	
	int nInverter;//inverter����
	int* inverterNumArray;//��������ű�ţ�����inverterTreatment

	int nTimeSwitch;//inverter����
	int* timeSwitchNumArray;//��������ű�ţ�����inverterTreatment

	/* ŵ�ٵ�ֵ����������� */
	int nBranch_NEC;//��Ҫ����ŵ�ٵ�ֵ������Ԫ������
	int* branchNumArray_NEC;//��Ҫ����ŵ�ٵ�ֵ������Ԫ�����

	/* �Զ���Ԫ����� */
	int nBranch_userDef;//�Զ���Ԫ������
	int* branchNumArray_userDef;//�Զ���Ԫ�����

	/*�ܿ�Դ���*/
	int nControlledBranches;//�ܿ�֧·����
	int* controlledBranches;//�ܿ�֧·�������

	/*ʱ��������*/
	int nTimedResistance;//ʱ�����Ԫ������
	int* timedResistanceNumArray;

	/*����ϵͳ���*/
	int nCtrlNodes;//����ϵͳ�ڵ�����
	int* nodeCalMark;//���ƽڵ��Ƿ�����־����
	int* branchCalMark;//����Ԫ���Ƿ�����־����
	double* ctrlNodeValue;//���ƽڵ�ֵ�洢����
	int nCtrlBranches;//����֧·����
	std::vector<CtrlComponent*> * ctrlBranches;//����֧·����
	TMatrixD caseCtrlDfnMatrix;//����Ԫ���������
	int nLoopNodes;//�ջ��ڵ����
	int* loopNodeNumber;//�ջ��ڵ��ż���
	double* loopNodeValue;//�ջ��ڵ�ֵ�洢����
	double* loopNodeValue_bak;//�ջ��ڵ�ֵ�洢���ϱ���

	/* PI��������� */
	int nPIController; // PI����������
	int* PIControllerBranches; // PI������֧·���
	double* PIControllerInitialValue; // PI��������ֵ

	/* �ӳ�ģ����� */
	int nDelay; // �ӳ�ģ������
	int* DelayBranches; // �ӳ�ģ��֧·���
	double** DelayInitialValue; // �ӳ�ģ���ֵ

	// xuyin, 20121013
	// double** ctrlResultMatrix; // �洢����ϵͳ�������
	// xuyin, 20130224
	double** ctrlResultMatrix_2; // ���洢��Ҫ�洢�Ŀ��ƽ����������У����Ҫ�õ���Ϣ

	/*����ϵͳ���*/
	std::vector<MeasureComponent*>* msrComps;//����ϵͳԪ������
	int nMsrComps;//����ϵͳԪ������
	TMatrixD caseMsrDfnMatrix;//����Ԫ���������

	/*PWM������ƽ��ģ�����*/
	int nPWMConverter; // PWM����������
	int* branchNumArray_PWMConverter; // PWM�������������
	int* ctrlNumArray_PWMConverter; // PWM��������Ӧ�Ŀ��ƽڵ���

	/*�첽�����ʼ�����*/
	int nIndGen;//�첽���̨��
	int* branchNumArray_IndGen;//�첽���֧·�������
};

#endif