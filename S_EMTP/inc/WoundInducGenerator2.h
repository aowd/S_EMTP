#ifndef WOUNDINDUCGENERATOR2_H
#define WOUNDINDUCGENERATOR2_H

#include "Component.h"

//�Ƿ���Ҫ��������ļ�
//#define WIND_VELOCITY_DATA_INPUT	//��Ҫ�������ļ��������������

class WoundInducGenerator2:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	WoundInducGenerator2(int id, int firstNode,int lastNode, int control,double Vw);//���캯��,a��b��cΪ�ڵ���
	~WoundInducGenerator2();

	/************************************************************************/
	/*                ��EMTP�ӿڵĺ���                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage();//����֧·��ѹ
	virtual void calculateBranchCurrent();//����֧·����
	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void calculateNortonEquivalentResistance(double time);//����֧·��ŵ�ٵ�Ч����
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����
	
	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���

	// ƽ��ģ����أ�xuyin��20121224
	virtual void storeInternalVariables(); // �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
	virtual void restoreInternalVariables(); // ����У��ʱ�ָ��ڲ�����

	// �����ʼ��ר�ú���, xuyin, 20121226
	virtual void initializeGen(double** GenInitialMatrix, int& ptr);
	//	��ÿ��ʱ���ķ������ݱ�����ÿ̨����� ,By Gao Haixiang
	virtual void getWindVelocityData(double** VwMatrix, int rows, int& ptr);

	//������ת��, By Gao Haixiang
	virtual void saveMachineWr(double** machineWrMatrix, int& ptr, int counter);

private:
	/************************************************************************/
	/*                ������ڲ����õĺ���                                  */
	/************************************************************************/
	//����ϵ������
	virtual void calculateCoefficientMatrix();
	//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ����
	//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
	virtual void calculateDq0Results();
	//����е���̣�����ת��ת��
	//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
	virtual void calculateOmega();
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
	//������������ת��
	virtual void calculateTwind();
	//�趨���� , By Gao Haixiang
	virtual void setWindVelocity();
public:
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	double Wbase;//����ٶȣ�����ֵ
	double w;//����ٶȣ�����ֵ

	//dq0����ϵ�µĵ������
	double Rs,Rr;
	//dq0����ϵ�µĵ翹����
	double Xm,Xls,Xlr;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	double Sbase;
	double Vbase,Ibase,Zbase;//����Zbase���Լ���õ�

	//��ת��������
	double turnratio;

	//control=1,ת�ٿ���
	//control=0,ת�ؿ���
	//control=2,���ٿ���
	int control;

	//��е�����
	double H;//��е����
	double D;//����ϵ����δʵ�֣�
	double TL;//����ת��
	
	//�������
	double vWind;//����
	double p_air;//Air Density
	double rArea;//Rotor Area
	double GR;//Gearbox Ratio
	double Pwind;//output Power
	double Twind;//output Torque
	double Cp;//Coefficient of Power
	double beta;//Pitch Angle
	double TSR;//Tip Speed Raio
	double GE;//Gearbox Efficiency

	double* WindVelocityVector;//����������ݵĴ洢����
	int WindVelocityCounter;//���������ȡָ��


	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//
	double T;//T=2.0/deltaT/Wbase

	//�ɿ˱任�����������
	TMatrixD Ps,Pr,P;
	TMatrixD invPs,invPr,invP;

	//
	TMatrixD AA;//(4,4)
	TMatrixD BB;//(4,4)
	TMatrixD AA_inv;//

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vsabc[3];
	double Vsabc_1[3];
	double Vsabc_2[3];
	double Vsabc_his[3];
	double Vsabc_his2[3];
	double Isabc[3];
	double Isabc_1[3];
	double Isabc_2[3];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vsdq0[3];
	double Vsdq0_his[3];
	double Vsdq0_his2[3];
	double Isdq0[3];
	double Isdq0_his[3];
	double Isdq0_his2[3];
	//abc����ϵ��ת����Ȧ�ĵ�ѹ����
	double Vrabc[3];
	double Vrabc_1[3];
	double Vrabc_2[3];
	double Vrabc_his[3];
	double Vrabc_his2[3];
	double Irabc[3];
	double Irabc_1[3];
	double Irabc_2[3];
	//dq0����ϵ��ת�Ӳ�ĵ�ѹ����
	double Vrdq0[3];
	double Vrdq0_his[3];
	double Vrdq0_his2[3];
	double Irdq0[3];
	double Irdq0_his[3];
	double Irdq0_his2[3];


	/*  ��е��������  */
	double wr;//ת��ת��
	double wr_his;//wr����ʷֵ
	double Te;// ����ת��
	double Te_his; //����ת����ʷֵ
	double curAngle_s;//d�����ȶ���a��ĽǶ�
	double curAngle_s_1;// ���ӽǶ���ʷֵ
	double curAngle_r;//d������ת��a��ĽǶ�
	double curAngle_r_1;	// ת�ӽǶ���ʷֵ

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	TMatrixD nortonEquivalentConductanceMatrix;
	double nortonEquivalentConductance[36];//����ֵ
	double nortonEquivalentConductance_pu[36];//����ֵ
	double nortonEquivalentCurrent[6];//abc�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[6];
	double nortonEquivalentCurrent_2[6];
	double nortonEquivalentCurrent_dq0[6];//dq�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_dq0_his[6];
	

	// ���ƽ��ģ�͵ĵ���У���㷨ʹ��, xuyin, 20121224
	// �����ڲ�����ʱʹ��
	double Vsabc_bak[3];
	double Vsabc_his_bak[3];
	double Vsabc_his2_bak[3];
	double Isabc_bak[3];
	double Vsdq0_bak[3];
	double Vsdq0_his_bak[3];
	double Vsdq0_his2_bak[3];
	double Isdq0_bak[3];
	double Isdq0_his_bak[3];
	double Isdq0_his2_bak[3];
	double Vrabc_bak[3];
	double Vrabc_his_bak[3];
	double Vrabc_his2_bak[3];
	double Irabc_bak[3];
	double Vrdq0_bak[3];
	double Vrdq0_his_bak[3];
	double Vrdq0_his2_bak[3];
	double Irdq0_bak[3];
	double Irdq0_his_bak[3];
	double Irdq0_his2_bak[3];
	double wr_bak;
	double wr_his_bak;
	double Te_bak;
	double Te_his_bak;
	double curAngle_s_bak;
	double curAngle_r_bak;
	int WindVelocityCounter_bak;//���������ȡָ�뱸��
	double nortonEquivalentCurrent_bak[6];
	double nortonEquivalentCurrent_dq0_bak[6];

};

#endif