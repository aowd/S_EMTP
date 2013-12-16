#ifndef INDUCGENERATOR_H
#define INDUCGENERATOR_H

#include "Component.h"

class InducGenerator:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	InducGenerator(int a,int b,int c);//���캯��,a��b��cΪ�ڵ���
	~InducGenerator();

	/************************************************************************/
	/*                ��EMTP�ӿڵĺ���                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage();//����֧·��ѹ
	virtual void calculateBranchCurrent();//����֧·����
	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����
	virtual void updateConductanceMatrix(TMatrixD &conductanceMatrix);//ת��ת�ٱ仯ʱ���µ��ɾ���

	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���

private:
	/************************************************************************/
	/*                ������ڲ����õĺ���                                  */
	/************************************************************************/
	//����ϵ������AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
	virtual void calculateCoefficientMatrix();
	//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
	//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
	virtual void calculateDq0Results();
	//����е���̣�����ת��ת��
	//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
	virtual void calculateOmega();
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
public:
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	double w;//����ٶȣ�����ֵ
	double w_pu;
	//*****************
	//��ν��г�ʼ��
	//*****************

	//dq0����ϵ�µĵ������
	double Rs,Rr;
	//dq0����ϵ�µĵ翹����
	double Xm,Xls,Xlr;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	double Vbase,Ibase,Zbase;//����Zbase���Լ���õ�

	//��ת��������
	double turnratio;

	//control=1,ת�ٿ���
	//control=0,ת�ؿ���
	int control;

	//��е�����
	double H;
	double D;
	double Tm;//����ת��

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//
	double coff;
	double T;//T=1.0/deltaT/w*(1+coff)

	//�ɿ˱任�����������
	TMatrixD Ps,Pr;
	TMatrixD invPs,invPr;

	//
	TMatrixD AA;//(4,4)
	TMatrixD BB;//(4,4)
	//
	TMatrixD Rss;//AA(1:2,1:2)
	TMatrixD Rsr;//AA(1:2,3:4)
	TMatrixD Rrs;//AA(3:4,1:2)
	TMatrixD Rrr;//AA(3:4,3:4)
	//
	TMatrixD Bss;//BB(1:2,1:2)
	TMatrixD Bsr;//BB(1:2,3:4)
	TMatrixD Brs;//BB(3:4,1:2)
	TMatrixD Brr;//BB(3:4,3:4)

	//����ϵ������ʹ��ת����abc����ϵʱ�ܽ���
	TMatrixD R_s;//
	TMatrixD R_ave;//
	TMatrixD R_res;//R_res=Rs-Rave;
	TMatrixD inv_Rs;//inv(R_ave)
	TMatrixD inv_Rr;//inv(Rrr)

	//
	TMatrixD Gs1;//
	TMatrixD Gs2;//
	TMatrixD Gs3;//
	TMatrixD Gs4;//
	
	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vsabc[3];
	double Vsabc_1[3];
	double Vsabc_2[3];
	double Isabc[3];
	double Isabc_1[3];
	double Isabc_2[3];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vsdq0[2];
	double Vsdq0_his[2];
	double Vsdq0_his2[2];
	double Isdq0[2];
	double Isdq0_his[2];//Ԥ����
	double Isdq0_his2[2];//Ԥ����
	double Isdq0_forcast[2];//Ԥ��ֵ
	//abc����ϵ��ת����Ȧ�ĵ�ѹ����
	double Vrabc[3];
	double Vrabc_1[3];
	double Vrabc_2[3];
	double Irabc[3];
	double Irabc_1[3];
	double Irabc_2[3];
	//dq0����ϵ��ת�Ӳ�ĵ�ѹ����
	double Vrdq0[2];
	double Vrdq0_his[2];
	double Vrdq0_his2[2];
	double Irdq0[2];
	double Irdq0_his[2];//Ԥ����
	double Irdq0_his2[2];//Ԥ����
	double Irdq0_forcast[2];//Ԥ��ֵ
	/*  ��е��������  */
	double wr;//ת��ת��
	double wr_his;//wr����ʷֵ
	double curAngle_s;//d�����ȶ���a��ĽǶ�
	double curAngle_r;//d������ת��a��ĽǶ�

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	double nortonEquivalentCurrent[3];//abc�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[3];
	double nortonEquivalentCurrent_2[3];
	double nortonEquivalentCurrent_dq0[2];//dq0�����µ�ŵ�ٵ�Ч����������ֵ
};

#endif