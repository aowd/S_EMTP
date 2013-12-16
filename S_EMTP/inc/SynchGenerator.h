#ifndef SYNCHGENERATOR_H
#define SYNCHGENERATOR_H

#include "Component.h"

class SynchGenerator:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	SynchGenerator(int a,int b,int c);//���캯��,a��b��cΪ�ڵ���
	~SynchGenerator();

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
	double Ef;//���ŵ�ѹ������ֵ
	double P0,Q0;//��������й����޹�������ֵ
	double Sm;//���ڹ��ʣ�����ֵ���ɼ���õ�
	double Vm;//���ѹ��ֵ������ֵ
	double Im;//�������ֵ������ֵ���ɼ���õ�
	double ph;//A���ѹ�ĳ�ʼ��λ�����ȱ�ʾ
	double Sangle;//�������ؽǣ��ɼ���õ�
	double Vangle;//A���ѹ��ʼ��ǣ��ɼ���õ�
	double Iangle;//A�������ʼ��ǣ��ɼ���õ�
	double Angle;//d������a��ĽǶȵĳ�ʼֵ���ɼ���õ�

	//dq0����ϵ�µĵ������
	double Ra,Rf,RD,RQ;
	//dq0����ϵ�µĵ翹����
	double Xad,Xaq,Xd,Xq,Xf,XD,XQ;
	double XfD,Xaf,XaD,XaQ;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	double Vbase,Ibase,Zbase;//����Zbase���Լ���õ�

	//��е�����
	double H;
	double D;
	double Tm;//����ת��

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//
	double coff;
	double T;//T=1.0/deltaT/w*0.5*(1+coff)

	//���ŵ�ѹ��ϵ��
	double Vf;//Vf=Rf/Xad

	//�ɿ˱任�����������
	TMatrixD P;
	TMatrixD invP;

	//
	TMatrixD AA;//(5,5)
	TMatrixD BB;//(5,5)
	//
	TMatrixD Rss;//AA(1:2,1:2)
	TMatrixD Rsr;//AA(1:2,3:5)
	TMatrixD Rrs;//AA(3:5,1:2)
	TMatrixD Rrr;//AA(3:5,3:5)
	//
	TMatrixD Bss;//BB(1:2,1:2)
	TMatrixD Bsr;//BB(1:2,3:5)
	TMatrixD Brs;//BB(3:5,1:2)
	TMatrixD Brr;//BB(3:5,3:5)

	//����ϵ������ʹ��ת����abc����ϵʱ�ܽ���
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
	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	TMatrixD Gr4;//
	TMatrixD Gr5;//
	TMatrixD Gr15;//Gr15=Gr1+Gr5

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vabc[3];
	double Vabc_1[3];
	double Vabc_2[3];
	double Iabc[3];
	double Iabc_1[3];
	double Iabc_2[3];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vdq0[2];
	double Vdq0_his[2];
	double Vdq0_his2[2];
	double Vdq0_emend[2];
	double Idq0[2];
	double Idq0_his[2];//Ԥ����
	double Idq0_his2[2];//Ԥ����
	double Idq0_forcast[2];//Ԥ��ֵ
	double Idq0_emend[2];
	//dq0���ϵ��ת�Ӳ�ĵ�ѹ����
	double VfDQ[3];
	double IfDQ[3];
	double IfDQ_his[3];
	/*  ��е��������  */
	double wr;//ת��ת��
	double wr_his;//wr����ʷֵ
	double curAngle;//d������a��ĽǶ�
	double curAngle_his;//curAngle����ʷֵ

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