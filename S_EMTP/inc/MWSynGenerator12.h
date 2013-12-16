#ifndef MWSYNGENERATOR12_H
#define MWSYNGENERATOR12_H

#include "Component.h"

class MWSynGenerator12:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	MWSynGenerator12(int firstNode,int lastNode);
	~MWSynGenerator12();

	/************************************************************************/
	/*                ��EMTP�ӿڵĺ���                                      */
	/************************************************************************/
	virtual void initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time);//��ʼ��֧·��ѹ����
	virtual void readNodeVoltage(TVectorD& nodeVoltageArray);//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
	virtual void calculateBranchVoltage();
	virtual void calculateBranchCurrent();
	virtual void calculateNortonEquivalentCurrent(double time);//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	virtual void formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray);//�γɽڵ�ŵ�ٵ�Ч��������
	virtual void formConductanceMatrix(TMatrixD &conductanceMatrix);//�γɽڵ㵼����
	virtual void saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter);//����֧·����
	virtual void saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter);//����֧·����

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
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
	
	void parkTransMatrix(double P[8][3],double angle);
	void invParkTransMatrix(double invP[12][2],double angle);

private:
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
	double ph;//A1���ѹ�ĳ�ʼ��λ�����ȱ�ʾ
	double Sangle;//�������ؽǣ��ɼ���õ�
	double Vangle;//A1���ѹ��ʼ��ǣ��ɼ���õ�
	double Iangle;//A1�������ʼ��ǣ��ɼ���õ�
	double Angle;//d������a1��ĽǶȵĳ�ʼֵ���ɼ���õ�//added by chenlj ���øñ���������Ϊ��ʼ�н�

	//dq0����ϵ�µĵ������
	double Ra,Rf,RD,RQ;
	//dq0����ϵ�µĵ翹����
	double Xad,Xaq,Xd,Xq,Xf,XD,XQ;
	double Xdm1,Xdm2,Xqm1,Xqm2,XfD,Xaf,XaD,XaQ;

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
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vabc[12];
	double Vabc_pu[12];
	double Vabc_1[12];
	double Vabc_2[12];
	double Iabc[12];
	double Iabc_pu[12];
	double Iabc_1[12];
	double Iabc_2[12];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vdq0[8];
	double Vdq0_his[8];
	double Vdq0_his2[8];
	double Idq0[8];
	double Idq0_his[8];//Ԥ����
	double Idq0_his2[8];//Ԥ����
	double Idq0_forcast[8];//Ԥ��ֵ
	//dq0���ϵ��ת�Ӳ�ĵ�ѹ����
	double VfDQ[3];
	double IfDQ[3];
	double IfDQ_his[3];
	/*  ��е��������  */
	double wr;//ת��ת��
	double wr_his;//wr����ʷֵ
	double curAngle;//d������a��ĽǶ�

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	double nortonEquivalentCurrent[12];//abc�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[12];
	double nortonEquivalentCurrent_2[12];
	double nortonEquivalentCurrent_dq0[8];//dq0�����µ�ŵ�ٵ�Ч����������ֵ
};

#endif