#ifndef ASYCH5_H
#define ASYCH5_H

#include "Component.h"

class Asych5:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	Asych5(int firstNode,int lastNode);//���캯��,����Ϊ5�����Ľڵ���
	~Asych5();

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
	virtual void interpolate(double ratio);//������ֵ����֧·�ĵ�ѹ�������в�ֵ
	virtual void updateResult(int updateTypeNum);//���¿��ش�����������ڴ洢����ı���

private:
	/************************************************************************/
	/*                �綯���ڲ����õĺ���                                  */
	/************************************************************************/
	//����ϵ������AA,BB,Rss...Brr,Reqs...inv_Rr,Gs1...Gs4,Gr1...Gr4,...
	virtual void calculateCoefficientMatrix();
	//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
	//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
	virtual void calculateDq0Results(double time);
	//����е���̣�����ת��ת��
	//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
	virtual void calculateOmega();
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &P,double Angle);
	virtual void invParkTransMatrix(TMatrixD &INRP,double Angle);
	//�ж��Ƿ��нڵ�ӵ�
	virtual bool isSeries(int from,int to);

private:
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//dq0����ϵ�µĵ������
	double Rp,Rr1,Rrt;
	//dq0����ϵ�µĵ翹����
	double Xl1,Xlr1,Xm1,Xlt1,Xlrt,Xmt,Xs0;
	//�м����
	double Xs1,Xr1,Xst,Xrt;

	//���Ӳ��ѹ�������ͽ�Ƶ�ʵĻ�ֵ
	double Vbase,Ibase,wbase;

	//��е�����
	double Je;//ת������

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/

	//�ɿ˱任�����������
	TMatrixD P;//(3,5)
	TMatrixD invP;//(5,3)

	/*     ���λ��ַ���Ҫʹ�õ�ϵ������      */
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

	//����ϵ������ʹ��ת����abcde����ϵʱ�ܽ���
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
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vnode[10];//�ڵ��ѹ
	double Vabcde[5];//֧·��ѹ
	double Vabcde_1[5];
	double Vabcde_2[5];
	double Iabcde[5];//֧·����
	double Iabcde_1[5];
	double Iabcde_2[5];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vdqs[3];
	double Idqs[3];
	double Idqs_his[3];//Ԥ����
	double Idqs_his2[3];//Ԥ����
	double Idqs_forcast[3];//Ԥ��ֵ
	//dq0����ϵ�¶��Ӳ�Ĵ���
	TVectorD psi;
	TVectorD Isr;//�������ʱ�õĵ���

	//double psi_1[6];
	//dq0����ϵ��ת�Ӳ�ĵ�ѹ����
	double Vr[2];
	double Ir[2];
	double Ir_his[2];
	/*  ��е��������  */
	double speed0;//ת��ת�ٳ�ʼֵ������ֵ��
	double wr0;//ת�ӽ�Ƶ�ʳ�ʼֵ	
	double wr;//ת��ת��
	double wr_his;//wr����ʷֵ
	double theta;//ͬ������ϵd������a��ĽǶȣ�Park�任ʱʹ�ã�
	double Te;//���ת��
	double Te_his;
	double Tm;//���ػ�еת��
	double Tm_his;

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	double nortonEquivalentConductance;
	double nortonEquivalentConductance_his;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	double nortonEquivalentCurrent[5];//abcde�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[5];
	double nortonEquivalentCurrent_2[5];
	double nortonEquivalentCurrent_dq0[3];//dq0�����µ�ŵ�ٵ�Ч����������ֵ

	//��ʱ����
	//double temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18;
};

#endif