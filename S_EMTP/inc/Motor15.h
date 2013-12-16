#ifndef MOTOR15_H
#define MOTOR15_H

#include "Component.h"

class Motor15:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	Motor15(int firstNode,int lastNode);//���캯��,����Ϊ15�����Ľڵ���
	~Motor15();

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
	void calculateCoefficientMatrix_1();
	//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
	//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
	virtual void calculateDq0Results(double time);
	//����е���̣�����ת��ת��
	//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
	virtual void calculateOmega();
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &P,double Angle);
	virtual void invParkTransMatrix(TMatrixD &INRP,double Angle);

	void parkTransMatrix(double P[9][5],double Angle);
	void invParkTransMatrix(double invP[15][3],double Angle);
	//�ж��Ƿ��нڵ�ӵ�
	virtual bool isSeries(int from,int to);

private:
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//dq0����ϵ�µĵ������
	double Rp,Rr1,Rrt;
	//dq0����ϵ�µĵ翹����
	double Xl1,Xl2,Xlr1,Xm1,Xlt1,Xlt2,Xlrt,Xmt,Xdqm1,Xdqmt,Xs0;
	//�м����
	double Xs1,Xsm1,Xr1,Xst,Xsmt,Xrt;

	//���Ӳ��ѹ�������ͽ�Ƶ�ʵĻ�ֵ
	double Vbase,Ibase,wbase;

	//��е�����
	double Je;//ת������
	double p;//������

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/

	//�ɿ˱任�����������
	TMatrixD P;//(9,15)
	TMatrixD invP;//(15,9)

	double P_1[9][5];
	double invP_1[15][3];

	/*     ���λ��ַ���Ҫʹ�õ�ϵ������      */
	double coff;//alpha
	double T;//T=1.0/(deltaT*wbase)
	//
	TMatrixD A;//(11,11)
	TMatrixD B;//(11,11)
	TMatrixD AA;//(11,11)
	TMatrixD BB;//(11,11)

	double B_1[6][6];

	//
	TMatrixD Rss;//AA(1:9,1:9)
	TMatrixD Rsr;//AA(1:9,10:11)
	TMatrixD Rrs;//AA(10:11,1:9)
	TMatrixD Rrr;//AA(10:11,10:11)
	//
	TMatrixD Bss;//BB(1:9,1:9)
	TMatrixD Bsr;//BB(1:9,10:11)
	TMatrixD Brs;//BB(10:11,1:9)
	TMatrixD Brr;//BB(10:11,10:11)

	//����ϵ������ʹ��ת����abcde����ϵʱ�ܽ���
	TMatrixD Rs;//
	TMatrixD Rs1;//
	TMatrixD Rs2;//Rs2=Rs-Rs1;

	//
	TMatrixD Gs1;//
	TMatrixD Gs2;//
	TMatrixD Gs3;//
	TMatrixD Gs4;//

	double Gs_1[9][9];
	double Gs_2[9][9];
	double Gs_3[9][9];
	double Gs_4[9][2];

	//
	TMatrixD Gr1;//
	TMatrixD Gr2;//
	TMatrixD Gr3;//
	TMatrixD Gr4;//

	double Gr_1[2][2];
	double Gr_2[2][9];
	double Gr_3[2][9];
	double Gr_4[2][2];

	/*    ����ŷ��������Ҫʹ�õ�ϵ������     */
	//����

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vnode[30];//�ڵ��ѹ
	double Vabcde[15];//֧·��ѹ
	double Vabcde_1[15];
	double Vabcde_2[15];
	double Iabcde[15];//֧·����
	double Iabcde_1[15];
	double Iabcde_2[15];
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	double Vdqs[9];
	double Idqs[9];
	double Idqs_his[9];//Ԥ����
	double Idqs_his2[9];//Ԥ����
	double Idqs_forcast[9];//Ԥ��ֵ
	//dq0����ϵ�¶��Ӳ�Ĵ���
	TVectorD psi;
	TVectorD Isr;//�������ʱ�õĵ���

	double psi_1[6];
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
	double nortonEquivalentCurrent[15];//abcde�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[15];
	double nortonEquivalentCurrent_2[15];
	double nortonEquivalentCurrent_dq0[9];//dq0�����µ�ŵ�ٵ�Ч����������ֵ

	//��ʱ����
	double temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18;
};

#endif