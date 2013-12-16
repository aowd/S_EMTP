#ifndef WINDINDUCGENERATOR_H
#define WINDINDUCGENERATOR_H

#include "Component.h"

class WindInducGenerator:public Component {
public:
	/************************************************************************/
	/*                ���캯������������                                    */
	/************************************************************************/
	WindInducGenerator(int firstNode,int lastNode);//���캯��,a��b��cΪ�ڵ���
	~WindInducGenerator();

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
	//����ϵ������
	virtual void calculateCoefficientMatrix();
	//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̵���
	//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
	virtual void calculateDq0Results();
	//����е���̣�����ת��ת��
	//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
	virtual void calculateOmega(double time);
	virtual void calculateTw(double time);
	//�����ɿ˱任�뷴�任����
	virtual void parkTransMatrix(TMatrixD &T,double angle);
	virtual void invParkTransMatrix(TMatrixD &T_,double angle);
public:
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	double w;//����ٶȣ�����ֵ
	double w_pu;//����ٶȣ�����ֵ
	//*****************
	//��ν��г�ʼ��
	//*****************

	//dq0����ϵ�µĵ������
	double Rs,Rr;
	//dq0����ϵ�µĵ翹����
	double Xm,Xls,Xlr;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	double Vbase,Ibase,Zbase,Sbase;//����Zbase���Լ���õ�

	//��ת��������
	double turnratio;

	//control=1,ת�ٿ���
	//control=0,ת�ؿ���
	int control;

	//��е�����
	double H_g;
	double D_g;
	double H_wt;
	double D_wt;
	double k;

	//�������
	double v;
	double p_air;
	double R;
	double gr;
	double Pw;
	double Tw;
	double Cp;


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

	//����ϵ������ʹ��ת����abc����ϵʱ�ܽ���
	TMatrixD AA_ave;//
	TMatrixD AA_res;//AA_res=AA-AAave;
	TMatrixD inv_Aa;//inv(AA_ave)
	
	//
	TMatrixD Gn1;//
	TMatrixD Gn2;//
	TMatrixD Gn3;//
	

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
	double curAngle_s_his;
	double curAngle_r;//d������ת��a��ĽǶ�
	double curAngle_r_his;
	double w_wt;
	double w_wt_his;
	double curAngle_wt;
	double curAngle_wt_his;
	double pitch_angle;

	TMatrixD Co;
	TMatrixD Co_inv;
	TMatrixD bo;
	TMatrixD co;

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	double nortonEquivalentConductance[2];
	double nortonEquivalentConductance_his[2];//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	double nortonEquivalentCurrent[6];//abc�����µ�ŵ�ٵ�Ч����������ֵ
	double nortonEquivalentCurrent_1[6];
	double nortonEquivalentCurrent_2[6];
	double nortonEquivalentCurrent_dq0[4];//dq0�����µ�ŵ�ٵ�Ч����������ֵ
};

#endif