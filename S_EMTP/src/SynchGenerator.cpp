#include "SynchGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
SynchGenerator::SynchGenerator(int a,int b,int c)
{
	isUserDef = 1;
	need_NEC=1;
	//���캯����a��b��cΪ�ڵ���
	nPort=3;
	nodeNumber=new int[3];
	type=6;

	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	w = 2*PI*60;//����ٶȣ�����ֵ,f=60Hz
	Ef = 1.3255;//���ŵ�ѹ������ֵ
	//��������й����޹�������ֵ
	P0 = 126391.44e6;//MW
	Q0 = 238242.25e6;//MVar
	Sm = 0;//���ڹ��ʣ�����ֵ���ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Vm = 1.0;//���ѹ��ֵ������ֵ
	Im = 0;//�������ֵ������ֵ���ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	ph = 1.642;//A���ѹ�ĳ�ʼ��λ�����ȱ�ʾ
	Sangle = 0;//�������ؽǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Vangle = 0;//A���ѹ��ʼ��ǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Iangle = 0;//A�������ʼ��ǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Angle = 0;//ת�ӽǵĳ�ʼֵ���ɼ���õ����������û������

	//dq0����ϵ�µĵ������,����ֵ
	Ra = 0.00621;
	Rf = 0.00025765;
	RD = 0.003033;
	RQ = 0.0030469;
	//dq0����ϵ�µĵ翹����,����ֵ
	Xad = 0.34288;
	Xaq = 0.15453;
	Xd = 0.40075;
	Xq = 0.2124;
	Xf = 0.37918;
	XD = 0.35411;
	XQ = 0.162953;
	XfD = Xad;
	Xaf = Xad;
	XaD = Xad;
	XaQ = Xaq;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	Vbase = 277e3*sqrt(2.0);//V
	Ibase = 366e3*sqrt(2.0);//A
	Zbase = Vbase/Ibase;

	//��е�����
	H = 1.0;
	D = 0;
	Tm = 0.4221;//����ת��

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//����coff��T
	coff = 99.0/101.0;
	T = 1.0/deltaT/w*0.5*(1+coff);

	//���ŵ�ѹ��ϵ��
	Vf=Rf/Xad;

	//�ɿ˱任�����������
	P.ResizeTo(1,2,1,3);
	invP.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,5,1,5);
	BB.ResizeTo(1,5,1,5);
	////
	Rss.ResizeTo(1,2,1,2);
	Rsr.ResizeTo(1,2,1,3);
	Rrs.ResizeTo(1,3,1,2);
	Rrr.ResizeTo(1,3,1,3);

	////
	Bss.ResizeTo(1,2,1,2);
	Bsr.ResizeTo(1,2,1,3);
	Brs.ResizeTo(1,3,1,2);
	Brr.ResizeTo(1,3,1,3);

	////����ϵ������ʹ��ת����abc����ϵʱ�ܽ���
	Rs.ResizeTo(1,2,1,2);
	R_ave.ResizeTo(1,2,1,2);
	R_res.ResizeTo(1,2,1,2);
	inv_Rs.ResizeTo(1,2,1,2);
	inv_Rr.ResizeTo(1,3,1,3);

	////
	Gs1.ResizeTo(1,2,1,2);
	Gs2.ResizeTo(1,2,1,3);
	Gs3.ResizeTo(1,2,1,2);
	Gs4.ResizeTo(1,2,1,3);
	Gs5.ResizeTo(1,2,1,3);
	Gs6.ResizeTo(1,2,1,2);
	Gs13.ResizeTo(1,2,1,2);
	Gs25.ResizeTo(1,2,1,3);

	////
	Gr1.ResizeTo(1,3,1,3);
	Gr2.ResizeTo(1,3,1,2);
	Gr3.ResizeTo(1,3,1,2);
	Gr4.ResizeTo(1,3,1,3);
	Gr5.ResizeTo(1,3,1,3);
	Gr15.ResizeTo(1,3,1,3);

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ�����������0����ʼ��ʱ������ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vabc[i-1] = 0;
		Vabc_1[i-1] = 0;
		Vabc_2[i-1] = 0;
		Iabc[i-1] = 0;
		Iabc_1[i-1] = 0;
		Iabc_2[i-1] = 0;
	}
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=2;i++)
	{
		Vdq0[i-1] = 0;
		Idq0[i-1] = 0;
		Idq0_his[i-1] = 0;
		Idq0_his2[i-1] = 0;
		Idq0_forcast[i-1] = 0;
	}
	//dq0����ϵ��ת�Ӳ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		VfDQ[i-1] = 0;
		IfDQ[i-1] = 0;
		IfDQ_his[i-1] = 0; 
	}

	/*  ��е��������  */
	wr = 1;//ת��ת��
	wr_his = 1;//wr����ʷֵ
	curAngle = 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_his = 0;//curAngle����ʷֵ,û����

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���

	for(int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
	}

	/************************************************************************/
	/*              ��EMTP�ӿ�����ı���                                    */
	/************************************************************************/
	nodeNumber[0] = a;
	nodeNumber[1] = b;
	nodeNumber[2] = c;

	/************************************************************************/
	/*              ����ϵ�������ŵ�ȵ�ֵ����                              */
	/************************************************************************/
	calculateCoefficientMatrix();
}
SynchGenerator::~SynchGenerator()
{//��������
	
}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/
void SynchGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time)
{//��ʼ��֧·��ѹ����
	//��ʼ��wr��curAngle
	Sm = sqrt(P0*P0+Q0*Q0);
	Sangle = atan(Q0/P0);
	Sm = Sm/(1.5*Vbase*Ibase);
	Im = Sm/Vm;
	Vangle = ph/180*PI;
	Iangle = Vangle - Sangle;
	curAngle = atan((Vm*sin(Vangle)+Xq*Im*cos(Iangle)+Ra*Im*sin(Iangle))/
		(Vm*cos(Vangle)-Xq*Im*sin(Iangle)+Ra*Im*cos(Iangle)))
		-0.5*PI;

	Iabc[0] = Im*cos(Iangle);
	Iabc[1] = Im*cos(Iangle-(2.0/3.0)*PI);
	Iabc[2] = Im*cos(Iangle+(2.0/3.0)*PI);

	Vabc[0] = Vm*cos(Vangle);
	Vabc[1] = Vm*cos(Vangle-(2.0/3.0)*PI);
	Vabc[2] = Vm*cos(Vangle+(2.0/3.0)*PI);

	parkTransMatrix(P,curAngle);

	for (int i=0;i<2;i++)
	{
		Vdq0[i] = P(i+1,1)*Vabc[0]+P(i+1,2)*Vabc[1]+P(i+1,3)*Vabc[2];
		Vdq0_his[i] = Vdq0[i];
		Vdq0_his2[i] = Vdq0[i];
		Vdq0_emend[i] = Vdq0[i];
		Idq0[i] = P(i+1,1)*Iabc[0]+P(i+1,2)*Iabc[1]+P(i+1,3)*Iabc[2];
		Idq0_his[i] = Idq0[i];
		Idq0_his2[i] = Idq0[i];
		Idq0_forcast[i] = Idq0[i];
		Idq0_emend[i] = Idq0[i];
	}

	IfDQ[0] = (Vdq0[1]+Ra*Idq0[1]+Xd*Idq0[0])/Xad;
	IfDQ[1] = 0;
	IfDQ[2] = 0;

	VfDQ[0] = IfDQ[0]*Rf;
	VfDQ[1] = 0;
	VfDQ[2] = 0;

	//���е���̣���ʼ��wr��curAngle
	calculateOmega();

	//ŵ�ٵ�ֵ������ʼ��
	//����dq0�����µ���ʷ������
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			- Gs1(i,1)*Idq0_forcast[0]-Gs1(i,2)*Idq0_forcast[1]
		- Gs3(i,1)*Idq0[0]-Gs3(i,2)*Idq0[1]
		+ Gs4(i,1)*IfDQ[0]+Gs4(i,2)*IfDQ[1]+Gs4(i,3)*IfDQ[2]
		+ Gs6(i,1)*Vdq0[0]+Gs6(i,2)*Vdq0[1]
		+ Gs25(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}

	//����park��任����
	invParkTransMatrix(invP,curAngle);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invP(i,1)*nortonEquivalentCurrent_dq0[0]
		+ invP(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}

	//��PSCAD���ݳ�ʼ�����ӵ�ѹ����
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<3;i++)
	{
		Iabc[i] = initialCurrentArray[ptr+i];
	}
	ptr+=3;

	calculateCoefficientMatrix();
}
void SynchGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//�ӽڵ��ѹ�����ж�֧·�ڵ��ѹ

	for (int i=0;i<3;i++)
	{
		Vabc_1[i] =Vabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vabc[i]=0;
		}
		else
		{
			Vabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}
}
void SynchGenerator::calculateBranchVoltage()
{//����֧·��ѹ
}

void SynchGenerator::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Iabc_1[i] = Iabc[i];
		Iabc[i] = -Vabc[i]*nortonEquivalentConductance-nortonEquivalentCurrent[i];
	}
}

void SynchGenerator::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	calculateCoefficientMatrix();
	calculateDq0Results();//����dq0����ϵ�µĽ����Ԥ�����ӵ���
	calculateOmega();
	
	//����dq0�����µ���ʷ������
	for (int i=1;i<=2;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			- Gs1(i,1)*Idq0_forcast[0]-Gs1(i,2)*Idq0_forcast[1]
			- Gs3(i,1)*Idq0[0]-Gs3(i,2)*Idq0[1]
			+ Gs4(i,1)*IfDQ[0]+Gs4(i,2)*IfDQ[1]+Gs4(i,3)*IfDQ[2]
			+ Gs6(i,1)*Vdq0[0]+Gs6(i,2)*Vdq0[1]
			+ Gs25(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}

	//����park��任����
	invParkTransMatrix(invP,curAngle);

	for (int k=0;k<3;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent[i-1] = invP(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invP(i,2)*nortonEquivalentCurrent_dq0[1];
		nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
	}
}

void SynchGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}
}
void SynchGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance;
	}
}

void SynchGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix(counter,ptr) = Iabc[0];
	branchCurrentMatrix(counter,ptr+1) = Iabc[1];
	branchCurrentMatrix(counter,ptr+2) = Iabc[2];
	ptr+=3;
}

void SynchGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix_1[counter-1][ptr] = Iabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Iabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Iabc[2];
	ptr+=3;
}

void SynchGenerator::interpolate(double ratio)
{//������ֵ����֧·�ĵ�ѹ�������в�ֵ

	//��ѹ����������ֵ
	for (int k=0;k<3;k++)
	{
		Vabc[k] = (1-ratio)*Vabc_1[k] + ratio*Vabc[k];
		Iabc[k] = (1-ratio)*Iabc_1[k] + ratio*Iabc[k];
	}

	//ŵ�ٵ�ֵ������ֵ
	double temp1,temp2;
	for (int k=0;k<3;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void SynchGenerator::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<3;i++)
		{
			Vabc_2[i] = Vabc_1[i];
			Iabc_2[i] = Iabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<3;i++)
		{
			Vabc_1[i] = Vabc_2[i];
			Iabc_1[i] = Iabc_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<3;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

void SynchGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
{//ת��ת�ٱ仯ʱ���µ��ɾ���
	int N;
	double tmpd;
	for (int i=0;i<3;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his+this->nortonEquivalentConductance;
	}
}

/************************************************************************/
/*                ������ڲ����õĺ���                                  */
/************************************************************************/
//����ϵ������AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
//ͬʱ����nortonEquivalentConductance
void SynchGenerator::calculateCoefficientMatrix()
{
	//����AA���󣬼��������
	AA(1,1) = Ra+2.0*T*Xd;
	AA(1,2) = -wr*Xq;
	AA(1,3) = 2.0*T*Xaf;
	AA(1,4) = 2.0*T*XaD;
	AA(1,5) = -wr*XaQ;

	AA(2,1) = wr*Xd;
	AA(2,2) = Ra+2.0*T*Xq;
	AA(2,3) = wr*Xaf;
	AA(2,4) = wr*XaD;
	AA(2,5) = 2.0*T*XaQ;

	AA(3,1) = 2.0*T*Xaf;
	AA(3,2) = 0;
	AA(3,3) = Rf+2.0*T*Xf;
	AA(3,4) = 2.0*T*XfD;
	AA(3,5) = 0;

	AA(4,1) = 2.0*T*XaD;
	AA(4,2) = 0;
	AA(4,3) = 2.0*T*XfD;
	AA(4,4) = RD+2.0*T*XD;
	AA(4,5) = 0;

	AA(5,1) = 0;
	AA(5,2) = 2.0*T*XaQ;
	AA(5,3) = 0;
	AA(5,4) = 0;
	AA(5,5) = RQ+2.0*T*XQ;

	//����BB����
	BB(1,1) = coff*Ra-2.0*T*Xd;
	BB(1,2) = -coff*wr*Xq;
	BB(1,3) = -2.0*T*Xaf;
	BB(1,4) = -2.0*T*XaD;
	BB(1,5) = -coff*wr*XaQ;

	BB(2,1) = coff*wr*Xd;
	BB(2,2) = coff*Ra-2.0*T*Xq;
	BB(2,3) = coff*wr*Xaf;
	BB(2,4) = coff*wr*XaD;
	BB(2,5) = -2.0*T*XaQ;

	BB(3,1) = -2.0*T*Xaf;
	BB(3,2) = 0;
	BB(3,3) = coff*Rf-2.0*T*Xf;
	BB(3,4) = -2.0*T*XfD;
	BB(3,5) = 0;

	BB(4,1) = -2.0*T*XaD;
	BB(4,2) = 0;
	BB(4,3) = -2.0*T*XfD;
	BB(4,4) = coff*RD-2.0*T*XD;
	BB(4,5) = 0;

	BB(5,1) = 0;
	BB(5,2) = -2.0*T*XaQ;
	BB(5,3) = 0;
	BB(5,4) = 0;
	BB(5,5) = coff*RQ-2.0*T*XQ;

	//�������Rss,Rsr,Rrs,Rrr,Bss,Bsr,Brs,Brr
	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=2;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rsr(i,j) = AA(i,j+2);
			Bsr(i,j) = BB(i,j+2);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=2;j++)
		{
			Rrs(i,j) = AA(i+2,j);
			Brs(i,j) = BB(i+2,j);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rrr(i,j) = AA(i+2,j+2);
			Brr(i,j) = BB(i+2,j+2);
		}
	}

	//����Rs,R_ave,R_res,inv_Rs,inv_Rr
	inv_Rr = Rrr;
	inv_Rr.Invert();

	Rs = Rss-Rsr*inv_Rr*Rrs;

	double resistance = (Rs(1,1)+Rs(2,2))/2.0;
	R_ave(1,1) = resistance;
	R_ave(1,2) = 0;
	R_ave(2,1) = 0;
	R_ave(2,2) = resistance;
	
	R_res = Rs - R_ave;

	//����ŵ�ٵ�Ч�絼
	nortonEquivalentConductance = 1/resistance;
	inv_Rs = R_ave;
	inv_Rs.Invert();
	nortonEquivalentConductance /= Zbase;

	//����Gs1...Gs6
	Gs1 = inv_Rs*R_res;
	Gs1 = (-1.0)*Gs1;
	Gs2 = inv_Rs*Rsr*inv_Rr;
	Gs2 = (-1.0)*Gs2;
	Gs3 = Gs2*Brs + inv_Rs*Bss;
	Gs3 = (-1.0)*Gs3;
	Gs4 = Gs2*Brr + inv_Rs*Bsr;
	Gs4 = (-1.0)*Gs4;
	Gs5 = coff*Gs2;
	Gs6 = coff*inv_Rs;
	Gs13 = Gs1+Gs3;
	Gs25 = Gs2+Gs5;

	//����Gr1...Gr5
	Gr1 = inv_Rr;
	Gr2 = inv_Rr*Rrs;
	Gr2 = (-1.0)*Gr2;
	Gr3 = inv_Rr*Brs;
	Gr3 = (-1.0)*Gr3;
	Gr4 = inv_Rr*Brr;
	Gr4 = (-1.0)*Gr4;
	Gr5 = coff*inv_Rr;
	Gr15 = Gr1+Gr5;
}

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
void SynchGenerator::calculateDq0Results()
{
	double Vabc_pu[3];
	double Iabc_pu[3];

	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<3;i++)
	{
		Vabc_pu[i] = Vabc[i]/Vbase;
		Iabc_pu[i] = Iabc[i]/Ibase;
	}

	//����park�任����
	parkTransMatrix(P,curAngle);

	//���㶨�ӵ�ѹ
	for (int i=0;i<2;i++)
	{
		Vdq0_his2[i] = Vdq0_his[i];
		Vdq0_his[i] = Vdq0[i];
		Vdq0[i] = P(i+1,1)*Vabc_pu[0]+P(i+1,2)*Vabc_pu[1]+P(i+1,3)*Vabc_pu[2];
	}

	//���㶨�ӵ���
	for (int i=0;i<2;i++)
	{
		Idq0_his2[i] = Idq0_his[i];
		Idq0_his[i] = Idq0[i];
		Idq0[i] = P(i+1,1)*Iabc_pu[0]+P(i+1,2)*Iabc_pu[1]+P(i+1,3)*Iabc_pu[2];
	}
	
	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<2;i++)
	{
		Idq0_forcast[i] = 1.3333333*Idq0[i]+0.3333333*Idq0_his[i]-0.6666666*Idq0_his2[i];
	}

	//����ת�ӵ�ѹ
	Vf = Rf/Xad;
	VfDQ[0] = Ef*Vf;
	VfDQ[1] = 0;
	VfDQ[2] = 0;

	//����ת�ӵ���
	for (int i=0;i<3;i++)
	{
		IfDQ_his[i] = IfDQ[i];
	}
	for (int i=1;i<=3;i++)
	{
		IfDQ[i-1] =- Gr2(i,1)*Idq0[0] - Gr2(i,2)*Idq0[1]
		- Gr3(i,1)*Idq0_his[0] - Gr3(i,2)*Idq0_his[1]
		+ Gr4(i,1)*IfDQ_his[0] + Gr4(i,2)*IfDQ_his[1] + Gr4(i,3)*IfDQ_his[2]
		+ Gr15(i,1)*VfDQ[0];//VfDQ[1]=VfDQ[2]=0
	}
}

//����е���̣�����ת��ת��
//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
void SynchGenerator::calculateOmega()
{
	double tt1,tt2,tt3;
	tt1 = (Xd-Xq)*Idq0[0]*Idq0[1];
	tt2 = -Xad*(IfDQ[0]*Idq0[1]+IfDQ[1]*Idq0[1]);
	tt3 = -Xaq*IfDQ[2]*Idq0[0];

	wr_his = wr;
	wr = (Tm+(tt1+tt2-tt3))*deltaT/2/H+wr_his;

	curAngle = curAngle + deltaT*wr*w/2 + deltaT*wr_his*w/2;
	if (curAngle>2*PI)
	{
		curAngle -= 2*PI;
	}
}

//�����ɿ˱任�뷴�任����
void SynchGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void SynchGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}
