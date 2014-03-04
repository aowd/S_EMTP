#include "WoundInducGenerator2.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
WoundInducGenerator2::WoundInducGenerator2(int id, int firstNode,int lastNode, int control,double Vw)
{
	this->id = id;
	isUserDef = 0;
	need_NEC=1;
	//���캯����a��b��cΪ�ڵ���
	nPort=6;
	nodeNumber=new int[nPort];
	type=29;

	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	Wbase = 2*PI*50;//����ٶȣ�����ֵ,f=60Hz
	w=1;

	//dq0����ϵ�µĵ������,����ֵ
	Rs = 0.048;
	Rr = 0.018;
	//dq0����ϵ�µĵ翹����,����ֵ
	Xm = 3.8;
	Xls = 0.075;
	Xlr = 0.12;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	Vbase=690;
	Sbase=2e6;
	Vbase = Vbase/sqrt(3.0);//V
	Ibase = Sbase/3/Vbase;//A
	Zbase = Vbase/Ibase;

	//��ת������������
	turnratio=0.3;

	this->control = control;

	//��е�����
	H = 0.1;
	D = 0;
	TL =-0.7;//����ת��

	//�������
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=Vw-0.1*(id-1);

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/

	//�ɿ˱任�����������
	Ps.ResizeTo(1,3,1,3);
	invPs.ResizeTo(1,3,1,3);
	Pr.ResizeTo(1,3,1,3);
	invPr.ResizeTo(1,3,1,3);
	P.ResizeTo(1,6,1,6);
	invP.ResizeTo(1,6,1,6);

	//
	AA.ResizeTo(1,6,1,6);
	BB.ResizeTo(1,6,1,6);
	AA_inv.ResizeTo(1,6,1,6);
	//
	nortonEquivalentConductanceMatrix.ResizeTo(1,6,1,6);

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ�����������0����ʼ��ʱ������ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vsabc[i-1] = 0;
		Vsabc_1[i-1] = 0;
		Vsabc_2[i-1] = 0;
		Isabc[i-1] = 0;
		Isabc_1[i-1] = 0;
		Isabc_2[i-1] = 0;
	}
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vsdq0[i-1] = 0;
		Isdq0[i-1] = 0;
	}

	//abc����ϵ��ת����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vrabc[i-1] = 0;
		Vrabc_1[i-1] = 0;
		Vrabc_2[i-1] = 0;
		Irabc[i-1] = 0;
		Irabc_1[i-1] = 0;
		Irabc_2[i-1] = 0;
	}
	//dq0����ϵ��ת����Ȧ�ĵ�ѹ����
	for (int i=1;i<=3;i++)
	{
		Vrdq0[i-1] = 0;
		Irdq0[i-1] = 0;
	}

	/*  ��е��������  */
	if (control==1)
	{
		wr = 1.15;//ת��ת��
		wr_his = 1.15;//wr����ʷֵ
	} 
	else
	{
		wr = 0;//ת��ת��
		wr_his = 0;//wr����ʷֵ
	}
	curAngle_s= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_r= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Te=0;
	Te_his=0;


	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	for (int i=0;i<36;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_pu[i]=0;
	}

	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_dq0_his[i-1] = 0;
	}

	/************************************************************************/
	/*              ��EMTP�ӿ�����ı���                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}
#ifdef WIND_VELOCITY_DATA_INPUT
	WindVelocityCounter=0;
	WindVelocityCounter_bak=0;
#endif
}
WoundInducGenerator2::~WoundInducGenerator2()
{//��������

}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/

void WoundInducGenerator2::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time){ptr+=6;}

void WoundInducGenerator2::readNodeVoltage(TVectorD &nodeVoltageArray)
{//�ӽڵ��ѹ�����ж�֧·�ڵ��ѹ

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
		Vsabc_his2[i] = Vsabc_his[i];
		Vsabc_his[i] = Vsabc[i];
	}

	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i]==0)
		{
			Vsabc[i]=0;
		}
		else
		{
			Vsabc[i]= nodeVoltageArray(nodeNumber[i]);
		}	
	}

	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
		Vrabc_his2[i] = Vrabc_his[i];
		Vrabc_his[i] = Vrabc[i];
	}


	for (int i=0;i<3;i++)
	{
		if (nodeNumber[i+3]==0)
		{
			Vrabc[i]=0;
		}
		else
		{
			Vrabc[i]= nodeVoltageArray(nodeNumber[i+3]);
		}	
	}
}
void WoundInducGenerator2::calculateBranchVoltage()
{//����֧·��ѹ
}

void WoundInducGenerator2::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i]=0;
		for (int j=0;j<3;j++)
		{
			Isabc[i]+=Vsabc[j]*nortonEquivalentConductance[6*i+j];
			Isabc[i]+=Vrabc[j]*nortonEquivalentConductance[6*i+j+3];
		}
		Isabc[i]+=nortonEquivalentCurrent[i];

		Irabc_1[i] = Irabc[i];
		Irabc[i]=0;
		for (int j=0;j<3;j++)
		{
			Irabc[i]+=Vsabc[j]*nortonEquivalentConductance[6*3+6*i+j];
			Irabc[i]+=Vrabc[j]*nortonEquivalentConductance[6*3+6*i+j+3];
		}
		Irabc[i]+=nortonEquivalentCurrent[i+3];
	}

	calculateDq0Results();
}

void WoundInducGenerator2::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����

	//calculateDq0Results();//����dq0����ϵ�µĽ��
	//calculateOmega();
	//calculateCoefficientMatrix();

	TMatrixD CC=AA_inv*BB;
	//����dq0�����µ���ʷ������
	for (int i=1;i<=3;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];

		nortonEquivalentCurrent_dq0[i-1] =-(CC(i,1)*Isdq0[0]+CC(i,2)*Isdq0[1]+CC(i,3)*Isdq0[2]+CC(i,4)*Irdq0[0]+CC(i,5)*Irdq0[1]+CC(i,6)*Irdq0[2])
			+ AA_inv(i,1)*Vsdq0[0]+AA_inv(i,2)*Vsdq0[1]+AA_inv(i,3)*Vsdq0[2]+AA_inv(i,4)*Vrdq0[0]+AA_inv(i,5)*Vrdq0[1]+AA_inv(i,6)*Vrdq0[2];


		nortonEquivalentCurrent_dq0_his[i+2] = nortonEquivalentCurrent_dq0[i+2];

		nortonEquivalentCurrent_dq0[i+2] =-(CC(i+3,1)*Isdq0[0]+CC(i+3,2)*Isdq0[1]+CC(i+3,3)*Isdq0[2]+CC(i+3,4)*Irdq0[0]+CC(i+3,5)*Irdq0[1]+CC(i+3,6)*Irdq0[2])
			+ AA_inv(i+3,1)*Vsdq0[0]+AA_inv(i+3,2)*Vsdq0[1]+AA_inv(i+3,3)*Vsdq0[2]+AA_inv(i+3,4)*Vrdq0[0]+AA_inv(i+3,5)*Vrdq0[1]+AA_inv(i+3,6)*Vrdq0[2];
	}

	//����park��任����
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]+ invPs(i,2)*nortonEquivalentCurrent_dq0[1]+invPs(i,3)*nortonEquivalentCurrent_dq0[2];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[3]+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[4]+invPr(i-3,3)*nortonEquivalentCurrent_dq0[5];			
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}		
	}
}

void WoundInducGenerator2::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) = nodeNortonEquivalentCurrentArray(N)- nortonEquivalentCurrent[i];
	}
}
void WoundInducGenerator2::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N,M;
	for (int i=0;i<nPort;i++)
	{
		N = nodeNumber[i];
		for (int j=0;j<nPort;j++)
		{
			M=nodeNumber[j];
			conductanceMatrix(N,M)+=this->nortonEquivalentConductance[nPort*i+j];
		}
	}
}

void WoundInducGenerator2::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix(counter,ptr) = Isabc[0];
	branchCurrentMatrix(counter,ptr+1) = Isabc[1];
	branchCurrentMatrix(counter,ptr+2) = Isabc[2];
	ptr+=3;
	branchCurrentMatrix(counter,ptr) = Irabc[0];
	branchCurrentMatrix(counter,ptr+1) = Irabc[1];
	branchCurrentMatrix(counter,ptr+2) = Irabc[2];
	ptr+=3;
}

void WoundInducGenerator2::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//����֧·����
	branchCurrentMatrix_1[counter-1][ptr] = Isabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Isabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Isabc[2];
	ptr+=3;
	branchCurrentMatrix_1[counter-1][ptr] = Irabc[0];
	branchCurrentMatrix_1[counter-1][ptr+1] = Irabc[1];
	branchCurrentMatrix_1[counter-1][ptr+2] = Irabc[2];
	ptr+=3;
}
void WoundInducGenerator2::interpolate(double ratio)
{//������ֵ����֧·�ĵ�ѹ�������в�ֵ

	//��ѹ����������ֵ
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
	}

	for (int k=0;k<3;k++)
	{
		Vsdq0[k] = (1-ratio)*Vsdq0_his[k] + ratio*Vsdq0[k];
		Isdq0[k] = (1-ratio)*Isdq0_his[k] + ratio*Isdq0[k];
		Vrdq0[k] = (1-ratio)*Vrdq0_his[k] + ratio*Vrdq0[k];
		Irdq0[k] = (1-ratio)*Irdq0_his[k] + ratio*Irdq0[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_dq0[k] = (1-ratio)*nortonEquivalentCurrent_dq0_his[k] + ratio*nortonEquivalentCurrent_dq0[k];
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent[k] = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	}

	curAngle_r = (1-ratio)*curAngle_r_1 + ratio*curAngle_r;
	curAngle_s = (1-ratio)*curAngle_s_1 + ratio*curAngle_s;
	wr = (1-ratio)*wr_his + ratio*wr;

	//ŵ�ٵ�ֵ������ֵ
	//double temp1,temp2;
	//for (int k=0;k<6;k++)
	//{
	//	temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
	//	temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

	//	nortonEquivalentCurrent_1[k] = temp1;
	//	nortonEquivalentCurrent[k] = temp2;
	//}
}

void WoundInducGenerator2::updateResult(int updateTypeNum)
{
	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<3;i++)
		{
			Vsabc_2[i] = Vsabc_1[i];
			Isabc_2[i] = Isabc_1[i];
			Vrabc_2[i] = Vrabc_1[i];
			Irabc_2[i] = Irabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_2[i+3] = nortonEquivalentCurrent_1[i+3];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<3;i++)
		{
			Vsabc_1[i] = Vsabc_2[i];
			Isabc_1[i] = Isabc_2[i];
			Vrabc_1[i] = Vrabc_2[i];
			Irabc_1[i] = Irabc_2[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
			nortonEquivalentCurrent_1[i+3] = nortonEquivalentCurrent_2[i+3];
		}
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

/************************************************************************/
/*                ������ڲ����õĺ���                                  */
/************************************************************************/
//����ϵ������AA,BB
//ͬʱ����nortonEquivalentConductance
void WoundInducGenerator2::calculateCoefficientMatrix()
{
	//����AA���󣬼��������
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w*(Xm+Xls);
	AA(1,4) = T*Xm;
	AA(1,5) = -w*Xm;

	AA(2,1) = w*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,4) = w*Xm;
	AA(2,5) = T*Xm;

	AA(3,3)=Rs+T*(Xm+Xls);
	AA(3,6)=T*Xm;

	AA(4,1) = T*Xm;
	AA(4,2) = (wr-w)*Xm;
	AA(4,4) = Rr+T*(Xm+Xlr);
	AA(4,5) = (wr-w)*(Xm+Xlr);

	AA(5,1) = (w-wr)*Xm;
	AA(5,2) = T*Xm;
	AA(5,4) = (w-wr)*(Xm+Xlr);
	AA(5,5) = Rr+T*(Xm+Xlr);

	AA(6,3)=T*Xm;
	AA(6,6)=Rr+T*(Xm+Xlr);

	//����BB����
	BB(1,1) = Rs-T*(Xm+Xls);
	BB(1,2) = -w*(Xm+Xls);
	BB(1,4) = -T*Xm;
	BB(1,5) = -w*Xm;

	BB(2,1) = w*(Xm+Xls);
	BB(2,2) = Rs-T*(Xm+Xls);
	BB(2,4) = w*Xm;
	BB(2,5) = -T*Xm;

	BB(3,3)=Rs-T*(Xm+Xls);
	BB(3,6)=-T*Xm;

	BB(4,1) = -T*Xm;
	BB(4,2) = (wr-w)*Xm;
	BB(4,4) = Rr-T*(Xm+Xlr);
	BB(4,5) = (wr-w)*(Xm+Xlr);

	BB(5,1) = (w-wr)*Xm;
	BB(5,2) = -T*Xm;
	BB(5,4) = (w-wr)*(Xm+Xlr);
	BB(5,5) = Rr-T*(Xm+Xlr);

	BB(6,3)=-T*Xm;
	BB(6,6)=Rr-T*(Xm+Xlr);

	AA_inv = AA;
	AA_inv.Invert();
}

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
void WoundInducGenerator2::calculateDq0Results()
{
	double Vsabc_pu[3];
	double Isabc_pu[3];
	double Vrabc_pu[3];
	double Irabc_pu[3];

	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<3;i++)
	{
		Vsabc_pu[i] = Vsabc[i]/Vbase;
		Isabc_pu[i] = Isabc[i]/Ibase;
		Vrabc_pu[i] = Vrabc[i]/Vbase*turnratio;
		Irabc_pu[i] = Irabc[i]/Ibase/turnratio;
	}

	//����park�任����
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	//�����ѹ
	for (int i=0;i<3;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//�������
	for (int i=0;i<3;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}
}
//����е���̣�����ת��ת��
//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
void WoundInducGenerator2::calculateOmega()
{
	//���÷������϶���ʽ
	if (control==2)
	{
		//���������ת��
		calculateTwind();
		//������ת��Ϊ�϶�ת��
		TL=-Twind;
	}
	wr_his = wr;
	curAngle_r_1= curAngle_r;
	curAngle_s_1=curAngle_s;

	if (control==0||control==2)
	{
		Te_his=Te;
		Te=0.5*Xm*(Irdq0[0]*Isdq0[1]-Irdq0[1]*Isdq0[0]);
		wr = ((Te+Te_his)/2-TL)*deltaT/2/H+wr_his;
		curAngle_r = curAngle_r + deltaT*Wbase*(w-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w-wr)*Wbase;
		curAngle_s = curAngle_s + deltaT*w*Wbase;
	}
}

//�����ɿ˱任�뷴�任����
void WoundInducGenerator2::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
	T(3,1)= 2.0/3.0*1.0/2.0;
	T(3,2)= 2.0/3.0*1.0/2.0;
	T(3,3)= 2.0/3.0*1.0/2.0;
}
void WoundInducGenerator2::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(1,3)= 1.0/2.0;
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(2,3)= 1.0/2.0;
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
	T_(3,3)= 1.0/2.0;
}


// ƽ��ģ����أ�xuyin��20121224
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
void WoundInducGenerator2::storeInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc_bak[k] = Vsabc[k];
		Vsabc_his_bak[k] = Vsabc_his[k];
		Vsabc_his2_bak[k] = Vsabc_his2[k];
		Isabc_bak[k] = Isabc[k];
		Vrabc_bak[k] = Vrabc[k];		
		Vrabc_his_bak[k] = Vrabc_his[k];
		Vrabc_his2_bak[k] = Vrabc_his2[k];
		Irabc_bak[k] = Irabc[k];
	}

	for (int k=0; k<3; k++) {
		Vsdq0_bak[k] = Vsdq0[k];
		Vsdq0_his_bak[k] = Vsdq0_his[k];
		Vsdq0_his2_bak[k] = Vsdq0_his2[k];
		Isdq0_bak[k] = Isdq0[k];
		Isdq0_his_bak[k] = Isdq0_his[k];
		Isdq0_his2_bak[k] = Isdq0_his2[k];
		Vrdq0_bak[k] = Vrdq0[k];
		Vrdq0_his_bak[k] = Vrdq0_his[k];
		Vrdq0_his2_bak[k] = Vrdq0_his2[k];
		Irdq0_bak[k] = Irdq0[k];
		Irdq0_his_bak[k] = Irdq0_his[k];
		Irdq0_his2_bak[k] = Irdq0_his2[k];
	}

	wr_bak = wr;
	wr_his_bak = wr_his;
	Te_bak = Te;
	Te_his_bak=Te_his;
	curAngle_s_bak = curAngle_s;
	curAngle_r_bak = curAngle_r;
	WindVelocityCounter_bak = WindVelocityCounter;
}

// ����У��ʱ�ָ��ڲ�����
void WoundInducGenerator2::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Vsabc_his[k] = Vsabc_his_bak[k];
		Vsabc_his2[k] = Vsabc_his2_bak[k];
		Isabc[k] = Isabc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Vrabc_his[k] = Vrabc_his_bak[k];
		Vrabc_his2[k] = Vrabc_his2_bak[k];
		Irabc[k] = Irabc_bak[k];
	}

	for (int k=0; k<3; k++) {
		Vsdq0[k] = Vsdq0_bak[k];
		Vsdq0_his[k] = Vsdq0_his_bak[k];
		Vsdq0_his2[k] = Vsdq0_his2_bak[k];
		Isdq0[k] = Isdq0_bak[k];
		Isdq0_his[k] = Isdq0_his_bak[k];
		Isdq0_his2[k] = Isdq0_his2_bak[k];
		Vrdq0[k] = Vrdq0_bak[k];
		Vrdq0_his[k] = Vrdq0_his_bak[k];
		Vrdq0_his2[k] = Vrdq0_his2_bak[k];
		Irdq0[k] = Irdq0_bak[k];
		Irdq0_his[k] = Irdq0_his_bak[k];
		Irdq0_his2[k] = Irdq0_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te = Te_bak;
	Te_his = Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
	WindVelocityCounter = WindVelocityCounter_bak;
}

void WoundInducGenerator2::calculateNortonEquivalentResistance(double time)
{
	T = 2.0/deltaT/Wbase;
	calculateOmega();
	calculateCoefficientMatrix();
	
	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<4;i++)
	{
		for (int j=1;j<4;j++)
		{
			P(i,j)=Ps(i,j);
			P(i+3,j+3)=Pr(i,j);
			invP(i,j)=invPs(i,j);
			invP(i+3,j+3)=invPr(i,j);
		}
	}

	nortonEquivalentConductanceMatrix=invP*AA_inv*P;

	for (int i=0;i<6;i++)
	{
		for (int j=0;j<6;j++)
		{
			nortonEquivalentConductance[6*i+j]=nortonEquivalentConductanceMatrix(i+1,j+1);
		}
	}

	//����Ч������ֵ���㵽����ֵ
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			nortonEquivalentConductance[6*i+j]=nortonEquivalentConductance[6*i+j]/Zbase;
			nortonEquivalentConductance[6*(i+3)+j]=nortonEquivalentConductance[6*(i+3)+j]/Zbase*turnratio;
			nortonEquivalentConductance[6*i+(j+3)]=nortonEquivalentConductance[6*i+(j+3)]/Zbase*turnratio;
			nortonEquivalentConductance[6*(i+3)+(j+3)]=nortonEquivalentConductance[6*(i+3)+(j+3)]/Zbase*turnratio*turnratio;
		}
	}
}

// �����ʼ��ר�ú���, xuyin, 20121226
void WoundInducGenerator2::initializeGen(double** GenInitialMatrix, int& ptr)
{
	// ���˵�һ��Ϊʱ��t�⣬ÿ�������ӦGenInitialMatrix��15��
	// 1-3�У�Vsabc��4-6�У�Vrabc��7-9�У�Isabc��10-12�У�Irabc��13-15�У�Tm��Theta��wr
	// GenInitialMatrix�����У��ֱ��Ӧ��ǰʱ���Լ�ǰ����ʱ����ֵ������ʱ��˳��

	// ��������ֵ*1000
	for (int i=0; i<3; i++)
		for (int j=0; j<12; j++)
			GenInitialMatrix[i][ptr+j] *= 1000;

	// ��ȡʱ��
	double t = GenInitialMatrix[2][0]; //��ǰʱ�������Ϊ��3�е�1��

	// ��ʼ����ǰʱ�̵�abc���ѹ����
	for (int k=0; k<3; k++) {
		Vsabc[k] = GenInitialMatrix[2][ptr+k];
		Vsabc_his[k] = GenInitialMatrix[1][ptr+k];
		Vsabc_his2[k] = GenInitialMatrix[0][ptr+k];
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Vrabc_his[k] = GenInitialMatrix[1][ptr+3+k];
		Vrabc_his2[k] = GenInitialMatrix[0][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// ��еת�ء�ת�ٺ�ת�ӽ�
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	curAngle_s = w*Wbase*t;
	//curAngle_s = 0; // curAngle_s�ĳ�ֵ��������ѡȡ
	double Theta = GenInitialMatrix[2][ptr+13];
	curAngle_r = curAngle_s - Theta;

	// ��ʼ����ѹ������dq����
	double dqMatrix[3][12];

	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<3;i++)
	{
		for (int j=ptr;j<ptr+3;j++)
			GenInitialMatrix[i][j] /= Vbase;

		for (int j=ptr+3;j<ptr+6;j++)
			GenInitialMatrix[i][j] /= (Vbase/turnratio);

		for (int j=ptr+6;j<ptr+9;j++)
			GenInitialMatrix[i][j] /= Ibase;

		for (int j=ptr+9;j<ptr+12;j++)
			GenInitialMatrix[i][j] /= (Ibase*turnratio);
	}

	// ����dq����
	for (int i=0; i<3; i++)
	{
		// ����park�任�Ƕ�
		double theta_s = w*Wbase*GenInitialMatrix[i][0];
		double theta_r = theta_s - GenInitialMatrix[i][ptr+13];

		// ����park�任����
		parkTransMatrix(Ps,theta_s);
		parkTransMatrix(Pr,theta_r);

		// park�任
		// ���ӵ�ѹ
		for (int j=0; j<3; j++)
		{
			dqMatrix[i][j] = Ps(j+1,1)*GenInitialMatrix[i][ptr] + Ps(j+1,2)*GenInitialMatrix[i][ptr+1]
			+ Ps(j+1,3)*GenInitialMatrix[i][ptr+2];
		}
		// ת�ӵ�ѹ
		for (int j=3; j<6; j++)
		{
			dqMatrix[i][j] = Pr(j-2,1)*GenInitialMatrix[i][ptr+3] + Pr(j-2,2)*GenInitialMatrix[i][ptr+4]
			+ Pr(j-2,3)*GenInitialMatrix[i][ptr+5];
		}
		// ���ӵ���
		for (int j=6; j<9; j++)
		{
			dqMatrix[i][j] = Ps(j-5,1)*GenInitialMatrix[i][ptr+6] + Ps(j-5,2)*GenInitialMatrix[i][ptr+7]
			+ Ps(j-5,3)*GenInitialMatrix[i][ptr+8];
		}
		// ת�ӵ�ѹ
		for (int j=9; j<12; j++)
		{
			dqMatrix[i][j] = Pr(j-8,1)*GenInitialMatrix[i][ptr+9] + Pr(j-8,2)*GenInitialMatrix[i][ptr+10]
			+ Pr(j-8,3)*GenInitialMatrix[i][ptr+11];
		}
	}

	// dq������ֵ
	for (int k=0; k<3; k++) {
		Vsdq0[k] = dqMatrix[2][k];
		Vsdq0_his[k] = dqMatrix[1][k];
		Vsdq0_his2[k] = dqMatrix[0][k];
		Isdq0[k] = dqMatrix[2][k+6];
		Isdq0_his[k] = dqMatrix[1][k+6];
		Isdq0_his2[k] = dqMatrix[0][k+6];
		Vrdq0[k] = dqMatrix[2][k+3];
		Vrdq0_his[k] = dqMatrix[1][k+3];
		Vrdq0_his2[k] = dqMatrix[0][k+3];
		Irdq0[k] = dqMatrix[2][k+9];
		Irdq0_his[k] = dqMatrix[1][k+9];
		Irdq0_his2[k] = dqMatrix[0][k+9];
	}

	// �ڵ㵼�ɾ������
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq0_his2[0]*Isdq0_his2[1]-Irdq0_his2[1]*Isdq0_his2[0]);
	Te=0.5*Xm*(Irdq0_his[0]*Isdq0_his[1]-Irdq0_his[1]*Isdq0_his[0]);

	ptr += 15;
}

// �ڷ���������ģʽ�£�������������ת��
void WoundInducGenerator2::calculateTwind()
{
#ifdef WIND_VELOCITY_DATA_INPUT
	setWindVelocity();
#endif
	TSR=2.237*vWind/(wr*Wbase/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}

//������ת��, By Gao Haixiang
void WoundInducGenerator2::saveMachineWr(double** machineWrMatrix, int& ptr, int counter)
{
	machineWrMatrix[counter-1][ptr] = wr;
	ptr++;
}

//	��ÿ��ʱ���ķ������ݱ�����ÿ̨����� ,By Gao Haixiang
void WoundInducGenerator2::getWindVelocityData(double** VwMatrix, int rows, int& ptr)
{
	WindVelocityVector = new double[rows];
	for (int i=0;i<rows;i++)
	{
		WindVelocityVector[i] = VwMatrix[i][ptr];
	}
	ptr++;
}

// �趨����, By Gao Haixiang
void WoundInducGenerator2::setWindVelocity()
{
	vWind = WindVelocityVector[WindVelocityCounter];
	WindVelocityCounter++;
}