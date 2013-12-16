#include "WoundInducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
WoundInducGenerator::WoundInducGenerator(int firstNode,int lastNode, int control,double vw)
{
	isUserDef = 0;
	need_NEC=1;
	//���캯����a��b��cΪ�ڵ���
	nPort=6;
	nodeNumber=new int[nPort];
	type=18;

	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	w = 2*PI*50;//����ٶȣ�����ֵ,f=60Hz
	w_pu=1;
	
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
	H = 3.0;
	//H = 0.25;
	D = 0;
	TL =-1.0;//����ת��

	//�������
	p_air=1.225;
	rArea=3177;
	GR=90.5;
	beta=0;
	GE=0.979;
	vWind=vw;

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//����coff��T
	coff = 99.0/101.0;
	//coff=1;
	//coff = 0.8;
	T = 1.0/deltaT/w*(1+coff);

	//�ɿ˱任�����������
	Ps.ResizeTo(1,2,1,3);
	invPs.ResizeTo(1,3,1,2);
	Pr.ResizeTo(1,2,1,3);
	invPr.ResizeTo(1,3,1,2);

	//
	AA.ResizeTo(1,4,1,4);
	BB.ResizeTo(1,4,1,4);
	////
	AA_ave.ResizeTo(1,4,1,4);
	AA_res.ResizeTo(1,4,1,4);
	inv_Aa.ResizeTo(1,4,1,4);
	////
	Gn1.ResizeTo(1,4,1,4);
	Gn2.ResizeTo(1,4,1,4);
	Gn3.ResizeTo(1,4,1,4);

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
	for (int i=1;i<=2;i++)
	{
		Vsdq0[i-1] = 0;
		Isdq0[i-1] = 0;
		Isdq0_his[i-1] = 0;
		Isdq0_his2[i-1] = 0;
		Isdq0_forcast[i-1] = 0;
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
	for (int i=1;i<=2;i++)
	{
		Vrdq0[i-1] = 0;
		Irdq0[i-1] = 0;
		Irdq0_his[i-1] = 0;
		Irdq0_his2[i-1] = 0;
		Irdq0_forcast[i-1] = 0;
	}

	/*  ��е��������  */
	if (control==1)
	{
		wr = 1.16715;//ת��ת��
		wr_his = 1.16715;//wr����ʷֵ
	} 
	else
	{
		wr = 0;//ת��ת��
		wr_his = 0;//wr����ʷֵ
	}
	curAngle_s= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_r= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ


	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	for (int i=0;i<2;i++)
	{
		nortonEquivalentConductance[i] = 0;
		nortonEquivalentConductance_his[i] = 0;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���
	}
	
	for(int i=1;i<=6;i++)
	{
		nortonEquivalentCurrent[i-1] = 0;//abc�����µ�ŵ�ٵ�Ч����������ֵ
		nortonEquivalentCurrent_1[i-1] = 0;
		nortonEquivalentCurrent_2[i-1] = 0;
	}
	for(int i=1;i<=4;i++)
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

	/************************************************************************/
	/*              ����ϵ�������ŵ�ȵ�ֵ����                              */
	/************************************************************************/
	//calculateCoefficientMatrix();
}
WoundInducGenerator::~WoundInducGenerator()
{//��������

}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/

void WoundInducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr)
{//��ʼ��֧·��ѹ����
	//��ʼ��wr��curAngle
    //��ʼ��**************************//
	//?????
	//**********************************
	
	T = 1.0/deltaT/w*(1+coff);
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
	}
	for (int i=0;i<3;i++)
	{
		Vrabc_1[i] =Vrabc[i];
	}

	for (int i=0;i<3;i++)
	{
		Isabc[i] = initialCurrentArray[ptr+i];
		Isabc_1[i] = Isabc[i];
	}
	ptr+=3;

	for (int i=0;i<3;i++)
	{
		Irabc[i] = initialCurrentArray[ptr+i];
		Irabc_1[i] = Irabc[i];
	}
	ptr+=3;

	calculateCoefficientMatrix();
	//���е���̣���ʼ��wr��curAngle

	parkTransMatrix(Ps,curAngle_s);
	parkTransMatrix(Pr,curAngle_r);

	for (int i=0;i<2;i++)
	{
		Vsdq0[i] = (Ps(i+1,1)*Vsabc[0]+Ps(i+1,2)*Vsabc[1]+Ps(i+1,3)*Vsabc[2])/Vbase;
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0_his2[i] = Vsdq0[i];
		Isdq0[i] = (Ps(i+1,1)*Isabc[0]+Ps(i+1,2)*Isabc[1]+Ps(i+1,3)*Isabc[2])/Ibase;
		Isdq0_his[i] = Isdq0[i];
		Isdq0_his2[i] = Isdq0[i];
		Isdq0_forcast[i] = Isdq0[i];
	}

	for (int i=0;i<2;i++)
	{
		Vrdq0[i] = (Pr(i+1,1)*Vrabc[0]+Pr(i+1,2)*Vrabc[1]+Pr(i+1,3)*Vrabc[2])/Vbase*turnratio;
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0_his2[i] = Vrdq0[i];
		Irdq0[i] = (Pr(i+1,1)*Irabc[0]+Pr(i+1,2)*Irabc[1]+Pr(i+1,3)*Irabc[2])/Ibase/turnratio;
		Irdq0_his[i] = Irdq0[i];
		Irdq0_his2[i] = Irdq0[i];
		Irdq0_forcast[i] = Irdq0[i];
	}

	//ŵ�ٵ�ֵ������ʼ��
	//����dq0�����µ���ʷ������
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 
			Gn1(i,1)*Isdq0_forcast[0]+Gn1(i,2)*Isdq0_forcast[1]
		+ Gn1(i,3)*Irdq0_forcast[0]+Gn1(i,4)*Irdq0_forcast[1]
		+Gn2(i,1)*Isdq0[0]+Gn2(i,2)*Isdq0[1]
		+Gn2(i,3)*Irdq0[0]+Gn2(i,4)*Irdq0[1]
		+ Gn3(i,1)*Vsdq0[0]+Gn3(i,2)*Vsdq0[1]
		+ Gn3(i,3)*Vrdq0[0]+Gn3(i,4)*Vrdq0[1];
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
	}

	//����park��任����
	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	//����abc�����µ���ʷ������
	for (int i=1;i<=6;i++)
	{
		if (i<=3)
		{
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[2]
			+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[3];
		    nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}
	}

	for (int k=0;k<6;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}
}

void WoundInducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
{//�ӽڵ��ѹ�����ж�֧·�ڵ��ѹ

	for (int i=0;i<3;i++)
	{
		Vsabc_1[i] =Vsabc[i];
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
void WoundInducGenerator::calculateBranchVoltage()
{//����֧·��ѹ
}

#ifdef PREDICT
void WoundInducGenerator::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance[0]+nortonEquivalentCurrent[i];
		Irabc_1[i] = Irabc[i];
		Irabc[i] = Vrabc[i]*nortonEquivalentConductance[1]+nortonEquivalentCurrent[i+3];
	}
}
#else
void WoundInducGenerator::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = nortonEquivalentCurrent[i];
		for (int j=0;j<3;j++) {
			Isabc[i] += Vsabc[j]*Yne_abc(i+1,j+1);
			Isabc[i] += Vrabc[j]*Yne_abc(i+1,j+4);
		}

		Irabc_1[i] = Irabc[i];
		Irabc[i] = nortonEquivalentCurrent[i+3];
		for (int j=0;j<3;j++) {
			Irabc[i] += Vsabc[j]*Yne_abc(i+4,j+1);
			Irabc[i] += Vrabc[j]*Yne_abc(i+4,j+4);
		}
	}
}
#endif

void WoundInducGenerator::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	//if (control==2)
	//{
	//	//������Ծ�仯�ķ���
	//	if (time<1.4-deltaT/2)
	//	{
	//		vWind=14.8;
	//	}
	//	else if (time<1.6-deltaT/2)
	//	{
	//		vWind=13;
	//	}
	//	else if (time<1.8-deltaT/2)
	//	{
	//		vWind=11;
	//	}
	//	else if (time<2.2-deltaT/2)
	//	{
	//		vWind=16.8;
	//	}
	//	else
	//	{
	//		vWind=14.8;
	//	}
	//}

	calculateCoefficientMatrix();
	calculateDq0Results();//����dq0����ϵ�µĽ����Ԥ�����ӵ���
	calculateOmega();

	//����dq0�����µ���ʷ������
#ifdef PREDICT
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
		nortonEquivalentCurrent_dq0[i-1] = 
			Gn1(i,1)*Isdq0_forcast[0]+Gn1(i,2)*Isdq0_forcast[1]
		+ Gn1(i,3)*Irdq0_forcast[0]+Gn1(i,4)*Irdq0_forcast[1]
		+ Gn2(i,1)*Isdq0[0]+Gn2(i,2)*Isdq0[1]
		+ Gn2(i,3)*Irdq0[0]+Gn2(i,4)*Irdq0[1]
		+ Gn3(i,1)*Vsdq0[0]+Gn3(i,2)*Vsdq0[1]
		+ Gn3(i,3)*Vrdq0[0]+Gn3(i,4)*Vrdq0[1];
	}
#else
	for (int i=1;i<=4;i++)
	{
		nortonEquivalentCurrent_dq0_his[i-1] = nortonEquivalentCurrent_dq0[i-1];
		nortonEquivalentCurrent_dq0[i-1] =
			  Gn2(i,1)*Isdq0[0] + Gn2(i,2)*Isdq0[1]
			+ Gn2(i,3)*Irdq0[0] + Gn2(i,4)*Irdq0[1]
			+ Gn3(i,1)*Vsdq0[0] + Gn3(i,2)*Vsdq0[1]
			+ Gn3(i,3)*Vrdq0[0] + Gn3(i,4)*Vrdq0[1];
	}
#endif

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
			nortonEquivalentCurrent[i-1] = invPs(i,1)*nortonEquivalentCurrent_dq0[0]
			+ invPs(i,2)*nortonEquivalentCurrent_dq0[1];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase;
		}
		else
		{
			nortonEquivalentCurrent[i-1] = invPr(i-3,1)*nortonEquivalentCurrent_dq0[2]
			+ invPr(i-3,2)*nortonEquivalentCurrent_dq0[3];
			nortonEquivalentCurrent[i-1] = nortonEquivalentCurrent[i-1]*Ibase*turnratio;
		}		
	}
}

void WoundInducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}

	// nodeNortonEquivalentCurrentArray.Print();
}
#ifdef PREDICT
void WoundInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N;
	double tmpd;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		if (i<3)
		{
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[0];
		} 
		else
		{
			conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance[1];
		}
	}
}
#else
void WoundInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N1, N2;
	for (int i=0;i<6;i++)
	{
		N1 = nodeNumber[i];
		for (int j=0; j<6; j++)
		{
			N2 = nodeNumber[j];
			conductanceMatrix(N1,N2) += Yne_abc(i+1,j+1);
		}
	}
}
#endif

void WoundInducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
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

void WoundInducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
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
void WoundInducGenerator::interpolate(double ratio)
{//������ֵ����֧·�ĵ�ѹ�������в�ֵ

	//��ѹ����������ֵ
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
	}

	for (int k=0;k<2;k++)
	{
		Vsdq0[k] = (1-ratio)*Vsdq0_his[k] + ratio*Vsdq0[k];
		Isdq0[k] = (1-ratio)*Isdq0_his[k] + ratio*Isdq0[k];
		Vrdq0[k] = (1-ratio)*Vrdq0_his[k] + ratio*Vrdq0[k];
		Irdq0[k] = (1-ratio)*Irdq0_his[k] + ratio*Irdq0[k];
	}

	for (int k=0;k<4;k++)
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

void WoundInducGenerator::updateResult(int updateTypeNum)
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
		}
		for (int i=0;i<6;i++)
		{
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
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

void WoundInducGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
{//ת��ת�ٱ仯ʱ���µ��ɾ���
	int N;
	double tmpd;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		if (i<3)
		{
			conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his[0]+this->nortonEquivalentConductance[0];
		} 
		else
		{	
			conductanceMatrix(N,N) = tmpd-this->nortonEquivalentConductance_his[1]+this->nortonEquivalentConductance[1];
		}
	}
}

/************************************************************************/
/*                ������ڲ����õĺ���                                  */
/************************************************************************/
//����ϵ������AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs4,Gr1...Gr3
//ͬʱ����nortonEquivalentConductance
#ifdef PREDICT
void WoundInducGenerator::calculateCoefficientMatrix()
{
	//// �ٶ�Ԥ��
	//double wr = 2*this->wr - this->wr_his;

	//����AA���󣬼��������
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w_pu*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w_pu*Xm;

	AA(2,1) = w_pu*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w_pu*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w_pu)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w_pu)*(Xm+Xlr);

	AA(4,1) = (w_pu-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w_pu-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//����BB����
	BB(1,1) = coff*Rs-T*(Xm+Xls);
	BB(1,2) = -coff*w_pu*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -coff*w_pu*Xm;

	BB(2,1) = coff*w_pu*(Xm+Xls);
	BB(2,2) = coff*Rs-T*(Xm+Xls);
	BB(2,3) = coff*w_pu*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = coff*(wr-w_pu)*Xm;
	BB(3,3) = coff*Rr-T*(Xm+Xlr);
	BB(3,4) = coff*(wr-w_pu)*(Xm+Xlr);

	BB(4,1) = coff*(w_pu-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = coff*(w_pu-wr)*(Xm+Xlr);
	BB(4,4) = coff*Rr-T*(Xm+Xlr);

	double Rsn=Rs+T*(Xm+Xls);
	double Rrn=Rr+T*(Xm+Xlr);
	AA_ave(1,1) = Rsn;
	AA_ave(2,2) = Rsn;
	AA_ave(3,3) = Rrn;
	AA_ave(4,4) = Rrn;
	for (int i=1;i<=4;i++)
	{
		for(int j=1;j<=4;j++)
		{
			if (i!=j)
			{
				AA_ave(i,j)=0;
			}
		}
	}

	AA_res = AA - AA_ave;

	//����ŵ�ٵ�Ч�絼
	nortonEquivalentConductance[0] = 1/Rsn;
	nortonEquivalentConductance[1] = 1/Rrn;
	inv_Aa = AA_ave;
	inv_Aa.Invert();
	//����Ч�ӱ���ֵ���㵽����ֵ
	nortonEquivalentConductance[0] /= Zbase;
	nortonEquivalentConductance[1] = nortonEquivalentConductance[1]/Zbase*turnratio*turnratio;

	//����Gn1...Gn3
	Gn1 = (-1.0)*inv_Aa*AA_res;
	Gn2 = (-1.0)*inv_Aa*BB;
	Gn3 = coff*inv_Aa;
}
#else
void WoundInducGenerator::calculateCoefficientMatrix()
{
	//����AA���󣬼��������
	AA(1,1) = Rs+T*(Xm+Xls);
	AA(1,2) = -w_pu*(Xm+Xls);
	AA(1,3) = T*Xm;
	AA(1,4) = -w_pu*Xm;

	AA(2,1) = w_pu*(Xm+Xls);
	AA(2,2) = Rs+T*(Xm+Xls);
	AA(2,3) = w_pu*Xm;
	AA(2,4) = T*Xm;

	AA(3,1) = T*Xm;
	AA(3,2) = (wr-w_pu)*Xm;
	AA(3,3) = Rr+T*(Xm+Xlr);
	AA(3,4) = (wr-w_pu)*(Xm+Xlr);

	AA(4,1) = (w_pu-wr)*Xm;
	AA(4,2) = T*Xm;
	AA(4,3) = (w_pu-wr)*(Xm+Xlr);
	AA(4,4) = Rr+T*(Xm+Xlr);

	//����BB����
	BB(1,1) = coff*Rs-T*(Xm+Xls);
	BB(1,2) = -coff*w_pu*(Xm+Xls);
	BB(1,3) = -T*Xm;
	BB(1,4) = -coff*w_pu*Xm;

	BB(2,1) = coff*w_pu*(Xm+Xls);
	BB(2,2) = coff*Rs-T*(Xm+Xls);
	BB(2,3) = coff*w_pu*Xm;
	BB(2,4) = -T*Xm;

	BB(3,1) = -T*Xm;
	BB(3,2) = coff*(wr-w_pu)*Xm;
	BB(3,3) = coff*Rr-T*(Xm+Xlr);
	BB(3,4) = coff*(wr-w_pu)*(Xm+Xlr);

	BB(4,1) = coff*(w_pu-wr)*Xm;
	BB(4,2) = -T*Xm;
	BB(4,3) = coff*(w_pu-wr)*(Xm+Xlr);
	BB(4,4) = coff*Rr-T*(Xm+Xlr);

	Yne_dq.ResizeTo(1,4,1,4);
	Yne_abc.ResizeTo(1,6,1,6);

	Yne_dq = AA;
	Yne_dq.Invert();

	//����Gn1...Gn3
	Gn2 = (-1.0)*Yne_dq*BB;
	Gn3 = coff*Yne_dq;
 }
#endif

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
#ifdef PREDICT
void WoundInducGenerator::calculateDq0Results()
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
	for (int i=0;i<2;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//�������
	for (int i=0;i<2;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}

	//Ԥ����һʱ�̵���
	for (int i=0;i<2;i++)
	{
		 Isdq0_forcast[i] = 1.3333333*Isdq0[i]+0.3333333*Isdq0_his[i]-0.6666666*Isdq0_his2[i];
		 Irdq0_forcast[i] = 1.3333333*Irdq0[i]+0.3333333*Irdq0_his[i]-0.6666666*Irdq0_his2[i];
		 //Isdq0_forcast[i] = 2*Isdq0[i] - Isdq0_his[i];
		 //Irdq0_forcast[i] = 2*Irdq0[i] - Irdq0_his[i];
		 //Isdq0_forcast[i] = 1.25*Isdq0[i]+0.5*Isdq0_his[i]-0.75*Isdq0_his2[i];
		 //Irdq0_forcast[i] = 1.25*Irdq0[i]+0.5*Irdq0_his[i]-0.75*Irdq0_his2[i];
	}
}
#else
void WoundInducGenerator::calculateDq0Results()
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
	for (int i=0;i<2;i++)
	{
		Vsdq0_his2[i] = Vsdq0_his[i];
		Vsdq0_his[i] = Vsdq0[i];
		Vsdq0[i] = Ps(i+1,1)*Vsabc_pu[0]+Ps(i+1,2)*Vsabc_pu[1]+Ps(i+1,3)*Vsabc_pu[2];
		Vrdq0_his2[i] = Vrdq0_his[i];
		Vrdq0_his[i] = Vrdq0[i];
		Vrdq0[i] = Pr(i+1,1)*Vrabc_pu[0]+Pr(i+1,2)*Vrabc_pu[1]+Pr(i+1,3)*Vrabc_pu[2];
	}

	//�������
	for (int i=0;i<2;i++)
	{
		Isdq0_his2[i] = Isdq0_his[i];
		Isdq0_his[i] = Isdq0[i];
		Isdq0[i] = Ps(i+1,1)*Isabc_pu[0]+Ps(i+1,2)*Isabc_pu[1]+Ps(i+1,3)*Isabc_pu[2];
		Irdq0_his2[i] = Irdq0_his[i];
		Irdq0_his[i] = Irdq0[i];
		Irdq0[i] = Pr(i+1,1)*Irabc_pu[0]+Pr(i+1,2)*Irabc_pu[1]+Pr(i+1,3)*Irabc_pu[2];
	}

	//����abc�����µĽڵ㵼�ɾ���
	TMatrixD pt(1,4,1,6);
	TMatrixD ipt(1,6,1,4);

	for (int i=1;i<=2;i++) {
		for (int j=1;j<=3;j++) {
			pt(i,j) = Ps(i,j);
			pt(i+2,j+3) = Pr(i,j);
		}
	}

	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<=3;i++) {
		for (int j=1;j<=2;j++) {
			ipt(i,j) = invPs(i,j);
			ipt(i+3,j+2) = invPr(i,j);
		}
	}

	Yne_abc = ipt*Yne_dq*pt;

	for (int i=1; i<=6; i++)
		for (int j=1; j<=6; j++)
			Yne_abc(i,j) /= Zbase;
}
#endif

//����е���̣�����ת��ת��
//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
void WoundInducGenerator::calculateOmega()
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
		curAngle_r = curAngle_r + deltaT*w*(w_pu-(wr_his+wr)/2);
		curAngle_s = curAngle_s + deltaT*w_pu*w;
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(w_pu-wr)*w;
		curAngle_s = curAngle_s + deltaT*w_pu*w;
	}
}

//�����ɿ˱任�뷴�任����
void WoundInducGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void WoundInducGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}


// ƽ��ģ����أ�xuyin��20121224
// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
void WoundInducGenerator::storeInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc_bak[k] = Vsabc[k];
		Isabc_bak[k] = Isabc[k];
		Vrabc_bak[k] = Vrabc[k];
		Irabc_bak[k] = Irabc[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq0_bak[k] = Vsdq0[k];
		Isdq0_bak[k] = Isdq0[k];
		Isdq0_his_bak[k] = Isdq0_his[k];
		Isdq0_his2_bak[k] = Isdq0_his2[k];
		Vrdq0_bak[k] = Vrdq0[k];
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
}

// ����У��ʱ�ָ��ڲ�����
void WoundInducGenerator::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Vsabc[k] = Vsabc_bak[k];
		Isabc[k] = Isabc_bak[k];
		Vrabc[k] = Vrabc_bak[k];
		Irabc[k] = Irabc_bak[k];
	}
	
	for (int k=0; k<2; k++) {
		Vsdq0[k] = Vsdq0_bak[k];
		Isdq0[k] = Isdq0_bak[k];
		Isdq0_his[k] = Isdq0_his_bak[k];
		Isdq0_his2[k] = Isdq0_his2_bak[k];
		Vrdq0[k] = Vrdq0_bak[k];
		Irdq0[k] = Irdq0_bak[k];
		Irdq0_his[k] = Irdq0_his_bak[k];
		Irdq0_his2[k] = Irdq0_his2_bak[k];
	}

	wr = wr_bak;
	wr_his = wr_his_bak;
	Te=Te_bak;
	Te_his=Te_his_bak;
	curAngle_s = curAngle_s_bak;
	curAngle_r = curAngle_r_bak;
}

// �����ʼ��ר�ú���, xuyin, 20121226
void WoundInducGenerator::initializeGen(double** GenInitialMatrix, int& ptr)
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
		Vrabc[k] = GenInitialMatrix[2][ptr+3+k];
		Isabc[k] = GenInitialMatrix[2][ptr+6+k];
		Irabc[k] = GenInitialMatrix[2][ptr+9+k];
	}

	// ��еת�ء�ת�ٺ�ת�ӽ�
	TL = GenInitialMatrix[2][ptr+12];
	wr = GenInitialMatrix[2][ptr+14];
	wr_his = GenInitialMatrix[1][ptr+14];
	curAngle_s = w*w_pu*t;
	//curAngle_s = 0; // curAngle_s�ĳ�ֵ��������ѡȡ
	double Theta = GenInitialMatrix[2][ptr+13];
	curAngle_r = curAngle_s - Theta;

	// ��ʼ����ѹ������dq����
	double dqMatrix[3][8];

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
		double theta_s = w*GenInitialMatrix[i][0];
		double theta_r = theta_s - GenInitialMatrix[i][ptr+13];

		// ����park�任����
		parkTransMatrix(Ps,theta_s);
		parkTransMatrix(Pr,theta_r);

		// park�任
		// ���ӵ�ѹ
		for (int j=0; j<2; j++)
		{
			dqMatrix[i][j] = Ps(j+1,1)*GenInitialMatrix[i][ptr] + Ps(j+1,2)*GenInitialMatrix[i][ptr+1]
				+ Ps(j+1,3)*GenInitialMatrix[i][ptr+2];
		}
		// ת�ӵ�ѹ
		for (int j=2; j<4; j++)
		{
			dqMatrix[i][j] = Pr(j-1,1)*GenInitialMatrix[i][ptr+3] + Pr(j-1,2)*GenInitialMatrix[i][ptr+4]
				+ Pr(j-1,3)*GenInitialMatrix[i][ptr+5];
		}
		// ���ӵ���
		for (int j=4; j<6; j++)
		{
			dqMatrix[i][j] = Ps(j-3,1)*GenInitialMatrix[i][ptr+6] + Ps(j-3,2)*GenInitialMatrix[i][ptr+7]
				+ Ps(j-3,3)*GenInitialMatrix[i][ptr+8];
		}
		// ת�ӵ�ѹ
		for (int j=6; j<8; j++)
		{
			dqMatrix[i][j] = Pr(j-5,1)*GenInitialMatrix[i][ptr+9] + Pr(j-5,2)*GenInitialMatrix[i][ptr+10]
				+ Pr(j-5,3)*GenInitialMatrix[i][ptr+11];
		}
	}
	
	// dq������ֵ
	for (int k=0; k<2; k++) {
		Vsdq0[k] = dqMatrix[2][k];
		Isdq0[k] = dqMatrix[2][k+4];
		Isdq0_his[k] = dqMatrix[1][k+4];
		Isdq0_his2[k] = dqMatrix[0][k+4];
		Vrdq0[k] = dqMatrix[2][k+2];
		Irdq0[k] = dqMatrix[2][k+6];
		Irdq0_his[k] = dqMatrix[1][k+6];
		Irdq0_his2[k] = dqMatrix[0][k+6];
	}

	// �ڵ㵼�ɾ������
	calculateCoefficientMatrix();
	Te_his=0.5*Xm*(Irdq0_his2[0]*Isdq0_his2[1]-Irdq0_his2[1]*Isdq0_his2[0]);
	Te=0.5*Xm*(Irdq0_his[0]*Isdq0_his[1]-Irdq0_his[1]*Isdq0_his[0]);

#ifndef PREDICT
	TMatrixD pt(1,4,1,6);
	TMatrixD ipt(1,6,1,4);

	for (int i=1;i<=2;i++) {
		for (int j=1;j<=3;j++) {
			pt(i,j) = Ps(i,j);
			pt(i+2,j+3) = Pr(i,j);
		}
	}

	invParkTransMatrix(invPs,curAngle_s);
	invParkTransMatrix(invPr,curAngle_r);

	for (int i=1;i<=3;i++) {
		for (int j=1;j<=2;j++) {
			ipt(i,j) = invPs(i,j);
			ipt(i+3,j+2) = invPr(i,j);
		}
	}

	Yne_abc = ipt*Yne_dq*pt;

	for (int i=1; i<=6; i++)
		for (int j=1; j<=6; j++)
			Yne_abc(i,j) /= Zbase;
#endif

	ptr += 15;
}

// �ڷ���������ģʽ�£�������������ת��
void WoundInducGenerator::calculateTwind()
{
	TSR=2.237*vWind/(wr*w/GR);
	Cp=0.5*(TSR-0.022*beta*beta-5.6)*exp(-0.17*TSR);
	Pwind=0.5*Cp*rArea*(vWind*vWind*vWind)*p_air*GE/Sbase;
	Twind=Pwind/wr;
}