#include "MWSynGenerator12.h"

#include <cmath>
#include <iostream>
using namespace std;

MWSynGenerator12::MWSynGenerator12(int firstNode,int lastNode)
{
	//�Զ���Ԫ��
	isUserDef = 1;
	need_NEC=1;
	/************************************************************************/
	/*              ��EMTP�ӿ�����ı���                                    */
	/************************************************************************/
	type=12;
	nPort=12;
	nodeNumber=new int[nPort];
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}

	//���캯��,ai��bi��ciΪ�ڵ���,i=1,2,3,4
	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	w = 2*PI*120;//����ٶȣ�����ֵ,f=120Hz
	Ef = 1.0;//���ŵ�ѹ������ֵ
	//��������й����޹�������ֵ
	P0 = 22000000e6;//W
	Q0 = 19034783e6;//Var

	P0 = P0/4 ;//W
	Q0 = Q0/4 ;//Var

	Sm = 0;//���ڹ��ʣ�����ֵ���ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Vm = 1.0;//���ѹ��ֵ������ֵ
	Im = 0;//�������ֵ������ֵ���ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	ph = 5;//A1���ѹ�ĳ�ʼ��λ�����ȱ�ʾ

	Sangle = 0;//�������ؽǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Vangle = 0;//A���ѹ��ʼ��ǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Iangle = 0;//A�������ʼ��ǣ��ɼ���õ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	Angle = 0;//ת�ӽǵĳ�ʼֵ���ɼ���õ���//added by chenlj ���øñ���������Ϊ��ʼ�н�

	//dq0����ϵ�µĵ������,����ֵ
	Ra = 0.00621;
	Rf = 0.00025765;
	RD = 0.003033;
	RQ = 0.0030469;
	//dq0����ϵ�µĵ翹����,����ֵ
	Xad = 0.34288;
	Xaq = 0.15453;
	Xd = 0.40075;
	Xdm1 = 0.342877;
	Xdm2 = 0.342877;
	Xq = 0.2124;
	Xqm1 = 0.15459;
	Xqm2 = 0.15459;
	Xf = 0.37918;
	XD = 0.35411;
	XQ = 0.162953;
	XfD = Xad;
	Xaf = Xad;
	XaD = Xad;
	XaQ = Xaq;

	//���Ӳ��ѹ���������迹�Ļ�ֵ
	Vbase = 466e3*sqrt(2.0);//V�޸�
	Ibase = 20820e3*sqrt(2.0);//A�޸�
	Zbase = Vbase/Ibase;

	//��е�����
	H = 3.0;//s
	D = 0.01;
	Tm = 1.0;//����ת��

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//����coff��T
	coff = 99.0/101.0;
	T = 1.0/deltaT/w*0.5*(1+coff);

	//���ŵ�ѹ��ϵ��
	Vf=Rf/Xad;

	//�ɿ˱任�����������
	P.ResizeTo(1,8,1,12);
	invP.ResizeTo(1,12,1,8);

	for (int k=0;k<8;k++)
	{
		for (int i=0;i<3;i++)
		{
			P_1[k][i] = 0;
		}
		
	}

	for (int k=0;k<12;k++)
	{
		for (int i=0;i<2;i++)
		{
			invP_1[k][i] = 0;
		}

	}

	//
	AA.ResizeTo(1,11,1,11);
	BB.ResizeTo(1,11,1,11);
	////
	Rss.ResizeTo(1,8,1,8);
	Rsr.ResizeTo(1,8,1,3);
	Rrs.ResizeTo(1,3,1,8);
	Rrr.ResizeTo(1,3,1,3);

	////
	Bss.ResizeTo(1,8,1,8);
	Bsr.ResizeTo(1,8,1,3);
	Brs.ResizeTo(1,3,1,8);
	Brr.ResizeTo(1,3,1,3);

	////����ϵ������ʹ��ת����abc����ϵʱ�ܽ���
	Rs.ResizeTo(1,8,1,8);
	R_ave.ResizeTo(1,8,1,8);
	R_res.ResizeTo(1,8,1,8);
	inv_Rs.ResizeTo(1,8,1,8);
	inv_Rr.ResizeTo(1,3,1,3);

	////
	Gs1.ResizeTo(1,8,1,8);
	Gs2.ResizeTo(1,8,1,3);
	Gs3.ResizeTo(1,8,1,8);
	Gs4.ResizeTo(1,8,1,3);
	Gs5.ResizeTo(1,8,1,3);
	Gs6.ResizeTo(1,8,1,8);
	Gs13.ResizeTo(1,8,1,8);
	Gs25.ResizeTo(1,8,1,3);

	////
	Gr1.ResizeTo(1,3,1,3);
	Gr2.ResizeTo(1,3,1,8);
	Gr3.ResizeTo(1,3,1,8);
	Gr4.ResizeTo(1,3,1,3);
	Gr5.ResizeTo(1,3,1,3);
	Gr15.ResizeTo(1,3,1,3);

	/************************************************************************/
	/*              ��ʱ����Ҫ�����������                                  */
	/************************************************************************/
	//���¸�����������Ϊ����ֵ�����������0����ʼ��ʱ������ֵ
	/*  ������������  */
	//abc����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=0;i<12;i++)
	{
		Vabc[i] = 0;
		Vabc_pu[i] = 0;
		Vabc_1[i] = 0;
		Vabc_2[i] = 0;

		Iabc[i] = 0;
		Iabc_pu[i] = 0;
		Iabc_1[i] = 0;
		Iabc_2[i] = 0;
	}
	//dq0����ϵ�¶�����Ȧ�ĵ�ѹ����
	for (int i=1;i<=8;i++)
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

	/************************************************************************/
	/*              ��Ч��·�еĲ���                                        */
	/************************************************************************/
	nortonEquivalentConductance = 0;
	nortonEquivalentConductance_his = 0;//��Ч���ɵ���ʷֵ�����ڸ��½ڵ㵼�ɾ���

	for(int i=1;i<=8;i++)
	{
		nortonEquivalentCurrent_dq0[i-1] = 0;//dq0�����µ�ŵ�ٵ�Ч����������ֵ
	}

	/************************************************************************/
	/*              ����ϵ�������ŵ�ȵ�ֵ����                              */
	/************************************************************************/
	calculateCoefficientMatrix();

	for (int k=0;k<8;k++)
	{
		for (int i=0;i<8;i++)
		{
			Gs_1[k][i] = Gs1(k+1,i+1);
			Gs_3[k][i] = Gs3(k+1,i+1);
			Gs_6[k][i] = Gs6(k+1,i+1);
		}
		for (int j=0;j<3;j++)
		{
			Gs_4[k][j] = Gs4(k+1,j+1);
			Gs_25[k][j] = Gs25(k+1,j+1);
		}
	}
	
	for (int k=0;k<3;k++)
	{
		for (int i=0;i<8;i++)
		{
			Gr_2[k][i] = Gr2(k+1,i+1);
			Gr_3[k][i] = Gr3(k+1,i+1);
		}
		for (int j=0;j<3;j++)
		{
			Gr_4[k][j] = Gr4(k+1,j+1);
			Gr_15[k][j] = Gr15(k+1,j+1);
		}
	}
}
MWSynGenerator12::~MWSynGenerator12()
{//��������
}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/
void MWSynGenerator12::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{//��ʼ��֧·��ѹ����
	//��ʼ��wr��curAngle
	ph = ph*PI/180.0;//A1���ѹ�ĳ�ʼ��λ�����ȱ�ʾ
	Sm = sqrt(P0*P0+Q0*Q0);
	Sangle = atan(Q0/P0);
	Sm = Sm/(1.5*Vbase*Ibase);
	Im = Sm/Vm;
	Vangle = ph - PI/2.0;
	Iangle = Vangle - Sangle;
	Angle = ph/180*PI;//chenlj

	for(int k=1;k<=4;k++)
	{
		Iabc[3*k-3] = Im*cos(Iangle-(k-1)*PI/12);
		Iabc[3*k-2] = Im*cos(Iangle-(2.0/3.0)*PI-(k-1)*PI/12);
		Iabc[3*k-1] = Im*cos(Iangle+(2.0/3.0)*PI-(k-1)*PI/12);

		Vabc[3*k-3] = Vm*cos(Vangle-(k-1)*PI/12);
		Vabc[3*k-2] = Vm*cos(Vangle-(2.0/3.0)*PI-(k-1)*PI/12);
		Vabc[3*k-1] = Vm*cos(Vangle+(2.0/3.0)*PI-(k-1)*PI/12);
	}

	curAngle = Angle;
	parkTransMatrix(P,curAngle);

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<2;j++)
		{
			Vdq0[2*i+j] = P(2*i+j+1,3*i+1)*Vabc[3*i]+P(2*i+j+1,3*i+2)*Vabc[3*i+1]+P(2*i+j+1,3*i+3)*Vabc[3*i+2];
			Vdq0_his[2*i+j] = Vdq0[2*i+j];
			Vdq0_his2[2*i+j] = Vdq0[2*i+j];
			Idq0[2*i+j] = P(2*i+j+1,3*i+1)*Iabc[3*i]+P(2*i+j+1,3*i+2)*Iabc[3*i+1]+P(2*i+j+1,3*i+3)*Iabc[3*i+2];
			Idq0_his[2*i+j] = Idq0[2*i+j];
			Idq0_his2[2*i+j] = Idq0[2*i+j];
			Idq0_forcast[2*i+j] = Idq0[2*i+j];
		}
	}


	IfDQ[0] = (Vdq0[1]+Ra*Idq0[1]+Xd*Idq0[0])/Xad;//�޸�
	IfDQ[1] = 0;
	IfDQ[2] = 0;

	VfDQ[0] = IfDQ[0]*Rf;
	VfDQ[1] = 0;
	VfDQ[2] = 0;

	//����ŵ�ٵ�ֵ����
	TVectorD temp1,temp2,temp3,temp4;
	TVectorD Idq0_forcast_V,Idq0_V,IfDQ_V,Vdq0_V,VfDQ_V;
	TVectorD Inew;
	temp1.ResizeTo(1,8);
	temp2.ResizeTo(1,8);
	temp3.ResizeTo(1,8);
	temp4.ResizeTo(1,8);
	Idq0_forcast_V.ResizeTo(1,8);
	Idq0_V.ResizeTo(1,8);
	IfDQ_V.ResizeTo(1,3);
	Vdq0_V.ResizeTo(1,8);
	VfDQ_V.ResizeTo(1,3);
	Inew.ResizeTo(1,8);
	for (int i=1;i<=8;i++)
	{
		Idq0_forcast_V(i) = -Idq0_forcast[i-1];
		Idq0_V(i) = -Idq0[i-1];
		Vdq0_V(i) = Vdq0[i-1];
	}
	for (int i=1;i<=3;i++)
	{
		IfDQ_V(i) = IfDQ[i-1];
		VfDQ_V(i) = VfDQ[i-1];
	}
	temp1 = Gs1*Idq0_forcast_V;
	temp2 = Gs3*Idq0_V;
	temp1 = temp1 + temp2;
	temp2 = Gs4*IfDQ_V;
	temp3 = Gs6*Vdq0_V;
	temp4 = Gs25*VfDQ_V;

	Inew = temp1+temp2+temp3+temp4;

	for (int i=0;i<8;i++)
	{
		nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	}

	//����park��任����
	curAngle = Angle + deltaT*w;
	invParkTransMatrix(invP,curAngle);

	for (int k=0;k<12;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	for (int i=0;i<4;i++)
	{
		for (int j=1;j<=3;j++)
		{
			nortonEquivalentCurrent[3*i+j-1] = invP(3*i+j,2*i+1)*nortonEquivalentCurrent_dq0[2*i]
			+ invP(3*i+j,2*i+2)*nortonEquivalentCurrent_dq0[2*i+1];
			nortonEquivalentCurrent[3*i+j-1] = nortonEquivalentCurrent[3*i+j-1]*Ibase;
		}
	}

	//��PSCAD���ݳ�ʼ�����ӵ�ѹ����
	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<12;i++)
	{
		Iabc[i] = initialCurrentArray[ptr+i];
	}
	ptr+=12;
}
void MWSynGenerator12::readNodeVoltage(TVectorD& nodeVoltageArray)
{//�ӽڵ��ѹ�����ж�֧·�ڵ��ѹ
	for (int i=0;i<12;i++)
	{
		Vabc_1[i] =Vabc[i];
	}

	for (int i=0;i<12;i++)
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
void MWSynGenerator12::calculateBranchVoltage()
{//����֧·��ѹ

}

void MWSynGenerator12::calculateBranchCurrent()
{//����֧·����

	for (int i=0;i<12;i++)
	{
		Iabc_1[i] =Iabc[i];
		Iabc[i] = -Vabc[i]*nortonEquivalentConductance-nortonEquivalentCurrent[i];
	}
}

void MWSynGenerator12::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	curAngle = Angle +time*w;
	calculateDq0Results();//����dq0����ϵ�µĽ����Ԥ�����ӵ���,�����������park�任

	//TVectorD temp1,temp2,temp3,temp4;
	//TVectorD Idq0_forcast_V,Idq0_V,IfDQ_V,Vdq0_V,VfDQ_V;
	//TVectorD Inew;
	//temp1.ResizeTo(1,8);
	//temp2.ResizeTo(1,8);
	//temp3.ResizeTo(1,8);
	//temp4.ResizeTo(1,8);
	//Idq0_forcast_V.ResizeTo(1,8);
	//Idq0_V.ResizeTo(1,8);
	//IfDQ_V.ResizeTo(1,3);
	//Vdq0_V.ResizeTo(1,8);
	//VfDQ_V.ResizeTo(1,3);
	//Inew.ResizeTo(1,8);
	//for (int i=1;i<=8;i++)
	//{
	//	Idq0_forcast_V(i) = -Idq0_forcast[i-1];
	//	Idq0_V(i) = -Idq0[i-1];
	//	Vdq0_V(i) = Vdq0[i-1];
	//}
	//for (int i=1;i<=3;i++)
	//{
	//	IfDQ_V(i) = IfDQ[i-1];
	//	VfDQ_V(i) = VfDQ[i-1];
	//}
	//temp1 = Gs1*Idq0_forcast_V;
	//temp2 = Gs3*Idq0_V;
	//temp1 = temp1 + temp2;
	//temp2 = Gs4*IfDQ_V;
	//temp3 = Gs6*Vdq0_V;
	//temp4 = Gs25*VfDQ_V;

	//Inew = temp1+temp2+temp3+temp4;
  //
  //
	//for (int i=0;i<8;i++)
	//{
	//	nortonEquivalentCurrent_dq0[i] = Inew(i+1);
	//}


	for (int k=0;k<8;k++)
	{
		nortonEquivalentCurrent_dq0[k] = 0;
		for (int i=0;i<8;i++)
		{
			nortonEquivalentCurrent_dq0[k]=nortonEquivalentCurrent_dq0[k]-Gs_1[k][i]*Idq0_forcast[i]-Gs_3[k][i]*Idq0[i]+Gs_6[k][i]*Vdq0[i];
		}
		nortonEquivalentCurrent_dq0[k]=nortonEquivalentCurrent_dq0[k]+Gs_4[k][0]*IfDQ[0]+Gs_4[k][1]*IfDQ[1]+Gs_4[k][2]*IfDQ[2]+Gs_25[k][0]*VfDQ[0];
	}


	//����park��任����
	curAngle = Angle +time*w +deltaT*w;
	//invParkTransMatrix(invP,curAngle);
	invParkTransMatrix(invP_1,curAngle);

	for (int k=0;k<12;k++)
	{
		nortonEquivalentCurrent_1[k] = nortonEquivalentCurrent[k];
	}

	//����abc�����µ���ʷ������
	//for (int i=0;i<4;i++)
	//{
	//	for (int j=1;j<=3;j++)
	//	{
	//		nortonEquivalentCurrent[3*i+j-1] = invP(3*i+j,2*i+1)*nortonEquivalentCurrent_dq0[2*i]
	//		+ invP(3*i+j,2*i+2)*nortonEquivalentCurrent_dq0[2*i+1];
	//		nortonEquivalentCurrent[3*i+j-1] = nortonEquivalentCurrent[3*i+j-1]*Ibase;
	//	}
	//}

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<3;j++)
		{
			nortonEquivalentCurrent[3*i+j] = invP_1[3*i+j][0]*nortonEquivalentCurrent_dq0[2*i]
			+ invP_1[3*i+j][1]*nortonEquivalentCurrent_dq0[2*i+1];
			nortonEquivalentCurrent[3*i+j] = nortonEquivalentCurrent[3*i+j]*Ibase;
		}
	}
}

void MWSynGenerator12::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<12;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}
}
void MWSynGenerator12::formConductanceMatrix(TMatrixD& conductanceMatrix)
{//�γɽڵ㵼����
	int N;
	double tmpd;
	for (int i=0;i<12;i++)
	{
		N = nodeNumber[i];
		tmpd = conductanceMatrix(N,N);
		conductanceMatrix(N,N) = tmpd+this->nortonEquivalentConductance;
	}
}
void MWSynGenerator12::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
{//����֧·����
	for (int k=0;k<12;k++)
	{
		branchCurrentMatrix(counter,ptr+k) = Iabc[k];
	}
	ptr+=12;
}

void MWSynGenerator12::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{//����֧·����
	for (int k=0;k<12;k++)
	{
		branchCurrentMatrix_1[counter-1][ptr+k] = Iabc[k];
	}
	ptr+=12;
}

void MWSynGenerator12::interpolate(double ratio)//������ֵ����֧·�ĵ�ѹ�������в�ֵ
{
	//��ѹ����������ֵ
	for (int k=0;k<12;k++)
	{
		Vabc[k] = (1-ratio)*Vabc_1[k] + ratio*Vabc[k];
		Iabc[k] = (1-ratio)*Iabc_1[k] + ratio*Iabc[k];
	}

	//ŵ�ٵ�ֵ������ֵ
	double temp1,temp2;
	for (int k=0;k<12;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void MWSynGenerator12::updateResult(int updateTypeNum)
{
	double V_temp[12],I_temp[12];//��ʱ��������������ʱʹ��

	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for (int i=0;i<12;i++)
		{
			Vabc_2[i] = Vabc_1[i];
			Iabc_2[i] = Iabc_1[i];
			nortonEquivalentCurrent_2[i] = nortonEquivalentCurrent_1[i];
		}
		break;
	case 2://��_2��������ֵ����_1������
		for (int i=0;i<12;i++)
		{
			Vabc_1[i] = Vabc_2[i];
			Iabc_1[i] = Iabc_2[i];
			nortonEquivalentCurrent[i] = nortonEquivalentCurrent_1[i];
			nortonEquivalentCurrent_1[i] = nortonEquivalentCurrent_2[i];
		}	
		break;
	case 3://ŵ�ٵ�ֵ�����洢����_1��������ֵ����_2������
		for (int i=0;i<12;i++)
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
//����ϵ������AA,BB,Rss...Brr,Rs...inv_Rr,Gs1...Gs6,Gr1...Gr5
//ͬʱ����nortonEquivalentConductance
void MWSynGenerator12::calculateCoefficientMatrix()
{
	//����AA���󣬼��������
	AA(1,1) = Ra+2.0*T*Xd;
	AA(1,2) = -Xq;
	AA(1,3) = 2.0*T*Xdm1;
	AA(1,4) = -Xqm1;
	AA(1,5) = 2.0*T*Xdm2;
	AA(1,6) = -Xqm2;
	AA(1,7) = 2.0*T*Xdm1;
	AA(1,8) = -Xqm1;
	AA(1,9) = 2.0*T*Xaf;
	AA(1,10) = 2.0*T*XaD;
	AA(1,11) = -XaQ;

	AA(2,1) = Xd;
	AA(2,2) = Ra+2.0*T*Xq;
	AA(2,3) = Xdm1;
	AA(2,4) = 2.0*T*Xqm1;
	AA(2,5) = Xdm2;
	AA(2,6) = 2.0*T*Xqm2;
	AA(2,7) = Xdm1;
	AA(2,8) = 2.0*T*Xqm1;
	AA(2,9) = Xaf;
	AA(2,10) = XaD;
	AA(2,11) = 2.0*T*XaQ;

	AA(3,1) = 2.0*T*Xdm1;
	AA(3,2) = -Xqm1;
	AA(3,3) = Ra+2.0*T*Xd;
	AA(3,4) = -Xq;
	AA(3,5) = 2.0*T*Xdm1;
	AA(3,6) = -Xqm1;
	AA(3,7) = 2.0*T*Xdm2;
	AA(3,8) = -Xqm2;
	AA(3,9) = 2.0*T*Xaf;
	AA(3,10) = 2.0*T*XaD;
	AA(3,11) = -XaQ;

	AA(4,1) = Xdm1;
	AA(4,2) = 2.0*T*Xqm1;
	AA(4,3) = Xd;
	AA(4,4) = Ra+2.0*T*Xq;
	AA(4,5) = Xdm1;
	AA(4,6) = 2.0*T*Xqm1;
	AA(4,7) = Xdm2;
	AA(4,8) = 2.0*T*Xqm2;
	AA(4,9) = Xaf;
	AA(4,10) = XaD;
	AA(4,11) = 2.0*T*XaQ;

	AA(5,1) = 2.0*T*Xdm2;
	AA(5,2) = -Xqm2;
	AA(5,3) = 2.0*T*Xdm1;
	AA(5,4) = -Xqm1;
	AA(5,5) = Ra+2.0*T*Xd;
	AA(5,6) = -Xq;
	AA(5,7) = 2.0*T*Xdm1;
	AA(5,8) = -Xqm1;
	AA(5,9) = 2.0*T*Xaf;
	AA(5,10) = 2.0*T*XaD;
	AA(5,11) = -XaQ;

	AA(6,1) = Xdm2;
	AA(6,2) = 2.0*T*Xqm2;
	AA(6,3) = Xdm1;
	AA(6,4) = 2.0*T*Xqm1;
	AA(6,5) = Xd;
	AA(6,6) = Ra+2.0*T*Xq;
	AA(6,7) = Xdm1;
	AA(6,8) = 2.0*T*Xqm1;
	AA(6,9) = Xaf;
	AA(6,10) = XaD;
	AA(6,11) = 2.0*T*XaQ;

	AA(7,1) = 2.0*T*Xdm1;
	AA(7,2) = -Xqm1;
	AA(7,3) = 2.0*T*Xdm2;
	AA(7,4) = -Xqm2;
	AA(7,5) = 2.0*T*Xdm1;
	AA(7,6) = -Xqm1;
	AA(7,7) = Ra+2.0*T*Xd;
	AA(7,8) = -Xq;
	AA(7,9) = 2.0*T*Xaf;
	AA(7,10) = 2.0*T*XaD;
	AA(7,11) = -XaQ;

	AA(8,1) = Xdm1;
	AA(8,2) = 2.0*T*Xqm1;
	AA(8,3) = Xdm2;
	AA(8,4) = 2.0*T*Xqm2;
	AA(8,5) = Xdm1;
	AA(8,6) = 2.0*T*Xqm1;
	AA(8,7) = Xd;
	AA(8,8) = Ra+2.0*T*Xq;
	AA(8,9) = Xaf;
	AA(8,10) = XaD;
	AA(8,11) = 2.0*T*XaQ;

	AA(9,1) = 2.0*T*Xaf;
	AA(9,2) = 0;
	AA(9,3) = 2.0*T*Xaf;
	AA(9,4) = 0;
	AA(9,5) = 2.0*T*Xaf;
	AA(9,6) = 0;
	AA(9,7) = 2.0*T*Xaf;
	AA(9,8) = 0;
	AA(9,9) = Rf+2.0*T*Xf;
	AA(9,10) = 2.0*T*XfD;
	AA(9,11) = 0;

	AA(10,1) = 2.0*T*XaD;
	AA(10,2) = 0;
	AA(10,3) = 2.0*T*XaD;
	AA(10,4) = 0;
	AA(10,5) = 2.0*T*XaD;
	AA(10,6) = 0;
	AA(10,7) = 2.0*T*XaD;
	AA(10,8) = 0;
	AA(10,9) = 2.0*T*XfD;
	AA(10,10) = RD+2.0*T*XD;
	AA(10,11) = 0;

	AA(11,1) = 0;
	AA(11,2) = 2.0*T*XaQ;
	AA(11,3) = 0;
	AA(11,4) = 2.0*T*XaQ;
	AA(11,5) = 0;
	AA(11,6) = 2.0*T*XaQ;
	AA(11,7) = 0;
	AA(11,8) = 2.0*T*XaQ;
	AA(11,9) = 0;
	AA(11,10) = 0;
	AA(11,11) = RQ+2.0*T*XQ;

	//����BB����
	BB(1,1) = coff*Ra-2.0*T*Xd;
	BB(1,2) = -coff*Xq;
	BB(1,3) = -2.0*T*Xdm1;
	BB(1,4) = -coff*Xqm1;
	BB(1,5) = -2.0*T*Xdm2;
	BB(1,6) = -coff*Xqm2;
	BB(1,7) = -2.0*T*Xdm1;
	BB(1,8) = -coff*Xqm1;
	BB(1,9) = -2.0*T*Xaf;
	BB(1,10) = -2.0*T*XaD;
	BB(1,11) = -coff*XaQ;

	BB(2,1) = coff*Xd;
	BB(2,2) = coff*Ra-2.0*T*Xq;
	BB(2,3) = coff*Xdm1;
	BB(2,4) = -2.0*T*Xqm1;
	BB(2,5) = coff*Xdm2;
	BB(2,6) = -2.0*T*Xqm2;
	BB(2,7) = coff*Xdm1;
	BB(2,8) = -2.0*T*Xqm1;
	BB(2,9) = coff*Xaf;
	BB(2,10) = coff*XaD;
	BB(2,11) = -2.0*T*XaQ;

	BB(3,1) = -2.0*T*Xdm1;
	BB(3,2) = -coff*Xqm1;
	BB(3,3) = coff*Ra-2.0*T*Xd;
	BB(3,4) = -coff*Xq;
	BB(3,5) = -2.0*T*Xdm1;
	BB(3,6) = -coff*Xqm1;
	BB(3,7) = -2.0*T*Xdm2;
	BB(3,8) = -coff*Xqm2;
	BB(3,9) = -2.0*T*Xaf;
	BB(3,10) = -2.0*T*XaD;
	BB(3,11) = -coff*wr*XaQ;

	BB(4,1) = coff*Xdm1;
	BB(4,2) = -2.0*T*Xqm1;
	BB(4,3) = coff*Xd;
	BB(4,4) = coff*Ra-2.0*T*Xq;
	BB(4,5) = coff*Xdm1;
	BB(4,6) = -2.0*T*Xqm1;
	BB(4,7) = coff*Xdm2;
	BB(4,8) = -2.0*T*Xqm2;
	BB(4,9) = coff*Xaf;
	BB(4,10) = coff*XaD;
	BB(4,11) = -2.0*T*XaQ;

	BB(5,1) = -2.0*T*Xdm2;
	BB(5,2) = -coff*Xqm2;
	BB(5,3) = -2.0*T*Xdm1;
	BB(5,4) = -coff*Xqm1;
	BB(5,5) = coff*Ra-2.0*T*Xd;
	BB(5,6) = -coff*Xq;
	BB(5,7) = -2.0*T*Xdm1;
	BB(5,8) = -coff*Xqm1;
	BB(5,9) = -2.0*T*Xaf;
	BB(5,10) = -2.0*T*XaD;
	BB(5,11) = -coff*XaQ;

	BB(6,1) = coff*Xdm2;
	BB(6,2) = -2.0*T*Xqm2;
	BB(6,3) = coff*Xdm1;
	BB(6,4) = -2.0*T*Xqm1;
	BB(6,5) = coff*Xd;
	BB(6,6) = coff*Ra-2.0*T*Xq;
	BB(6,7) = coff*Xdm1;
	BB(6,8) = -2.0*T*Xqm1;
	BB(6,9) = coff*Xaf;
	BB(6,10) = coff*XaD;
	BB(6,11) = -2.0*T*XaQ;

	BB(7,1) = -2.0*T*Xdm1;
	BB(7,2) = -coff*Xqm1;
	BB(7,3) = -2.0*T*Xdm2;
	BB(7,4) = -coff*Xqm2;
	BB(7,5) = -2.0*T*Xdm1;
	BB(7,6) = -coff*Xqm1;
	BB(7,7) = coff*Ra-2.0*T*Xd;
	BB(7,8) = -coff*Xq;
	BB(7,9) = -2.0*T*Xaf;
	BB(7,10) = -2.0*T*XaD;
	BB(7,11) = -coff*XaQ;

	BB(8,1) = coff*Xdm1;
	BB(8,2) = -2.0*T*Xqm1;
	BB(8,3) = coff*Xdm2;
	BB(8,4) = -2.0*T*Xqm2;
	BB(8,5) = coff*Xdm1;
	BB(8,6) = -2.0*T*Xqm1;
	BB(8,7) = coff*Xd;
	BB(8,8) = coff*Ra-2.0*T*Xq;
	BB(8,9) = coff*Xaf;
	BB(8,10) = coff*XaD;
	BB(8,11) = -2.0*T*XaQ;

	BB(9,1) = -2.0*T*Xaf;
	BB(9,2) = 0;
	BB(9,3) = -2.0*T*Xaf;
	BB(9,4) = 0;
	BB(9,5) = -2.0*T*Xaf;
	BB(9,6) = 0;
	BB(9,7) = -2.0*T*Xaf;
	BB(9,8) = 0;
	BB(9,9) = coff*Rf-2.0*T*Xf;
	BB(9,10) = -2.0*T*XfD;
	BB(9,11) = 0;

	BB(10,1) = -2.0*T*XaD;
	BB(10,2) = 0;
	BB(10,3) = -2.0*T*XaD;
	BB(10,4) = 0;
	BB(10,5) = -2.0*T*XaD;
	BB(10,6) = 0;
	BB(10,7) = -2.0*T*XaD;
	BB(10,8) = 0;
	BB(10,9) = -2.0*T*XfD;
	BB(10,10) = coff*RD-2.0*T*XD;
	BB(10,11) = 0;

	BB(11,1) = 0;
	BB(11,2) = -2.0*T*XaQ;
	BB(11,3) = 0;
	BB(11,4) = -2.0*T*XaQ;
	BB(11,5) = 0;
	BB(11,6) = -2.0*T*XaQ;
	BB(11,7) = 0;
	BB(11,8) = -2.0*T*XaQ;
	BB(11,9) = 0;
	BB(11,10) = 0;
	BB(11,11) = coff*RQ-2.0*T*XQ;

	//�������Rss,Rsr,Rrs,Rrr,Bss,Bsr,Brs,Brr
	for (int i=1;i<=8;i++)
	{
		for (int j=1;j<=8;j++)
		{
			Rss(i,j) = AA(i,j);
			Bss(i,j) = BB(i,j);
		}
	}

	for (int i=1;i<=8;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rsr(i,j) = AA(i,j+8);
			Bsr(i,j) = BB(i,j+8);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=8;j++)
		{
			Rrs(i,j) = AA(i+8,j);
			Brs(i,j) = BB(i+8,j);
		}
	}

	for (int i=1;i<=3;i++)
	{
		for (int j=1;j<=3;j++)
		{
			Rrr(i,j) = AA(i+8,j+8);
			Brr(i,j) = BB(i+8,j+8);
		}
	}

	//����Rs,R_ave,R_res,inv_Rs,inv_Rr
	inv_Rr = Rrr;
	inv_Rr.Invert();

	Rs = Rss-Rsr*inv_Rr*Rrs;

	R_ave(1,1) = (Rs(1,1)+Rs(2,2))/2.0;
	R_ave(2,2) = (Rs(1,1)+Rs(2,2))/2.0;
	R_ave(3,3) = (Rs(3,3)+Rs(4,4))/2.0;
	R_ave(4,4) = (Rs(3,3)+Rs(4,4))/2.0;
	R_ave(5,5) = (Rs(5,5)+Rs(6,6))/2.0;
	R_ave(6,6) = (Rs(5,5)+Rs(6,6))/2.0;
	R_ave(7,7) = (Rs(7,7)+Rs(8,8))/2.0;
	R_ave(8,8) = (Rs(7,7)+Rs(8,8))/2.0;

	R_res = Rs - R_ave;

	//����ŵ�ٵ�Ч�絼
	double resistance = (Rs(1,1)+Rs(2,2))/2.0;
	nortonEquivalentConductance = 1/resistance;
	nortonEquivalentConductance /= Zbase;

	inv_Rs = R_ave;
	inv_Rs.Invert();

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
void MWSynGenerator12::calculateDq0Results()
{
	//��abc�����µĶ��ӵ�ѹ����ת��Ϊ����ֵ
	for (int i=0;i<12;i++)
	{
		Vabc_pu[i] = Vabc[i]/Vbase;
		Iabc_pu[i] = Iabc[i]/Ibase;
	}

	//����park�任����
	//parkTransMatrix(P,curAngle);
	parkTransMatrix(P_1,curAngle);

	//���㶨�ӵ�ѹ
	//for (int i=0;i<4;i++)
	//{
	//	for (int j=1;j<=2;j++)
	//	{
	//		Vdq0_his2[2*i+j-1] = Vdq0_his[2*i+j-1];
	//		Vdq0_his[2*i+j-1] = Vdq0[2*i+j-1];
	//		Vdq0[2*i+j-1] = P(2*i+j,3*i+1)*Vabc_pu[3*i] + P(2*i+j,3*i+2)*Vabc_pu[3*i+1] + P(2*i+j,3*i+3)*Vabc_pu[3*i+2];
	//	}
	//}

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<2;j++)
		{
			Vdq0_his2[2*i+j] = Vdq0_his[2*i+j];
			Vdq0_his[2*i+j] = Vdq0[2*i+j];
			Vdq0[2*i+j] = P_1[2*i+j][0]*Vabc_pu[3*i] + P_1[2*i+j][1]*Vabc_pu[3*i+1] + P_1[2*i+j][2]*Vabc_pu[3*i+2];
		}
	}

	//���㶨�ӵ���
	//for (int i=0;i<4;i++)
	//{
	//	for (int j=1;j<=2;j++)
	//	{
	//		Idq0_his2[2*i+j-1] = Idq0_his[2*i+j-1];
	//		Idq0_his[2*i+j-1] = Idq0[2*i+j-1];
	//		Idq0[2*i+j-1] = P(2*i+j,3*i+1)*Iabc_pu[3*i] + P(2*i+j,3*i+2)*Iabc_pu[3*i+1] + P(2*i+j,3*i+3)*Iabc_pu[3*i+2];
	//	}
	//}
	for (int i=0;i<4;i++)
	{
		for (int j=0;j<2;j++)
		{
			Idq0_his2[2*i+j] = Idq0_his[2*i+j];
			Idq0_his[2*i+j] = Idq0[2*i+j];
			Idq0[2*i+j] = P_1[2*i+j][0]*Iabc_pu[3*i] + P_1[2*i+j][1]*Iabc_pu[3*i+1] + P_1[2*i+j][2]*Iabc_pu[3*i+2];
		}
	}




	//Ԥ����һʱ�̶��ӵ���
	for (int i=0;i<8;i++)
	{
		//Idq0_forcast[i] = (4.0/3.0)*Idq0[i]+(1.0/3.0)*Idq0_his[i]-(2.0/3.0)*Idq0_his2[i];
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

	//TVectorD Idq0_V,Idq0_his_V,IfDQ_his_V,VfDQ_V,temp5,temp6,temp7,temp8,Ifdq_V;
	//Idq0_V.ResizeTo(1,8);
	//Idq0_his_V.ResizeTo(1,8);
	//IfDQ_his_V.ResizeTo(1,3);
	//VfDQ_V.ResizeTo(1,3);
	//temp5.ResizeTo(1,3);
	//temp6.ResizeTo(1,3);
	//temp7.ResizeTo(1,3);
	//temp8.ResizeTo(1,3);
	//Ifdq_V.ResizeTo(1,3);
	//for (int i=1;i<=8;i++)
	//{
	//	Idq0_V(i) = -Idq0[i-1];
	//	Idq0_his_V(i) = -Idq0_his[i-1];
	//}
	//for (int i=1;i<=3;i++)
	//{
	//	IfDQ_his_V(i) = IfDQ_his[i-1];
	//	VfDQ_V(i) = VfDQ[i-1];
	//}
	//temp5 = Gr2*Idq0_V;
	//temp6 = Gr3*Idq0_his_V;
	//temp7 = Gr4*IfDQ_his_V;
	//temp8 = Gr15*VfDQ_V;
  //
	//Ifdq_V = temp5 + temp6 + temp7 + temp8;
  //
	//for (int i=0;i<3;i++)
	//{
	//	IfDQ[i] = Ifdq_V(i+1);
	//}
	
	for (int k=0;k<3;k++)
	{
		IfDQ[k] = 0;
		for (int i=0;i<8;i++)
		{
			IfDQ[k] = IfDQ[k]-Gr_2[k][i]*Idq0[i]-Gr_3[k][i]*Idq0_his[i];
		}
		IfDQ[k] = IfDQ[k]+Gr_4[k][0]*IfDQ_his[0]+Gr_4[k][1]*IfDQ_his[1]+Gr_4[k][2]*IfDQ_his[2]+Gr_15[k][0]*VfDQ[0];
	}
}


//�����ɿ˱任�뷴�任����
void MWSynGenerator12::parkTransMatrix(TMatrixD &T,double angle)
{
	for (int i=1;i<=8;i++)
	{
		for (int j=1;j<=12;j++)
		{
			T(i,j) = 0;
		}
	}

	for (int k=0;k<4;k++)
	{
		T(1+2*k,1+3*k) = 2.0/3.0*cos(angle-k*PI/12.0);
		T(1+2*k,2+3*k) = 2.0/3.0*cos(angle-2.0/3.0*PI-k*PI/12.0);
		T(1+2*k,3+3*k) = 2.0/3.0*cos(angle+2.0/3.0*PI-k*PI/12.0);
		T(2+2*k,1+3*k) = -2.0/3.0*sin(angle-k*PI/12.0);
		T(2+2*k,2+3*k) = -2.0/3.0*sin(angle-2.0/3.0*PI-k*PI/12.0);
		T(2+2*k,3+3*k) = -2.0/3.0*sin(angle+2.0/3.0*PI-k*PI/12.0);
	}

}
void MWSynGenerator12::invParkTransMatrix(TMatrixD &T_,double angle)
{
	for (int i=1;i<=12;i++)
	{
		for (int j=1;j<=8;j++)
		{
			T_(i,j) = 0;
		}
	}

	for (int k=0;k<4;k++)
	{
		T_(1+3*k,1+2*k) = cos(angle-k*PI/12.0);
		T_(1+3*k,2+2*k) = -sin(angle-k*PI/12.0);
		T_(2+3*k,1+2*k) = cos(angle-2.0/3.0*PI-k*PI/12.0);
		T_(2+3*k,2+2*k) = -sin(angle-2.0/3.0*PI-k*PI/12.0);
		T_(3+3*k,1+2*k) = cos(angle+2.0/3.0*PI-k*PI/12.0);
		T_(3+3*k,2+2*k) = -sin(angle+2.0/3.0*PI-k*PI/12.0);
	}
}

void MWSynGenerator12::parkTransMatrix(double P[8][3],double angle)
{
	for (int k=0;k<4;k++)
	{
		P[2*k][0] = 2.0/3.0*cos(angle-k*PI/12.0);
		P[2*k][1] = 2.0/3.0*cos(angle-2.0/3.0*PI-k*PI/12.0);
		P[2*k][2] = 2.0/3.0*cos(angle+2.0/3.0*PI-k*PI/12.0);
		P[1+2*k][0] = -2.0/3.0*sin(angle-k*PI/12.0);
		P[1+2*k][1] = -2.0/3.0*sin(angle-2.0/3.0*PI-k*PI/12.0);
		P[1+2*k][2] = -2.0/3.0*sin(angle+2.0/3.0*PI-k*PI/12.0);
	}
}
void MWSynGenerator12::invParkTransMatrix(double invP[12][2],double angle)
{
	for (int k=0;k<4;k++)
	{
		invP[3*k][0] = cos(angle-k*PI/12.0);
		invP[3*k][1] = -sin(angle-k*PI/12.0);
		invP[1+3*k][0] = cos(angle-2.0/3.0*PI-k*PI/12.0);
		invP[1+3*k][1] = -sin(angle-2.0/3.0*PI-k*PI/12.0);
		invP[2+3*k][0] = cos(angle+2.0/3.0*PI-k*PI/12.0);
		invP[2+3*k][1] = -sin(angle+2.0/3.0*PI-k*PI/12.0);
	}
}