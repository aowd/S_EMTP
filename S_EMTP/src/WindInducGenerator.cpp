#include "WindInducGenerator.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
WindInducGenerator::WindInducGenerator(int firstNode,int lastNode)
{
	isUserDef = 0;
	need_NEC=1;
	//���캯����a��b��cΪ�ڵ���
	nPort=6;
	nodeNumber=new int[nPort];
	type=19;

	/************************************************************************/
	/*                ��Ҫ���ⲿ����Ĳ���                                  */
	/************************************************************************/
	//���ڳ�ʼ���Ĳ���
	w = 2*PI*60;//����ٶȣ�����ֵ,f=60Hz
	w_pu=1;
	
	//dq0����ϵ�µĵ������,����ֵ
	Rs = 0.00621;
	Rr = 0.00025765;
	//dq0����ϵ�µĵ翹����,����ֵ
	Xm = 0.34288;
	Xls = 0.15453;
	Xlr = 0.40075;
	
	//���Ӳ��ѹ���������迹�Ļ�ֵ
	Vbase=10e3;
	Sbase=5e6;
	Vbase = (Vbase/sqrt(3.0))*sqrt(2.0);//V
	Ibase = Sbase/1.5/Vbase;//A
	Zbase = Vbase/Ibase;
	
	//��ת������������
	turnratio=1;

	//�������ת�ؿ���ģʽ
	control=0;

	//��е�����
	H_g = 0.5;
	D_g= 0.0;
	H_wt=2.5;//�ۺϺ����ֵ
	D_wt=0.0;//�ۺϺ����ֵ
	k=0.3;

	//�������
	p_air=1.229;
	R=40;
	//gr=1;
	gr=2.0;
	v=4;
	pitch_angle=2;//��λ����
	
	//H_wt=H_wt/gr/gr;
	//D_wt=D_wt/gr/gr;

	Co.ResizeTo(1,4,1,4);
	Co_inv.ResizeTo(1,4,1,4);
	Co(1,1)=2*H_g/deltaT+D_g/2;
	Co(1,2)=0;
	Co(1,3)=-k/2;
	Co(1,4)=-k/2;
	Co(2,1)=0;
	Co(2,2)=2*H_wt/deltaT+D_wt/2;
	Co(2,3)=k/2;
	Co(2,4)=k/2;
	Co(3,1)=w*deltaT/2;
	Co(3,2)=0;
	Co(3,3)=1;
	Co(3,4)=0;
	Co(4,1)=0;
	Co(4,2)=-w*deltaT/2;
	Co(4,3)=0;
	Co(4,4)=1;
	Co_inv=Co;
	Co_inv.Invert();
	bo.ResizeTo(1,4,1,1);
	co.ResizeTo(1,4,1,1);

	/************************************************************************/
	/*				ϵ����ϵ������,����������������ĵ�                   */
	/************************************************************************/
	//����coff��T
	//coff = 99/101;
	coff =0;
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
		wr = 1.05;//ת��ת��
		wr_his = 1.05;//wr����ʷֵ
	} 
	else
	{
		wr = 0;//ת��ת��
		wr_his = 0;//wr����ʷֵ
	}
	w_wt=0;
	//w_wt=1;
	//w_wt=w_wt*gr;
	w_wt_his=0;
	//w_wt_his=1;
	//w_wt_his=w_wt_his*gr;
	curAngle_s= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_r= 0;//ת�ӽ�,�ȸ�Ϊ0����ʼ��ʱ�����¸�ֵ
	curAngle_wt=0;//�ۺϺ�Ƕ�
	//curAngle_wt=curAngle_wt*gr;
	/*double rat=v/(w_wt*w/gr*0.44704);
	Cp=0.5*(rat-0.022*pitch_angle*pitch_angle-5.6)*exp(-0.17*rat);
	Pw=0.5*Cp*(PI*R*R)*(v*v*v)*p_air/Sbase;
	Tw=Pw/w_wt*/;
	Tw=0.8;


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
	}

	/************************************************************************/
	/*              ��EMTP�ӿ�����ı���                                    */
	/************************************************************************/
	for (int i=0;i<nPort;i++)
	{
		nodeNumber[i] = i+firstNode;
	}

	/************************************************************************/
	/*              ����ϵ�������ŵ�ٵ�ֵ����                              */
	/************************************************************************/
	calculateCoefficientMatrix();
}
WindInducGenerator::~WindInducGenerator()
{//��������

}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/

void WindInducGenerator::initializeBranch(TVectorD &initialVoltageArray,TVectorD &initialCurrentArray,int& ptr, double time)
{//��ʼ��֧·��ѹ����
	//��ʼ��wr��curAngle
    //��ʼ��**************************//
	//?????
	//**********************************

	readNodeVoltage(initialVoltageArray);
	calculateBranchVoltage();

	for (int i=0;i<3;i++)
	{
		Isabc[i] = initialCurrentArray[ptr+i];
	}
	ptr+=3;

	for (int i=0;i<3;i++)
	{
		Irabc[i] = initialCurrentArray[ptr+i];
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

void WindInducGenerator::readNodeVoltage(TVectorD &nodeVoltageArray)
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
void WindInducGenerator::calculateBranchVoltage()
{//����֧·��ѹ
}

void WindInducGenerator::calculateBranchCurrent()
{//����֧·����
	for (int i=0;i<3;i++)
	{
		Isabc_1[i] = Isabc[i];
		Isabc[i] = Vsabc[i]*nortonEquivalentConductance[0]+nortonEquivalentCurrent[i];
		Irabc_1[i] = Irabc[i];
		Irabc[i] = Vrabc[i]*nortonEquivalentConductance[1]+nortonEquivalentCurrent[i+3];
	}
}

void WindInducGenerator::calculateNortonEquivalentCurrent(double time)
{//����֧·��ŵ�ٵ�Ч��·�еĵ�����
	calculateCoefficientMatrix();
	calculateDq0Results();//����dq0����ϵ�µĽ����Ԥ�����ӵ���
	calculateOmega(time);

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

void WindInducGenerator::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{//�γɽڵ�ŵ�ٵ�Ч��������
	int N;
	for (int i=0;i<6;i++)
	{
		N = nodeNumber[i];
		nodeNortonEquivalentCurrentArray(N) -= nortonEquivalentCurrent[i];
	}
}
void WindInducGenerator::formConductanceMatrix(TMatrixD& conductanceMatrix)
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

void WindInducGenerator::saveBranchCurrent(TMatrixD &branchCurrentMatrix,int& ptr,int counter)
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

void WindInducGenerator::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
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
void WindInducGenerator::interpolate(double ratio)
{//������ֵ����֧·�ĵ�ѹ�������в�ֵ

	//��ѹ����������ֵ
	for (int k=0;k<3;k++)
	{
		Vsabc[k] = (1-ratio)*Vsabc_1[k] + ratio*Vsabc[k];
		Isabc[k] = (1-ratio)*Isabc_1[k] + ratio*Isabc[k];
		Vrabc[k] = (1-ratio)*Vrabc_1[k] + ratio*Vrabc[k];
		Irabc[k] = (1-ratio)*Irabc_1[k] + ratio*Irabc[k];
	}

	//ŵ�ٵ�ֵ������ֵ
	double temp1,temp2;
	for (int k=0;k<6;k++)
	{
		temp1 = (1-ratio)*nortonEquivalentCurrent_1[k] + ratio*nortonEquivalentCurrent[k];
		temp2 = (-ratio)*nortonEquivalentCurrent_1[k] + (1+ratio)*nortonEquivalentCurrent[k];

		nortonEquivalentCurrent_1[k] = temp1;
		nortonEquivalentCurrent[k] = temp2;
	}
}

void WindInducGenerator::updateResult(int updateTypeNum)
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

void WindInducGenerator::updateConductanceMatrix(TMatrixD & conductanceMatrix)
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
void WindInducGenerator::calculateCoefficientMatrix()
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
	nortonEquivalentConductance[0] /= Zbase;
	nortonEquivalentConductance[1] = nortonEquivalentConductance[1]/Zbase*turnratio*turnratio;

	//����Gn1...Gn3
	Gn1 = (-1.0)*inv_Aa*AA_res;
	Gn2 = (-1.0)*inv_Aa*BB;
	Gn3 = coff*inv_Aa;
}

//����dq0����ϵ�¶��Ӻ�ת�ӵĵ�ѹ��������Ԥ����һʱ�̶��ӵ���
//Ϊ��һʱ������ŵ�ٵ�Ч������׼��
void WindInducGenerator::calculateDq0Results()
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
	}
}

//����е���̣�����ת��ת��
//Ϊ��һʱ�����ϵ�����󡢼���ŵ�ٵ�ֵ������׼��
void WindInducGenerator::calculateOmega(double time)
{
	calculateTw(time);
	wr_his = wr;
	w_wt_his=w_wt;
	curAngle_r_his=curAngle_r;
	curAngle_s_his=curAngle_s;
	curAngle_wt_his=curAngle_wt;
	if (control==0)
	{
		curAngle_s=curAngle_s_his+ deltaT*w;
		bo(1,1)=2*H_g/deltaT*wr_his+k/2*curAngle_wt_his-k/2*curAngle_s_his+k/2*curAngle_r_his-Xm*(Irdq0[0]*Isdq0[1]-Irdq0[1]*Isdq0[0])-D_g/2*wr_his-k/2*curAngle_s;
		//bo(1,1)=2*H_g/deltaT*wr_his+k/2*curAngle_wt_his-k/2*curAngle_s_his+k/2*curAngle_r_his-0.1-D_g/2*wr_his-k/2*curAngle_s;
		bo(2,1)=2*H_wt/deltaT*w_wt_his+Tw-D_wt/2*w_wt_his-k/2*curAngle_wt_his+k/2*curAngle_s_his-k/2*curAngle_r_his+k/2*curAngle_s;
		bo(3,1)=curAngle_r_his+w*deltaT-w*deltaT/2*wr_his;
		bo(4,1)=curAngle_wt_his+w*deltaT/2*w_wt_his;
		co=Co_inv*bo;
		wr=co(1,1);
		w_wt=co(2,1);
		curAngle_r=co(3,1);
		curAngle_wt=co(4,1);
	}
	if (control==1)
	{
		curAngle_r = curAngle_r + deltaT*(1-wr)*w;
		curAngle_s=curAngle_s+ deltaT*w;
	}

/*	if (curAngle_s>2*PI)
	{
		curAngle_s -= 2*PI;
	}
	if (curAngle_r>2*PI)
	{
		curAngle_r -= 2*PI;
	}
	if (curAngle_r<0)
	{
		curAngle_r += 2*PI;
	}
	if (curAngle_wt>2*PI)
	{
		curAngle_wt -= 2*PI;
	}*/
}

//�����ɿ˱任�뷴�任����
void WindInducGenerator::parkTransMatrix(TMatrixD &T,double angle)
{
	T(1,1) = 2.0/3.0*cos(angle);
	T(1,2) = 2.0/3.0*cos(angle-2.0/3.0*PI);
	T(1,3) = 2.0/3.0*cos(angle+2.0/3.0*PI);
	T(2,1) = -2.0/3.0*sin(angle);
	T(2,2) = -2.0/3.0*sin(angle-2.0/3.0*PI);
	T(2,3) = -2.0/3.0*sin(angle+2.0/3.0*PI);
}
void WindInducGenerator::invParkTransMatrix(TMatrixD &T_,double angle)
{
	T_(1,1) = cos(angle);
	T_(1,2) = -sin(angle);
	T_(2,1) = cos(angle-2.0/3.0*PI);
	T_(2,2) = -sin(angle-2.0/3.0*PI);
	T_(3,1) = cos(angle+2.0/3.0*PI);
	T_(3,2) = -sin(angle+2.0/3.0*PI);
}
void WindInducGenerator::calculateTw(double time)
{
	if (time>1.0)
	{
		double rat=v/(w_wt*w/gr*0.44704);
		Cp=0.5*(rat-0.022*pitch_angle*pitch_angle-5.6)*exp(-0.17*rat);
		Pw=0.5*Cp*(PI*R*R)*(v*v*v)*p_air/Sbase;
		Tw=Pw/w_wt;
	} 
	else
	{
		Tw=0.8;
	}
}