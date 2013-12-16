#include "PWMConverter.h"

#include <cmath>
#include <iostream>
using namespace std;

#define PI 3.141592653589793238462643383279

/************************************************************************/
/*                ���캯������������                                    */
/************************************************************************/
// ���캯��
//�����ڵ�˳��: a,b,c,dc1,dc2�����ƽڵ�˳��:a,b,c,tri
PWMConverter::PWMConverter(int firstACNode,int firstDCNode, double R, double L, double C, int firstCtrlNode, int subType)
{
	this->subType = subType;

	isUserDef = 0; // �Ƿ�μӲ�ֵ��Ӧ��ƽ��ģ�͵�ʱ�򲻱�Ҫ
	need_NEC=1; // ��Ҫ����ŵ�ٵ�ֵ����

	nPort = 5;
	nodeNumber = new int[nPort];
	type = 24; // ���

	// �����˿ڱ��
	nodeNumber[0] = firstACNode;
	nodeNumber[1] = firstACNode + 1;
	nodeNumber[2] = firstACNode + 2;
	nodeNumber[3] = firstDCNode;
	nodeNumber[4] = firstDCNode + 1;

	for (int i=0; i<4; i++)
	{
		ctrlNodeNumber[i] = i+firstCtrlNode;
	}

	// Ԫ������
	this->R = R + 0.001; //���ǿ��ص�ͨ����
	this->L = L;
	this->C = C;

	// ŵ�ٵ�ֵ���
	// ���ɾ���
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Y11[i][j] = 0;

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Y12[i][j] = 0;

	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Y21[i][j] = 0;

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			Y22[i][j] = 0;

	for (int i=0; i<5; i++)
		for (int j=0; j<5; j++)
			Yne[i][j] = 0;

	// ��ֵ����
	for (int k=0; k<3; k++)
		Ine1[k] = 0;

	for (int k=0; k<2; k++)
		Ine2[k] = 0;

	for (int k=0; k<5; k++)
		Ine[k] = 0;

	// �˿ڵ�ѹ����
	// ��ѹ
	for (int k=0; k<3; k++)
		Ub1[k] = 0;

	for (int k=0; k<2; k++)
		Ub2[k] = 0;

	for (int k=0; k<5; k++)
		Ub[k] = 0;

	// ����
	for (int k=0; k<3; k++)
		Ib1[k] = 0;

	for (int k=0; k<2; k++)
		Ib2[k] = 0;

	for (int k=0; k<5; k++)
		Ib[k] = 0;

	// ƽ��ģ�Ͳ���
	for (int k = 0; k < 3; k++) {
		tao_s1[k] = 0;
		tao_s2[k] = 0;
		tao_s[k] = 0;
	}

	for (int k = 0; k < 3; k++) {
		tao_s1_1[k] = 0;
		tao_s2_1[k] = 0;
	}

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			D[i][j] = 0;
}

// ��������
PWMConverter::~PWMConverter()
{
}

/************************************************************************/
/*                ��EMTP�ӿڵĺ���                                      */
/************************************************************************/
//��ʼ��֧·��ѹ����
void PWMConverter::initializeBranch(TVectorD& initialVoltageArray,TVectorD& initialCurrentArray,int& ptr, double time)
{
	// ��ʼ��֧·��ѹ�͵���
	// ��ѹ
	readNodeVoltage(initialVoltageArray);
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub[k];
	}
	for (int k=0; k<2; k++) {
		Ub2[k] = Ub[k+3];
	}

	// ����
	for (int k=0; k<5; k++) {
		Ib[k] = initialCurrentArray(ptr+k);
	}
	ptr += 5;
	for (int k=0; k<3; k++) {
		Ib1[k] = Ib[k];
	}
	for (int k=0; k<2; k++) {
		Ib2[k] = Ib[k+3];
	}
}

//�ӽڵ��ѹ�����ж�֧·���ڵ��ѹ
void PWMConverter::readNodeVoltage(TVectorD& nodeVoltageArray)
{
	for (int i=0;i<5;i++)
	{
		if (nodeNumber[i]==0)
		{
			Ub[i] = 0;
		}
		else
		{
			Ub[i] = nodeVoltageArray(nodeNumber[i]);
		}	
	}
}

// ����֧·��ѹ
void PWMConverter::calculateBranchVoltage()
{
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub[k];
	}
	for (int k=0; k<2; k++) {
		Ub2[k] = Ub[k+3];
	}
}

// ����֧·����
void PWMConverter::calculateBranchCurrent()
{
	// Ib1 = Y11*Ub1 + Y12*Ub2 + Ine1
	for (int i=0; i<3; i++) {
		Ib1[i] = Ine1[i];
		for (int j=0; j<3; j++)
			Ib1[i] += Y11[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ib1[i] += Y12[i][j]*Ub2[j];
	}

	// Ib2 = Y21*Ub1 + Y22*Ub2 + Ine2
	for (int i=0; i<2; i++) {
		Ib2[i] = Ine2[i];
		for (int j=0; j<3; j++)
			Ib2[i] += Y21[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ib2[i] += Y22[i][j]*Ub2[j];
	}

	// Ib
	for (int k=0; k<3; k++) {
		Ib[k] = Ib1[k];
	}

	for (int k=0; k<2; k++) {
		Ib[k+3] = Ib2[k];
	}

	//// debug
	//for (int k=0; k<5; k++) {
	//	cout << Ib[k] << "  ";
	//}
	//cout << endl;
}

// ����Ԫ��ŵ�ٵ�ֵ����
void PWMConverter::calculateNortonEquivalentCurrent(double time)
{
	// Ine1 = Y11*Ub1 + Y12*Ub2 + (1-(deltaT/2)*(R/L))/(1+(deltaT/2)*(R/L))*Ib1
	double tmp = (1-(deltaT/2)*(R/L))/(1+(deltaT/2)*(R/L));
	for (int i=0; i<3; i++) {
		Ine1[i] = tmp*Ib1[i];
		for (int j=0; j<3; j++)
			Ine1[i] += Y11[i][j]*Ub1[j];
		for (int j=0; j<2; j++)
			Ine1[i] += Y12[i][j]*Ub2[j];
	}

	if ( subType == 1 )
	{
		// Ine2 = - DT*Ine1 - (2*C/deltaT)*Ub2 - DT*Ib1 - Ib2
		tmp = -2*C/deltaT;
		for (int i=0; i<2; i++) {
			Ine2[i] = tmp*Ub2[i] - Ib2[i];
			for (int j=0; j<3; j++)
				Ine2[i] -= D[j][i]*(Ine1[j]+Ib1[j]);
		}
	}
	else
	{
		// Ine2 = - DT*Ine1 - (2*C/deltaT)*M2*Ub2 - DT*Ib1 - Ib2
		// M2 = [1,-1;-1,1]
		tmp = -2*C/deltaT;
		Ine2[0] = tmp*(Ub2[0]-Ub2[1]);
		Ine2[1] = tmp*(Ub2[1]-Ub2[0]);
		for (int i=0; i<2; i++) {
			Ine2[i] -= Ib2[i];
			for (int j=0; j<3; j++)
				Ine2[i] -= D[j][i]*(Ine1[j]+Ib1[j]);
		}
	}

	// Ine
	for (int k=0; k<3; k++) {
		Ine[k] = Ine1[k];
	}

	for (int k=0; k<2; k++) {
		Ine[k+3] = Ine2[k];
	}
}

// �γɽڵ�ŵ�ٵ�Ч��������
void PWMConverter::formNodeNortonEquivalentCurrentArray(TVectorD &nodeNortonEquivalentCurrentArray)
{
	int N;
	for (int k=0; k<5; k++)
	{
		N = nodeNumber[k];
		if ( N != 0 ) { 
			nodeNortonEquivalentCurrentArray(N) -= Ine[k];
		}
	}
}

// �γɽڵ㵼����
void PWMConverter::formConductanceMatrix(TMatrixD &conductanceMatrix)
{
	// û�п��Ƕ˿ڽӵص����
	int N1, N2;
	double tmpd;
	for (int i=0; i<5; i++)
	{
		N1 = nodeNumber[i];
		for (int j=0; j<5; j++)
		{
			N2 = nodeNumber[j];
			tmpd = conductanceMatrix(N1,N2);
			conductanceMatrix(N1,N2) = tmpd + Yne[i][j];
		}
	}
}

// ����֧·����
void PWMConverter::saveBranchCurrent(double** branchCurrentMatrix_1,int& ptr,int counter)
{
	for (int k=0; k<5; k++) {
		branchCurrentMatrix_1[counter-1][ptr+k] = Ib[k];
	}
	ptr += 5;
}

// ����ڵ㵼�ɾ���
void PWMConverter::calculateYne()
{
	// Y11
	Y11[0][0] = deltaT/(2*L)/(1+deltaT/2*R/L);
	Y11[1][1] = Y11[0][0];
	Y11[2][2] = Y11[0][0];

	// Y12
	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Y12[i][j] = -Y11[0][0]*D[i][j];

	// Y21
	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Y21[i][j] = Y12[j][i];

	// Y22
	Y22[0][0] = 2*C/deltaT;
	Y22[1][1] = Y22[0][0];
	if ( subType == 1 ) // xuyin,20121221,subType=1��ʾֱ�����е�ӵأ�subType=2�򲻽ӵ�
	{
		Y22[0][1] = 0;
		Y22[1][0] = 0;
	}
	else
	{
		Y22[0][1] = -Y22[0][0];
		Y22[1][0] = -Y22[0][0];
	}

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			for (int k=0; k<3; k++)
				Y22[i][j] += Y11[0][0]*D[k][i]*D[k][j];

	// ��Y11,Y12,Y21,Y22�õ�Yne
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			Yne[i][j] = Y11[i][j];

	for (int i=0; i<3; i++)
		for (int j=0; j<2; j++)
			Yne[i][j+3] = Y12[i][j];

	for (int i=0; i<2; i++)
		for (int j=0; j<3; j++)
			Yne[i+3][j] = Y21[i][j];

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++)
			Yne[i+3][j+3] = Y22[i][j];
}

#ifndef GENERAL
// ��ʼ������ʱ��,version 1, ��ͨ�ð�
void PWMConverter::initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr)
{
	// PWMInitialMatrix��ÿ3��һ�飬��Ӧһ��PWM��������Ϊa,b,c������һ�����������ڵĿ��������ź�
	// nStep��ʾÿ���������ڰ�����ʱ��������ΪPWMInitialMatrix��������1

	// ��ʼ������ʱ��tao_s1,tao_s2,tao_s
	for (int k=0; k<3; k++) {
		tao_s1[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMInitialMatrix[i][ptr+k] == 1) {
				tao_s1[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMInitialMatrix[i][ptr+k] == 0) {
				tao_s2[k] = i/(double)nStep;
				break;
			}
		}

		tao_s[k] = tao_s2[k] - tao_s1[k];
	}
	ptr += 3;

	for (int k=0; k<3; k++) {
		tao_s1_1[k] = tao_s1[k];
		tao_s2_1[k] = tao_s2[k];
	}

	// �������D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}

	//// debug
	//for (int k=0; k<3; k++) {
	//	cout<<tao_s1[k]<<"\t"<<tao_s2[k]<<"\t"<<tao_s[k]<<endl;
	//}

	// ����ŵ�ٵ�ֵ���ɾ���
	calculateYne();
}

#else
// ��ʼ������ʱ��,version 2, ͨ�ð�, xuyin 20121211
void PWMConverter::initializeSwitchingInstants(double** PWMInitialMatrix, int nStep, int& ptr)
{
	// PWMInitialMatrix��ÿ3��һ�飬��Ӧһ��PWM��������Ϊa,b,c������һ�����������ڵĿ��������ź�
	// nStep��ʾÿ���������ڰ�����ʱ��������ΪPWMInitialMatrix��������1

	// ��ʼ��tao_s
	for (int k=0; k<3; k++) {
		tao_s[k] = 0.0;
		for (int i=1; i<=nStep; i++) {
			if ( PWMInitialMatrix[i-1][ptr+k] == 1 && PWMInitialMatrix[i][ptr+k] == 1 )
			{
					tao_s[k] += 1;
			}
			else if ( PWMInitialMatrix[i-1][ptr+k] == 0 && PWMInitialMatrix[i][ptr+k] == 0 )
			{
					// tao_s_new[k] += 0;
			}
			else
			{
				tao_s[k] += 0.5;
			}
		}
		tao_s[k] /= nStep;
	}
	ptr += 3;

	// �������D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}

	// ����ŵ�ٵ�ֵ���ɾ���
	calculateYne();
}
#endif

void PWMConverter::ctrlNodeNumber_PWM(int* ctrlNodeNumber, int k)
{
	for (int i=0; i<4; i++)
	{
		ctrlNodeNumber[4*k+i] = this->ctrlNodeNumber[i];
	}
}

// Ԥ�⿪��ʱ��1���������
void PWMConverter::predictSwithcingInstants()
{
	double tao_s1_new[3];
	double tao_s2_new[3];

	for (int k=0; k<3; k++) {
		tao_s1_new[k] = 2*tao_s1[k] - tao_s1_1[k];
		tao_s2_new[k] = 2*tao_s2[k] - tao_s2_1[k];
	}

	for (int k=0; k<3; k++) {
		tao_s1_1[k] = tao_s1[k];
		tao_s2_1[k] = tao_s2[k];
	}

	for (int k=0; k<3; k++) {
		tao_s1[k] = tao_s1_new[k];
		tao_s2[k] = tao_s2_new[k];
	}

	// �������D
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}
	
	// ����ŵ�ٵ�ֵ���ɾ���
	calculateYne();
}

// Ԥ�⿪��ʱ��2������������ϵͳ
int PWMConverter::predictSwithcingInstants(double ** PWMPredictMatrix, int nStep, int& ptr)
{
	// PWMPredictMatrix��ÿ4��һ�飬��Ӧһ��PWM������
	// ǰ3��Ϊa,b,c����ĵ��Ʋ��źţ���4��Ϊ�ز��ź�
	// nStep��ʾÿ���������ڰ�����ʱ��������ΪPWMPredictMatrix��������1

	// Ԥ��tao_s1,tao_s2,tao_s
	double temp1, temp2;
	double tao_s1_pre[3];
	double tao_s2_pre[3];
	for (int k=0; k<3; k++) {
		tao_s1_pre[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMPredictMatrix[i][ptr+k] >= PWMPredictMatrix[i][ptr+3]) {
				if (i > 0) {
					temp1 = fabs(PWMPredictMatrix[i][ptr+k]-PWMPredictMatrix[i][ptr+3]);
					temp2 = fabs(PWMPredictMatrix[i-1][ptr+k]-PWMPredictMatrix[i-1][ptr+3]);
					tao_s1_pre[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s1_pre[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2_pre[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMPredictMatrix[i][ptr+k] <= PWMPredictMatrix[i][ptr+3]) {
				if (i > nStep/2) {
					temp1 = fabs(PWMPredictMatrix[i][ptr+k]-PWMPredictMatrix[i][ptr+3]);
					temp2 = fabs(PWMPredictMatrix[i-1][ptr+k]-PWMPredictMatrix[i-1][ptr+3]);
					tao_s2_pre[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s2_pre[k] = i/(double)nStep;
				break;
			}
		}
	}
	ptr += 4;

	// �Ƚ�ǰ������Ԥ��ֵ
	double maxErr = 0.0;
	for (int k=0; k<3; k++) {
		if ( fabs(tao_s1[k]-tao_s1_pre[k]) > maxErr )
			maxErr = fabs(tao_s1[k]-tao_s1_pre[k]);
		if ( fabs(tao_s2[k]-tao_s2_pre[k]) > maxErr )
			maxErr = fabs(tao_s2[k]-tao_s2_pre[k]);
	}

	if ( maxErr < 0.005 )
	{
		return 0;
	}
	else
	{
		for (int k=0; k<3; k++) {
			tao_s1[k] = tao_s1_pre[k];
			tao_s2[k] = tao_s2_pre[k];
			tao_s[k] = tao_s2[k] - tao_s1[k];
		}

		// �������D
		for (int k=0; k<3; k++) {
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}
	
		// ����ŵ�ٵ�ֵ���ɾ���
		calculateYne();

		return 1;
	}
}

#ifndef GENERAL
// У������ʱ��version1,��ͨ�ð�
// ����������ʼʱ��Ϊ�ز����ֵʱ��
// ƽ��ʱ��γ�����Ҫǡ�õ���һ����������
int PWMConverter::correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr)
{
	// PWMCorrectMatrix��ÿ4��һ�飬��Ӧһ��PWM������
	// ǰ3��Ϊa,b,c����ĵ��Ʋ��źţ���4��Ϊ�ز��ź�
	// nStep��ʾÿ���������ڰ�����ʱ��������ΪPWMCorrectMatrix��������1

	// debug
	//for (int k=0; k<=nStep; k++){
	//	cout << k+1 << "\t";
	//	for (int j=0; j<4; j++)
	//		cout << PWMCorrectMatrix[k][j] << "\t";
	//	cout << endl;
	//}
	//cout << endl;

	// �������õ��Ŀ���ʱ��,version1,��ͨ�ð�
	// ����������ʼʱ��Ϊ�ز����ֵʱ��
	// ƽ��ʱ��γ�����Ҫǡ�õ���һ����������
	double temp1, temp2;
	double tao_s1_new[3], tao_s2_new[3];
	for (int k=0; k<3; k++) {
		tao_s1_new[k] = 0.5;
		for (int i=0; i<=nStep/2; i++) {
			if (PWMCorrectMatrix[i][ptr+k] >= PWMCorrectMatrix[i][ptr+3]) {
				if (i > 0) {
					temp1 = fabs(PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3]);
					temp2 = fabs(PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3]);
					tao_s1_new[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s1_new[k] = i/(double)nStep;
				break;
			}
		}

		tao_s2_new[k] = 1.0;
		for (int i=nStep/2; i<=nStep; i++) {
			if (PWMCorrectMatrix[i][ptr+k] <= PWMCorrectMatrix[i][ptr+3]) {
				if (i > nStep/2) {
					temp1 = fabs(PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3]);
					temp2 = fabs(PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3]);
					tao_s2_new[k] = i/(double)nStep - temp1/(temp1+temp2)/(double)nStep;
				}
				else
					tao_s2_new[k] = i/(double)nStep;
				break;
			}
		}
	}
	ptr += 4;

	//// debug
	//cout << "tao_s1:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s1[k] << "  ";
	//cout << endl;

	//cout << "tao_s2:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s2[k] << "  ";
	//cout << endl;

	//cout << "tao_s1_new:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s1_new[k] << "  ";
	//cout << endl;

	//cout << "tao_s2_new:  ";
	//for (int k=0; k<3; k++)
	//	cout << tao_s2_new[k] << "  ";
	//cout << endl << endl;

	// ����������
	double maxError = 0;
	for (int k=0; k<3; k++) {
		if ( fabs(tao_s1[k]-tao_s1_new[k]) > maxError )
			maxError = fabs(tao_s1[k]-tao_s1_new[k]);
		if ( fabs(tao_s2[k]-tao_s2_new[k]) > maxError )
			maxError = fabs(tao_s2[k]-tao_s2_new[k]);
	}

	// У��
	if ( maxError >= tol ) // �������趨ֵʱ��ҪУ��
	{
		for (int k=0; k<3; k++) {
			tao_s1[k] = tao_s1_new[k];
			tao_s2[k] = tao_s2_new[k];
			tao_s[k] = tao_s2[k] - tao_s1[k];
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}		

		calculateYne();

		return 1;
	}
	else
		return 0;
}

#else
// У������ʱ��version2, ͨ�ð�, xuyin 20121211
int PWMConverter::correctSwithcingInstants(double** PWMCorrectMatrix, double tol, int nStep, int& ptr)
{
	// PWMCorrectMatrix��ÿ4��һ�飬��Ӧһ��PWM������
	// ǰ3��Ϊa,b,c����ĵ��Ʋ��źţ���4��Ϊ�ز��ź�
	// nStep��ʾÿ���������ڰ�����ʱ��������ΪPWMCorrectMatrix��������1

	//cout << "PWMCorrectMatrix:" << endl;
	//for (int i=0; i<=nStep; i++) {
	//	for (int j=0; j<4; j++)
	//		cout << PWMCorrectMatrix[i][j] << "   ";
	//	cout << endl;
	//}

	// �������õ��Ŀ���ʱ��,version2,ͨ�ð�
	double tao_s_new[3];
	for (int k=0; k<3; k++)
	{
		tao_s_new[k] = 0.0;
		for (int i=1; i<=nStep; i++)
		{
			if (PWMCorrectMatrix[i-1][ptr+k] >= PWMCorrectMatrix[i-1][ptr+3]
				&& PWMCorrectMatrix[i][ptr+k] >= PWMCorrectMatrix[i][ptr+3]) 
			{
				tao_s_new[k] += 1;
			}
			else if (PWMCorrectMatrix[i-1][ptr+k] < PWMCorrectMatrix[i-1][ptr+3]
				&& PWMCorrectMatrix[i][ptr+k] < PWMCorrectMatrix[i][ptr+3]) 
			{
				// tao_s_new[k] += 0;
			}
			else
			{
				double temp1 = PWMCorrectMatrix[i][ptr+k]-PWMCorrectMatrix[i][ptr+3];
				double temp2 = PWMCorrectMatrix[i-1][ptr+k]-PWMCorrectMatrix[i-1][ptr+3];
				tao_s_new[k] += 0.5*(1+(temp1+temp2)/(fabs(temp1)+fabs(temp2)));
			}
		}
		tao_s_new[k] /= nStep;
	}

	ptr += 4;

	// ����������
	double maxError = 0;
	for (int k=0; k<3; k++)
		if ( fabs(tao_s[k]-tao_s_new[k]) > maxError )
			maxError = fabs(tao_s[k]-tao_s_new[k]);

	if (ptr == 4)
		int mydebug = 0;

	// У��
	if ( maxError >= tol ) // �������趨ֵʱ��ҪУ��
	{
		for (int k=0; k<3; k++) {
			tao_s[k] = tao_s[k] + 0.4*(tao_s_new[k]-tao_s[k]);
			D[k][0] = tao_s[k];
			D[k][1] = 1 - tao_s[k];
		}		

		calculateYne();

		return 1;
	}
	else
	{
		return 0;
	}
}
#endif

// �洢�ڲ��������Ա��ڵ���У��ʱ�ָ�
void PWMConverter::storeInternalVariables()
{
	// ����PWM�����������ݵ�ѹ��������
	for (int k=0; k<3; k++) {
		Ub1_bak[k] = Ub1[k];
		Ib1_bak[k] = Ib1[k];
	}

	for (int k=0; k<2; k++) {
		Ub2_bak[k] = Ub2[k];
		Ib2_bak[k] = Ib2[k];
	}

	for (int k=0; k<5; k++) {
		Ub_bak[k] = Ub[k];
		Ib_bak[k] = Ib[k];
	}
}

// ����У��ʱ�ָ��ڲ�����
void PWMConverter::restoreInternalVariables()
{
	for (int k=0; k<3; k++) {
		Ub1[k] = Ub1_bak[k];
		Ib1[k] = Ib1_bak[k];
	}

	for (int k=0; k<2; k++) {
		Ub2[k] = Ub2_bak[k];
		Ib2[k] = Ib2_bak[k];
	}

	for (int k=0; k<5; k++) {
		Ub[k] = Ub_bak[k];
		Ib[k] = Ib_bak[k];
	}
}

void  PWMConverter::correctYne(double* tao_s)
{
	for (int k=0; k<3; k++) {
		D[k][0] = tao_s[k];
		D[k][1] = 1 - tao_s[k];
	}		

	calculateYne();
}