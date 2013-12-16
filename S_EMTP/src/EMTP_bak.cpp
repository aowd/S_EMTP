#include "EMTP.h"
#include <set>
#include <sstream>

EMTP::EMTP(){
	branches = new std::vector<Component*>();
	ctrlBranches = new std::vector<CtrlComponent*>();
	msrComps = new std::vector<MeasureComponent*>();
}

EMTP::~EMTP(){

	if (branches != NULL){
		for (int i = 0; i < branches->size(); i ++){
			Component * tempCom = branches->at(i);
			if (tempCom != NULL){
				delete tempCom;
			}
		}
	}
	branches->clear();
	delete  branches;

	if (nSwitch != 0)
	{
		delete[] diodeNumArray;
		delete[] switchNumArray;
		delete[] switchRatioArray;
		delete[] switchModeArray;
	}

	if (nInverter != 0)
	{
		delete[] inverterNumArray;
	}

	if (nBranch_NEC != 0)
	{
		delete[] branchNumArray_NEC;
	}

	if (nBranch_userDef != 0)
	{
		delete[] branchNumArray_userDef;
	}

	if (ctrlBranches != NULL){
		for (int i = 0; i < ctrlBranches->size(); i ++){
			CtrlComponent * tempCtrlCom = ctrlBranches->at(i);
			if (tempCtrlCom != NULL){
				delete tempCtrlCom;
			}
		}
	}
	ctrlBranches->clear();
	delete  ctrlBranches;

	if (nCtrlBranches!=0)
	{
		delete[] branchCalMark;
		delete[] nodeCalMark;
		delete[] ctrlNodeValue;
	}
	if (msrComps != NULL){
		for (int i = 0; i < msrComps->size(); i ++){
			MeasureComponent * tempMsrCom = msrComps->at(i);
			if (tempMsrCom != NULL){
				delete tempMsrCom;
			}
		}
	}
	msrComps->clear();
	delete  msrComps;

	if (branchNumArray_PWMConverter != NULL) {
		delete[] branchNumArray_PWMConverter;
		delete[] ctrlNumArray_PWMConverter;
	}

	if (branchNumArray_IndGen != NULL) {
		delete[] branchNumArray_IndGen;
	}
}

void EMTP::initializeSystem() // ��ʼ��������0�������PSCAD����
{
/////////////////////////////////////////
// edited by xuyin, 20130119
// ����INITIAL_VIʱ���ô������PSCAD�����ļ����г�ʼ��
// �����������ļ���2�����ֱ�洢��ѹ����������Ϊcase***_V.dat��case***_I.dat
/////////////////////////////////////////
#ifdef INITIAL_VI
	TVectorD initialVoltageVec(1,nNodes), initialCurrentVec(1,nColumns);

	if (initializeOption==0)
	{
		initialVoltageVec *= 0.0;
		initialCurrentVec *= 0.0;
		cout<<"S_EMTP begin with 0 state"<<endl;
	} 
	else
	{
		cout<<"S_EMTP begin with PSCAD state"<<endl;

		// ��ȡ��ѹ��ֵ�ļ�
		double t; // ��һ��Ϊʱ��
		char chDir[100];   
		string caseDir;
#ifdef WIN32
		sprintf(chDir,".\\result\\pscadData\\case%d_V.dat",caseNumber);   
#else
		sprintf(chDir,"..//result//pscadData//case%d_V.dat",caseNumber);   
#endif
		caseDir=chDir;
		ifstream infile(caseDir.c_str(),ios::in);

		// ���˵�startTime��ǰ����
		infile >> t;
		double tmp = 0;
		while ( fabs(t-startTime) > initDeltaT/2 ) {
			for(int k=0;k<nNodes;k++)
				infile >> tmp;
			infile >> t;
		}

		// �γɵ�ѹ��ֵ����
		for(int k=1;k<=nNodes;k++)
			infile>>initialVoltageVec(k);

		infile.close();
		infile.clear();

		// ��ȡ������ֵ�ļ�
#ifdef WIN32
		sprintf(chDir,".\\result\\pscadData\\case%d_I.dat",caseNumber);   
#else
		sprintf(chDir,"..//result//pscadData//case%d_I.dat",caseNumber);   
#endif
		caseDir=chDir;
		infile.open(caseDir.c_str(),ios::in);

		// ���˵�startTime��ǰ����
		infile >> t;
		while ( fabs(t-startTime) > initDeltaT/2 ) {
			for(int k=0;k<nColumns;k++)
				infile >> tmp;
			infile >> t;
		}

		// �γɵ�����ֵ����
		for(int k=1;k<=nColumns;k++)
			infile>>initialCurrentVec(k);

		infile.close();
		infile.clear();

		initialVoltageVec *= 1000.0;
		initialCurrentVec *= 1000.0;
	}

	// ��ʼ���ڵ��ѹ
	nodeVoltageVec = initialVoltageVec;

	// ��ʼ��֧·
	int ptr = 1;//ָ��
	for (int i = 1; i <= nBranch; i ++)
	{
		Component * tempCom = branches->at(i-1);
		tempCom->initializeBranch(initialVoltageVec,initialCurrentVec,ptr);
	}

	// �洢��ʼ�����
	for (int i=1; i<=nNodes; i++)
		nodeVoltageMatrix_1[0][i-1] = initialVoltageVec(i);

	for (int i=1; i<=nColumns; i++)
		branchCurrentMatrix_1[0][i-1] = initialCurrentVec(i);

#else
	int nCol=0;//��ֵ������ά��
	nCol = nNodes + nColumns;

	int units=0;//��ֵ����ά���ĸ�λ��
	units=nCol%10;

	int nFiles=0;//PSCAD.out�ļ��ĸ�����ÿ��,out�ļ��д���11��������һ����ʱ�䣬����10����10������������
	nFiles=(nCol-units)/10;

	int nVecIndex=1;
	TVectorD initiaValVec,branchCurrentVec;
	initiaValVec.ResizeTo(1,nCol);
	branchCurrentVec.ResizeTo(1,nColumns);

	if(initializeOption==0)//���option=0
	{
		initiaValVec*=0.0;
		cout<<"S_EMTP begin with 0 state"<<endl;
	}
	else
	{
		cout<<"S_EMTP begin with PSCAD state"<<endl;
		//	step1:cal number of *.out files 
		if(units!=0)//�����λ��Ϊ0
		{
			nFiles=nFiles+1;
		}

		//	step2:scan all *.out files and read init data
		for (int fileIndex=1;fileIndex<=nFiles;fileIndex++)
		{
			int num1;
			int num2;
			double t;//��һ��Ϊʱ��
			char chDir[100];   
			string caseDir;
			num2=fileIndex%10;
			num1=(fileIndex-num2)/10;
#ifdef WIN32
			sprintf(chDir,".\\result\\pscadData\\case%d_%d%d.out",caseNumber,num1,num2);   
#else
			sprintf(chDir,"..//result//pscadData//case%d_%d%d.out",caseNumber,num1,num2);   
#endif
			caseDir=chDir;

			ifstream infile(caseDir.c_str(),ios::in);//open cur file stream

			////////////////////////////////////////////////////////////////////////////////////
			// edited by xuyin 20121109
			// ���˵�startTime��ǰ����
			if (fileIndex<nFiles)//�������һ���ļ���ֱ�Ӷ��ߵ�һ�е�10�����ݣ���ȥʱ�������
			{
				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}
				
				for(int k=0;k<10;k++)//read 10 initdata
				{
					infile>>initiaValVec(nVecIndex);
					nVecIndex++;
				}
			}
			else//�������һ���ļ���ֻ��Ҫ����units����������
			{
				if (units==0)
				{
					units=10;
				}

				infile>>t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0;k<units;k++)
						infile >> tmp;
					infile >> t;
				}
				
				for(int k=0;k<units;k++)//read units initdata
				{
					infile>>initiaValVec(nVecIndex);
					nVecIndex++;
				}
			}
			infile.close();    //close file stream
			////////////////////////////////////////////////////////////////////////////////
		}
		
		initiaValVec*=1000.0;
	}

	//step3:distribute init data to branches
	for (int i = 1; i <=nNodes; i ++)
	{
		//nodeVoltageMatrix(1, i) = initiaValVec(i);
		nodeVoltageMatrix_1[0][i-1] = initiaValVec(i);
		nodeVoltageVec(i)=initiaValVec(i);
	}

	// nodeVoltageVec.Print();

	//// debug
	//cout << "counter:" << counter << endl;
	//cout << "nodeVoltageVec:" << endl;
	//for (int k=1; k<=nNodes; k++)
	//	cout << nodeVoltageVec(k) << "  ";
	//cout << endl;

	for (int i = 1; i <= nColumns; i ++)
	{
		int idxOfI = nNodes + i;
		//branchCurrentMatrix(1, i) = initiaValVec(idxOfI);
		branchCurrentMatrix_1[0][i-1] = initiaValVec(idxOfI);
		branchCurrentVec(i) = initiaValVec(idxOfI);
	}

	// debug
	//nodeVoltageVec.Print();
	//branchCurrentVec.Print();

	int ptr = 1;//ָ��
	for (int i = 1; i <= nBranch; i ++)
	{
		Component * tempCom = branches->at(i-1);
		tempCom->initializeBranch(nodeVoltageVec,branchCurrentVec,ptr);
	}
#endif

	////////////////////////////////////////////////////////////////////
	// PWM������ƽ��ģ�Ϳ��ض���ʱ�̳�ʼ��,xuyin,20121011
	// edit by xuyin, 20130119, �޸�ʹ֮��Ӧ���������ʼ���ļ������
	if (nPWMConverter != 0) {
		int nCol_PWM = 3*nPWMConverter;
		int units_PWM = nCol_PWM%10;
		int nFiles_PWM = (nCol_PWM-units_PWM)/10;
		if(units_PWM!=0)
		{
			nFiles_PWM = nFiles_PWM + 1;
		}

		// ��������������
		int c = nCol_PWM + 1;
		int r = (len_subinterval+initDeltaT/2)/initDeltaT+1;

		// �γ�PWM������ƽ��ģ�ͳ�ʼ������
		double** PWMInitialMatrix = new double* [r];
		for (int k=0; k<r; k++) {
			PWMInitialMatrix[k] = new double [c];
		}

		// ��������
		for (int fileIndex=1;fileIndex<=nFiles_PWM;fileIndex++)
		{
			int num1 = fileIndex / 10;
			int num2 = fileIndex % 10;
			char chDir[100];
			string caseDir;
#ifdef WIN32
			sprintf(chDir,".\\result\\pscadData\\case%d_PWM_%d%d.out",caseNumber,num1,num2);   
#else
			sprintf(chDir,"..//result//pscadData//case%d_PWM_%d%d.out",caseNumber,num1,num2);   
#endif
			caseDir = chDir;
			ifstream infile(caseDir.c_str(),ios::in);

			double t; // ʱ��
			if (fileIndex < nFiles_PWM)//�������һ���ļ���ֱ�Ӷ��ߵ�һ�е�10�����ݣ���ȥʱ�������
			{
				// ���˵�startTime��ǰ������
				double tmp = 0;
				infile >> t;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// ��PWM��ʼ�����鸳ֵ
				for (int i=0; i<r; i++) {
					PWMInitialMatrix[i][0] = t;
					for (int j=1; j<=10; j++)
						infile >> PWMInitialMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			else
			{
				if (units_PWM == 0)
				{
					units_PWM = 10;
				}

				// ���˵�startTime��ǰ������
				double tmp = 0;
				infile >> t;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0; k<units_PWM; k++)
						infile >> tmp;
					infile >> t;
				}

				// ��PWM��ʼ�����鸳ֵ
				for (int i=0; i<r; i++) {
					PWMInitialMatrix[i][0] = t;
					for (int j=1; j<=units_PWM; j++)
						infile >> PWMInitialMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}

			infile.close();
			infile.clear();
		}

		// ��ʼ�����ض���ʱ��
		ptr = 1;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			tempCom->initializeSwitchingInstants(PWMInitialMatrix, r-1, ptr);
		}

		// ��ʼ��tau_s, xuyin, 20130224
		for (int k=0; k<nPWMConverter; k++)
		{
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			for (int j=0; j<3; j++)
			{
				tau_s[0][3*k+j] = tempCom->tao_s[j];
			}
		}

		// ɾ����ʱ����
		for (int k=0; k<r; k++) {
			delete[] PWMInitialMatrix[k];
		}
		delete[] PWMInitialMatrix;
	}
	//////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////
	// �첽�����ʼ��,xuyin,20121226
	if (nIndGen != 0)
	{
		int nCol_Gen = 15*nIndGen;
		int units_Gen = nCol_Gen%10;
		int nFiles_Gen = (nCol_Gen-units_Gen)/10;
		if(units_Gen!=0)
		{
			nFiles_Gen = nFiles_Gen + 1;
		}

		// �����ʼ������
		double* GenInitialMatrix[3];
		for (int k=0; k<3; k++)
		{
			GenInitialMatrix[k] = new double [nCol_Gen+1];
		}

		// ��������
		for (int fileIndex=1;fileIndex<=nFiles_Gen;fileIndex++)
		{
			int num1 = fileIndex / 10;
			int num2 = fileIndex % 10;
			char chDir[100];   
			string caseDir;
#ifdef WIN32
			sprintf(chDir,".\\result\\pscadData\\case%d_Gen_%d%d.out",caseNumber,num1,num2);   
#else
			sprintf(chDir,"..//result//pscadData//case%d_Gen_%d%d.out",caseNumber,num1,num2);   
#endif
			caseDir=chDir;
			ifstream infile(caseDir.c_str(),ios::in);

			// ���˵�startTime-2*deltaT��ǰ����
			double t; // ʱ��
			if (fileIndex < nFiles_Gen)//�������һ���ļ���ֱ�Ӷ��ߵ�һ�е�10�����ݣ���ȥʱ�������
			{
				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// �������ʼ�����鸳ֵ 
				for(int i=0;i<3;i++)//read 10 initdata
				{
					GenInitialMatrix[i][0] = t;
					for(int j=1; j<=10; j++)
						infile >> GenInitialMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			else//�������һ���ļ���ֻ��Ҫ����units����������
			{
				if (units_Gen == 0)
				{
					units_Gen = 10;
				}

				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<units_Gen;k++)
						infile >> tmp;
					infile >> t;
				}

				// �������ʼ�����鸳ֵ 
				for(int i=0;i<3;i++)//read 10 initdata
				{
					GenInitialMatrix[i][0] = t;
					for(int j=1; j<=units_Gen; j++)
						infile >> GenInitialMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			infile.close();
		}

		// ��ʼ�����
		ptr = 1; // ��һ��Ϊʱ�䣬�ʴ�1��ʼ
		for (int k=0; k<nIndGen; k++) {
			Component* tempCom = branches->at(branchNumArray_IndGen[k]);
			tempCom->initializeGen(GenInitialMatrix, ptr);
		}

		// ɾ����ʱ����
		for (int k=0; k<3; k++) {
			delete[] GenInitialMatrix[k];
		}

#ifdef WIND_VELOCITY_DATA_INPUT
		// ��������ļ�
		double** VwMatrix = new double*[rows];
		for (int k=0; k<rows; k++)
		{
			VwMatrix[k] = new double [nIndGen];
		}
		int nCol_GenVw = nIndGen;
		int units_GenVw = nCol_GenVw%10;
		int nFiles_GenVw = (nCol_GenVw-units_GenVw)/10;
		if(units_GenVw!=0)
		{
			nFiles_GenVw = nFiles_GenVw + 1;
		}
		for (int fileIndex=1;fileIndex<=nFiles_GenVw;fileIndex++)
		{
			int num1 = fileIndex / 10;
			int num2 = fileIndex % 10;
			char chDir[100];   
			string caseDir;
#ifdef WIN32
			sprintf(chDir,".\\result\\pscadData\\case%d_GenVwData_%d%d.out",caseNumber,num1,num2);   
#else
			sprintf(chDir,"..//result//pscadData//case%d_GenVwData_%d%d.out",caseNumber,num1,num2);   
#endif
			caseDir=chDir;
			ifstream infile(caseDir.c_str(),ios::in);

			// ���˵�startTime-2*deltaT��ǰ����
			double t; // ʱ��
			if (fileIndex < nFiles_Gen)//�������һ���ļ���ֱ�Ӷ��ߵ�һ�е�10�����ݣ���ȥʱ�������
			{
				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// �������ʼ�����鸳ֵ 
				for(int i=0;i<rows;i++)//read VwData
				{
					for(int j=0; j<10; j++)
						infile >> VwMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			else//�������һ���ļ���ֻ��Ҫ����units����������
			{
				if (units_GenVw == 0)
				{
					units_GenVw = 10;
				}

				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<units_Gen;k++)
						infile >> tmp;
					infile >> t;
				}

				// �������ʼ�����鸳ֵ 
				for(int i=0;i<rows;i++)//read 10 initdata
				{
					for(int j=0; j<units_GenVw; j++)
						infile >> VwMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			infile.close();
		}

		// ��ʼ�����
		ptr = 1; // ��һ��Ϊʱ�䣬�ʴ�1��ʼ
		for (int k=0; k<nIndGen; k++) {
			Component* tempCom = branches->at(branchNumArray_IndGen[k]);
			tempCom->getWindVelocityData(VwMatrix,rows,ptr);
		}

		// ɾ����ʱ����
		for (int k=0; k<rows; k++) {
			delete[] VwMatrix[k];
		}
#endif
	}	
	//////////////////////////////////////////////////////////////////////

	//����ϵͳ��ʼ��
	MeasureComponent* tempMsrCom;
	for(int i=0;i<nMsrComps;i++)
	{
		tempMsrCom = msrComps->at(i);
		tempMsrCom->initializeBranch();
	}

	//����ϵͳ��ʼ��
	///////////////////////////////////////////////////////////////////////
	// edited by xuyin, 20121109
	// ��ȡPI��������ֵ�ļ�
	if ( nPIController != 0 )
	{
		int nCol_PI = 2*nPIController;
		int units_PI = nCol_PI%10;
		int nFiles_PI = (nCol_PI-units_PI)/10;
		if(units_PI!=0)
		{
			nFiles_PI = nFiles_PI + 1;
		}

		// PI��������ʼ������
		PIControllerInitialValue = new double [nCol_PI+1];

		// ��������
		for (int fileIndex=1;fileIndex<=nFiles_PI;fileIndex++)
		{
			int num1 = fileIndex / 10;
			int num2 = fileIndex % 10;
			char chDir[100];   
			string caseDir;
#ifdef WIN32
			sprintf(chDir,".\\result\\pscadData\\case%d_PI_%d%d.out",caseNumber,num1,num2);   
#else
			sprintf(chDir,"..//result//pscadData//case%d_PI_%d%d.out",caseNumber,num1,num2);   
#endif
			caseDir=chDir;
			ifstream infile(caseDir.c_str(),ios::in);

			// ���˵�startTime��ǰ����
			double t; // ʱ��
			if (fileIndex < nFiles_PI)//�������һ���ļ���ֱ�Ӷ��ߵ�һ�е�10�����ݣ���ȥʱ�������
			{
				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// ��PI��������ʼ�����鸳ֵ 
				PIControllerInitialValue[0] = t;
				for(int j=1; j<=10; j++)
					infile >> PIControllerInitialValue[j+10*(fileIndex-1)];
			}
			else//�������һ���ļ���ֻ��Ҫ����units����������
			{
				if (units_PI == 0)
				{
					units_PI = 10;
				}

				infile >> t;
				// ���˵�startTime��ǰ������
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 )
				{
					for(int k=0;k<units_PI;k++)
						infile >> tmp;
					infile >> t;
				}

				// ��PI��������ʼ�����鸳ֵ 
				PIControllerInitialValue[0] = t;
				for(int j=1; j<=units_PI; j++)
					infile >> PIControllerInitialValue[j+10*(fileIndex-1)];
			}
			infile.close();
		}
	}
	
	// ��ʼ��
	CtrlComponent* tempCtrlCom;
	for(int i=0;i<nCtrlBranches;i++)
	{
		tempCtrlCom = ctrlBranches->at(i);
		tempCtrlCom->initializeCtrlBranch();
		tempCtrlCom->calculateCtrlEquivalentParameter();
	}
	
	// ��PI����������ʼֵ 
	ptr = 1;
	for (int i=0; i<nPIController; i++)
	{
		// ����ֵ
		tempCtrlCom = ctrlBranches->at(PIControllerBranches[i]);
		tempCtrlCom->initializePICtrl(PIControllerInitialValue, ptr);
	}
	////////////////////////////////////////////////////////////////////////

	formConductanceMatrix();//�γɽڵ㵼�ɾ���
	transferMeasurands();
	solveInitCtrlSystem();//����ʼ�����̣��õ������źŵĳ�ʼ��ֵ
	transferControlVariables();//�������ź�ֵ���䵽����ϵͳ

	if (!initializeOption)
	{
		calculateBranchNortonEquivalentCurrent(curTime);//����֧·ŵ�ٵ�Ч����
		formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
		solveNodeVoltageEquation();//���ڵ��ѹ����
		saveNodeVoltage();//����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�
		calculateBranchVoltage();//����֧·��ѹ������
		calculateBranchCurrent();//����֧·����
		saveBranchCurrent();//����֧·����
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// edit by xuyin, 20121123
	// ����У��ר��
	// �洢��ǰʱ�̵���Ԫ��״̬
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->storeInternalVariables();
	}
	// �洢����Ԫ��״̬
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->storeInternalVariables();
	}
	/////////////////////////////////////////////////////////////////////////////////////////
}

void EMTP::saveBranchCurrent()//����֧·����
{
	if (counter<finishTime/deltaT+2 && (counter%Nsave==1 || Nsave==1))
	{
		//int ptr = 1;
		int ptr = 0;
		int savecounter = (int)(counter/Nsave)+1;
		if (Nsave==1)
		{
			savecounter-=1;
		}
		for(int k=0;k<nBranch;k++)
		{
			Component * tempCom = branches->at(k);
			// tempCom->saveBranchCurrent(branchCurrentMatrix,ptr,counter);
			tempCom->saveBranchCurrent(branchCurrentMatrix_1,ptr,savecounter);
		}
	}
}

void EMTP::saveMahineWr()//������ת��
{
	if (counter<finishTime/deltaT+2 && (counter%Nsave==1 || Nsave==1))
	{
		int ptr = 0;
		int savecounter = (int)(counter/Nsave)+1;
		if (Nsave==1)
		{
			savecounter-=1;
		}
		for (int k=0; k<nIndGen; k++) {
			Component* tempCom = branches->at(branchNumArray_IndGen[k]);
			tempCom->saveMachineWr(machineWrMatrix, ptr,savecounter);
		}
	}
}

void EMTP::formConductanceMatrix()//�γɽڵ㵼����
{
	for(int i=1;i<=nNodes;i++)
	{
		for(int j=1;j<=nNodes;j++)
		{
			conductanceMatrix(i,j) = 0;
		}
	}
	for(int k=0;k<nBranch;k++)
	{
		Component * tempCom = branches->at(k);
		tempCom->formConductanceMatrix(conductanceMatrix);
	}

	lu_conductanceMatrix->SetMatrix(conductanceMatrix);
	lu_conductanceMatrix->Decompose();
}

void EMTP::calculateBranchNortonEquivalentCurrent(double time)//����֧·ŵ�ٵ�Ч����
{
	for(int k=0;k<nBranch_NEC;k++)
	{
		Component * tempCom = branches->at(branchNumArray_NEC[k]);
		tempCom->calculateNortonEquivalentCurrent(time);
	}
}

void EMTP::formNodeNortonEquivalentCurrentArray()
{
	//�γɽڵ�ŵ�ٵ�Ч��������
	for(int i=1;i<=nNodes;i++)
		nodeNortonEquivalentCurrentArray(i) = 0;

	for(int k=0;k<nBranch_NEC;k++)
	{
		Component * tempCom = branches->at(branchNumArray_NEC[k]);
		tempCom->formNodeNortonEquivalentCurrentArray(nodeNortonEquivalentCurrentArray);
	}
}

void EMTP::calculateBranchVoltage(){
	//����֧·��ѹ
	for(int k=0;k<nBranch;k++)
	{
		Component * tempCom = branches->at(k);
		tempCom->readNodeVoltage(nodeVoltageVec);
		tempCom->calculateBranchVoltage();
	}
}

void EMTP::calculateBranchCurrent(){
	//����֧·����
	for(int k=0;k<nBranch;k++)
	{
		Component * tempCom = branches->at(k);
		tempCom->calculateBranchCurrent();
	}
}

//////////////////////////////////////////////////////////////////////////
//                      ���ش�����غ�����ʼ��                          //
//////////////////////////////////////////////////////////////////////////

//���ش���������
void EMTP::switchTreatment()
{
	//ȷ�����ض���ʱ��
	bool switchOrNot = 0;//��־λ���Ƿ��п��ض���
	bool checkSwitchOrNot = 1;//��־λ���Ƿ���Ҫ���п��ع��̵��жϺʹ�����ʼֵΪ1��ÿ��ʱ����ʼʱ����Ҫ��⣩
	double ratio;//��¼���綯���Ŀ��صĲ�ֵ��
	double timeToSwitch = curTime-2*deltaT;//��¼���ض���ʱ��
	int switchID;//��¼���綯���Ŀ��صı��
	int switchCounter = 0;//ͳ����Ҫ�����Ŀ�����
	int switchMode = 0;//0Ϊ��Ȼ���ϣ�1Ϊ���ƿ���
	double maxRatio = 1;//Ϊȷ�����ض���ʱ�̶���curTime-deltaTʱ��֮ǰ�����õ���ֵ

	while (checkSwitchOrNot)
	{
		//ÿ�ο������̴���Ŀ�ʼ����ʼ������
		switchID = -1;checkSwitchOrNot=0;switchMode=0;

		//ȷ����Ҫ�����Ŀ��أ�ȷ���׸��������صĲ�ֵ���Լ����ض���ʱ��
		checkSwitch(ratio,switchID,maxRatio,switchCounter);
		
		//����п��ض�����������Ӧ����
		if (switchID!=-1)
		{
			timeToSwitch = timeToSwitch + ratio*deltaT;//����п��ض���������timeToSwitch
			maxRatio = (curTime-deltaT-timeToSwitch)/deltaT;

			switchOrNot = 1;
			checkSwitchOrNot = 1;//��Ȼ��Ҫ��һ���жϿ��ض���
			doSwitch(ratio,switchCounter, switchMode);//������Ҫ�����Ŀ��أ������γɵ��ɾ���
	
			linearInt(ratio);//(1)���Բ�ֵ�����timeToSwitchʱ�̵Ľ�
			if (switchMode!=0)//���ƿ��Ϲ��̣����ֶ�����ŵ�ٵ�Ч����ֵ���䣬��⿪�ض���t+ʱ��״̬
			{
				networkSolution();//���t+ʱ��״̬
				//cout<<"�ڲ�ѭ����ʼ"<<endl;
				Component* tempCom;
				int needReCal = 1;
				int nCount=0;
				while(needReCal)
				{
					needReCal = 0;
					for (int k=0;k<nSwitch;k++)
					{
						tempCom = branches->at(diodeNumArray[k]);
						if (tempCom->checkSwitch(curTime)==1)
						{
							nCount++;
							if (nCount>100)
							{
								int pauseit=1;
							}
							tempCom->switchIt();
							tempCom->modifyConductanceMatrix(conductanceMatrix);
							lu_conductanceMatrix->SetMatrix(conductanceMatrix);
							lu_conductanceMatrix->Decompose();
							networkSolution();
							needReCal = 1;
						}
					}
				}
				//cout<<"�ڲ�ѭ������"<<endl;
			}
			calculateBranchNortonEquivalentCurrent_forBasicComp(timeToSwitch+deltaT);//����֧·ŵ�ٵ�Ч����
			networkSolution();//(2)���timeToSwitch+deltaTʱ�̵Ľ�
		}		
	}

	//����Chatter�ʹ洢���
	if (switchOrNot)
	{
		//removeChatter(timeToSwitch);

		linearInt((curTime-deltaT-timeToSwitch)/deltaT);

		transferMeasurands();//������ϵͳ����ֵ���䵽����ϵͳ

		solveCtrlSystem();

		transferControlVariablesForSwitch();//�������ź�ֵ����������ϵͳ

		saveNewResult(counter-1);
	}
	else
	{
		/************************************************************************/
		/*    ����Ƿ���Inverter״̬�����仯�����б仯�������γɵ��ɾ���        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/
	}
} 

//�ж���Ҫ�����Ŀ��أ��Լ����ȶ����Ŀ��ض�Ӧ��ratio
void EMTP::checkSwitch(double& ratio,int& switchID,double maxRatio,int& switchCounter)
{
	Component* tempCom;
	double temp;
	int tempNum;
	ratio = 1;switchCounter = 0;switchID = -1;
	for (int k=0;k<nSwitch;k++)
	{
		tempNum = diodeNumArray[k];
		tempCom = branches->at(tempNum);
		if (tempCom->checkSwitch(curTime)==1 && (temp = tempCom->getSwitchRatio())<maxRatio)
		{
			switchNumArray[switchCounter] = tempNum;
			switchRatioArray[switchCounter] = temp;
			switchModeArray[switchCounter] = tempCom->getSwitchMode();
			switchCounter++;
			if (temp<ratio)
			{
				ratio = temp;
				switchID = tempNum;
			}
		}
	}
}

//������Ҫ�����Ŀ��أ������γɵ��ɾ���
void EMTP::doSwitch(double ratio,int switchCounter, int& switchMode)
{
	
	//�任��Ҫ�����Ŀ���
	Component* tempCom;
	for (int k=0;k<switchCounter;k++)
	{
		if ((switchRatioArray[k]-ratio)<0.0001)
		{
			tempCom = branches->at(switchNumArray[k]);
			tempCom->switchIt();
			tempCom->modifyConductanceMatrix(conductanceMatrix);
			switchMode += switchModeArray[k];
			//if (switchMode!=0)
			//{
			//	break;
			//}
		}
	}

	//���¶Ե��ɾ������lu�ֽ�
	lu_conductanceMatrix->SetMatrix(conductanceMatrix);
	lu_conductanceMatrix->Decompose();
}

void EMTP::removeChatter(double timeToSwitch)
{
	if(timeToSwitch<curTime-3.0/2.0*deltaT)//������Ŀ��ض�������ǰ���������
	{
		//���ð벽����ֵ�㷨���õ�curTime-deltaTʱ�̵Ľ�
		linearInt(0.5);
		calculateBranchNortonEquivalentCurrent_forBasicComp(timeToSwitch+3.0/2.0*deltaT);//����֧·ŵ�ٵ�Ч����
		networkSolution();
		linearInt((curTime-timeToSwitch-deltaT*3.0/2.0)/deltaT);//
		//�������������curTime-deltaTʱ�̵Ĵ洢���
		saveNewResult(counter-1);
		updateResult(3);

		/************************************************************************/
		/*    ����Ƿ���Inverter״̬�����仯�����б仯�������γɵ��ɾ���        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/
	}
	else//������Ŀ��ض������ں���������
	{
		//��1�����Ȳ�ֵ�ص���ʱ����curTime-deltaT����������������еĴ洢���
		linearInt((curTime-deltaT-timeToSwitch)/deltaT);
		saveNewResult(counter-1);

		/************************************************************************/
		/*    ����Ƿ���Inverter״̬�����仯�����б仯�������γɵ��ɾ���        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/

		//��2����ǰ����һ�Σ��õ�curTimeʱ�̵Ľ�
		updateResult(3);
		calculateBranchNortonEquivalentCurrent(curTime);//����֧·ŵ�ٵ�Ч����
		networkSolution();
		//����Chatterǰ��_1������ֵ����_2�����н��б���
		updateResult(1);
		//��3�����ð벽����ֵ�㷨��ȥ��curTimeʱ�̵Ľ��е�Chatter
		linearInt(0.5);
		calculateBranchNortonEquivalentCurrent_forBasicComp(curTime+deltaT/2.0);//����֧·ŵ�ٵ�Ч����
		networkSolution();
		//cout<<"��2�ΰ벽����ֵ��"<<endl;
		linearInt(0.5);
		//����curTimeʱ�̵Ľ�������������
		saveNewResult(counter);
		//����Chatter����_2������ֵ��ԭ_1������ֵ
		updateResult(2);
		//��4��ʱ����ǰ�ƽ�һ��ʱ����������һ�ֵĿ����ж�
		advanceTime();
		switchTreatment();
	}
}

//һ�����Բ�ֵ
void EMTP::linearInt(double ratio)
{
	Component* tempCom;

	//�ڵ��ѹ�Ĳ�ֵ
	for (int k=1;k<=nNodes;k++)
	{
		nodeVoltageVec(k) = nodeVoltageVec_1(k) + ratio*(nodeVoltageVec(k)-nodeVoltageVec_1(k));
	}

	//��֧·��ѹ�����Ĳ�ֵ
	for(int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		tempCom->interpolate(ratio);
	}
}

//����֧·ŵ�ٵ�Ч����
void EMTP::calculateBranchNortonEquivalentCurrent_forBasicComp(double time)
{
	Component * tempCom;
	for(int k=0;k<nBranch_NEC;k++)
	{
		tempCom = branches->at(branchNumArray_NEC[k]);
		if (!tempCom->isUserDef)//�������Զ���Ԫ��
		{
			tempCom->calculateNortonEquivalentCurrent(time);
		}
	}
}

//һ���������
void EMTP::networkSolution()
{
	formNodeNortonEquivalentCurrentArray();//�γɽڵ�ŵ�ٵ�Ч��������
	solveNodeVoltageEquation();//���ڵ��ѹ����
	calculateBranchVoltage();//����֧·��ѹ
	calculateBranchCurrent();//����֧·����
}

//���¿��ش�����������ڴ洢����ı���
void EMTP::updateResult(int updateTypeNum)
{
	Component * tempCom;//ָ��Ԫ����ָ��
	double temp;//��ʱ�������������ݽ���

	switch (updateTypeNum)
	{
	case 1://��_1��������ֵ����_2������
		for(int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			tempCom->updateResult(1);
		}
		for (int k=1;k<=nNodes;k++)
		{
			nodeVoltageVec_2(k) = nodeVoltageVec_1(k);
		}
		break;
	case 2://��_2��������ֵ����_1������
		for(int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			tempCom->updateResult(2);
		}
		for (int k=1;k<=nNodes;k++)
		{
			nodeVoltageVec_1(k) = nodeVoltageVec_2(k);
		}
		break;
	case 3:
		for(int k=0;k<nBranch_userDef;k++)
		{
			tempCom = branches->at(branchNumArray_userDef[k]);
			tempCom->updateResult(3);
		}
		break;
	default:
		cerr<<"��������ȷ�ĸ������ͱ�ţ�"<<endl;
		exit(1);
	}
}

//����ֵ�õ��Ľ��������������
void EMTP::saveNewResult(int counter)
{
	if (counter<finishTime/deltaT+2 && (counter%Nsave==1 || Nsave==1))
	{
		Component* tempCom;//��ʱ������ɨ��֧·ʱ��
		int savecounter = (int)(counter/Nsave)+1;
		if (Nsave==1)
		{
			savecounter-=1;
		}

		//���½ڵ��ѹ������
		for(int i=1;i<=nNodes;i++)
		{
			//nodeVoltageMatrix(counter,i)=nodeVoltageVec(i);
			nodeVoltageMatrix_1[savecounter-1][i-1] = nodeVoltageVec(i);
		}
	
		//����֧·����������
		//int ptr = 1;
		int ptr = 0;
		for(int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			//tempCom->saveBranchCurrent(branchCurrentMatrix,ptr,counter);
			tempCom->saveBranchCurrent(branchCurrentMatrix_1,ptr,savecounter);
		}
	}
	
}

//ʱ��ǰ��һ������
void EMTP::advanceTime()
{
	curTime += deltaT;
	counter++;
	if (counter<finishTime/deltaT+2)
	{
		timeArray(counter) = curTime;
	}
}

//����Ƿ���Inverter״̬�����仯�����б仯�������γɵ��ɾ���
void EMTP::inverterTreatment()
{
	//�ж��Ƿ���Inverter״̬�仯
	bool changeOrNot = 0;
	int tempNum;
	Component* tempCom;
	for (int k=0;k<nInverter;k++)
	{
		tempNum = inverterNumArray[k];
		tempCom = branches->at(tempNum);
		if (tempCom->checkSwitch(counter,conductanceMatrix))
		{
			changeOrNot = 1;
		}
	}

	for (int k=0;k<nTimeSwitch;k++)
	{
		tempNum = timeSwitchNumArray[k];
		tempCom = branches->at(tempNum);
		if (tempCom->checkSwitch(curTime))
		{
			changeOrNot = 1;
			formConductanceMatrix();
		}
	}

	//�����Inverter״̬�����仯�������γɵ��ɾ���
	if (changeOrNot)
	{
		lu_conductanceMatrix->SetMatrix(conductanceMatrix);
		lu_conductanceMatrix->Decompose();
	}
}

//////////////////////////////////////////////////////////////////////////
//              ���ش�����غ�����ֹ��                                  //
//////////////////////////////////////////////////////////////////////////

void EMTP::getCaseDifinitionMatrix(int caseNumber)
{
	caseDfnMatrix.ResizeTo(2000,10);
	char ch[100];     
	string dir1;
#ifdef WIN32
	sprintf(ch,".\\caseDfn\\case%d.txt",caseNumber);    
#else
	sprintf(ch,"..//caseDfn//case%d.txt",caseNumber);    
#endif
	dir1=ch;
	ifstream infile(dir1.c_str(),ios::in);
	nBranch=0;
	while(!infile.eof())
	{
		infile>>caseDfnMatrix(nBranch,0);
		if(infile.fail()) //ֻ���� >> ֮������жϣ�������Ч
		{				
			caseExistFlag=0;
			cout<<"case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//ֱ�ӷ���
			break;
		}

		for(int j=1;j<10;j++)
		{
			infile>>caseDfnMatrix(nBranch,j);

		}
		nBranch++;
	}
	caseDfnMatrix.ResizeTo(nBranch,10);
	infile.close();      

}

void EMTP::	specifySystem()
{
	caseExistFlag=1;
#ifdef  WIN32
	ifstream infile("configure.txt",ios::in);
#else
	ifstream infile("../configure.txt",ios::in);
#endif
	infile>>startTime>>finishTime>>deltaT>>saveDeltaT>>initDeltaT>>len_subinterval>>caseNumber>>initializeOption;
	infile.close(); 

	//========================================
	// xuyin, 20130224
	Np = int((finishTime-startTime+deltaT/2)/len_subinterval);
	Ns = int((len_subinterval+deltaT/2)/deltaT);
	Nsave=int((saveDeltaT+deltaT/2)/deltaT);
	//========================================

	getCaseDifinitionMatrix(caseNumber);
	nNodes=(int)caseDfnMatrix.GetSub(0,caseDfnMatrix.GetRowUpb(),1,2).Max();//�ڵ�����

	rows=1;//ʱ�䣬��ѹ������ȫ������ĳ���
	double t = startTime;
	while(t<finishTime-saveDeltaT/2)
	{
		t=t+saveDeltaT;
		rows++;
	}
	nColumns=nBranch;

	int type=0;
	int resistanceCount=0;
	int inductanceCount=0;
	int capacitanceCount=0;
	int acVoltageSourceCount=0;
	int TimedAcVoltageSourceCount=0;
	int transformerCount=0;
	int synchronousMachineCount=0;
	int inductionGeneratorCount=0;
	int WoundInducGeneratorCount=0;
	int WindInducGeneratorCount=0;
	int diodeCount=0;
	int timeSwitchCount=0;
	int mwSynGenerator12Count=0;
	int Motor15Count=0;
	int diodebridgeCount=0;
	int IGBTCount=0;
	int RectifierBridgeCount=0;
	int ThreeRectifierBridgeCount=0;
	int ThreePhaseRectifierBridgeCount=0;
	int InverterCount=0;
	int diodeWithSCCount=0;
	int Asych5Count=0;
	int ImpedanceCount=0;
	int newIGBTCount=0;
	int ControlledvoltageSourceCount=0;
	int PWMConverterCount = 0;
	int ThreePhaseTransformerCount = 0;

	for(int i=0;i<nBranch;i++)
	{   
		type=(int)caseDfnMatrix(i,0);
		Component * com = NULL;
		int nodeNumberArray[4] = {(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3),(int)caseDfnMatrix(i,4)};
		double apparentPowerRating = caseDfnMatrix(i,5);
		double voltageRating[2] = {caseDfnMatrix(i,6),caseDfnMatrix(i,7)};

		//for diodebridge,thyristorbridge
		int in_ACNode[3]={(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3)};
		int in_DCNode[2]={(int)caseDfnMatrix(i,4),(int)caseDfnMatrix(i,5)};
		int in_state[6]={0,0,0,0,0,0};
		//for thyristorbridge
		double in_alpha=1.0;//09.03.13
		//for threephasetransformer2
		int TPTnodeNumberArray[12] = {(int)caseDfnMatrix(i,1),0,(int)caseDfnMatrix(i,2),0,caseDfnMatrix(i,3),0,caseDfnMatrix(i,4),0,caseDfnMatrix(i,5),0,caseDfnMatrix(i,6),0};
		double TPTvoltageRating[2] = {caseDfnMatrix(i,8),caseDfnMatrix(i,9)};
		//for PWM
		/*int in_DCNode_2[2] = {(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2)};
		int in_ACNode_2[3] = {(int)caseDfnMatrix(i,3),(int)caseDfnMatrix(i,4),(int)caseDfnMatrix(i,5)};*/
		//for SPWM
		//int SPWMnodeNumberArray[5] = {(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5)};

		switch (type)
		{
		case 1: //resistance(1,'R',7,10,pscad1_1_3ph_R);
			resistanceCount=resistanceCount+1;   
			com =new Resistance(resistanceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3));
			branches->push_back(com);
			break;
		case 2:
			inductanceCount=inductanceCount+1;
			com =new Inductance(inductanceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3));
			branches->push_back(com);
			break;
		case 3:
			capacitanceCount=capacitanceCount+1;
			com =new Capacitance(capacitanceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3));
			branches->push_back(com);
			break;
		case 4:
			acVoltageSourceCount=acVoltageSourceCount+1;
			com =new AcVoltageSource(acVoltageSourceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6));
			branches->push_back(com);
			break;
		case 5:
			transformerCount=transformerCount+1;
			com =new SinglePhaseTransformer(transformerCount,nodeNumberArray,apparentPowerRating,voltageRating);
			branches->push_back(com);
			nColumns=nColumns+1;
			break;
		case 6:
			synchronousMachineCount=synchronousMachineCount+1;
			com =new SynchGenerator((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3));
			branches->push_back(com);
			nColumns=nColumns+2;
			break;
		case 7:
			diodeCount += 1;
			com =new Diode(diodeCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5));
			branches->push_back(com);
			break;
		case 8:
			timeSwitchCount += 1;
			com = new TimeSwitch(timeSwitchCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6));
			branches->push_back(com);
			break;
			//case 9:
			//	diodebridgeCount += 1;
			//	//Diodebridge(int id,int *ACNode,int *DCNode,int *state,double onValue,double offValue);

			//	com = new Diodebridge(diodebridgeCount,in_ACNode,in_DCNode,in_state,1.0e-2,1.0e6);
			//	branches->push_back(com);
			//	nColumns=nColumns+5;
			//	break;
		case 10:
			inductionGeneratorCount+=1;
			com=new InducGenerator((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3));
			branches->push_back(com);
			nColumns=nColumns+2;
			break;
		case 11:
			InverterCount += 1;
			com = new Inverter(InverterCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3),(int)caseDfnMatrix(i,4),(int)caseDfnMatrix(i,5),caseDfnMatrix(i,6),caseDfnMatrix(i,7));
			branches->push_back(com);
			nColumns=nColumns+19;
			break;
		case 12:
			mwSynGenerator12Count = mwSynGenerator12Count + 1;
			com =new MWSynGenerator12((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2));
			branches->push_back(com);
			nColumns=nColumns+11;
			break;
		case 13:
			IGBTCount += 1;
			com = new IGBT(IGBTCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6));
			branches->push_back(com);
			break;
			//case 14:	//���Ŀ��Ʒ���
			//	RectifierBridgeCount += 1;
			//	com = new RectifierBridge(RectifierBridgeCount,in_ACNode,in_DCNode,in_state,1.0e-2,1.0e6);
			//	branches->push_back(com);
			//	nColumns=nColumns+5;
			//	break;
		case 15:
			Motor15Count += 1;
			com = new Motor15((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2));
			branches->push_back(com);
			nColumns = nColumns + 14;
			break;
			//case 16://���ݽڵ��ѹ�жϵ�����
			//	ThreeRectifierBridgeCount+= 1;
			//	com = new ThreeRectifierBridge(ThreeRectifierBridgeCount,in_ACNode,in_DCNode,in_state,1.0e-2,1.0e6);
			//	branches->push_back(com);
			//	nColumns=nColumns+5;
			//	break;
		case 17:
			ThreePhaseRectifierBridgeCount += 1;
			com = new ThreePhaseRectifierBridge(ThreePhaseRectifierBridgeCount ,in_ACNode,in_DCNode,in_state,1.0e-2,1.0e6);
			branches->push_back(com);
			nColumns=nColumns+5;
			break;
		case 18:
			WoundInducGeneratorCount+=1;
			com=new WoundInducGenerator(WoundInducGeneratorCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),(int)caseDfnMatrix(i,3),caseDfnMatrix(i,4));
			branches->push_back(com);
			nColumns=nColumns+5;
			break;
		case 19:
			WindInducGeneratorCount+=1;
			com=new WindInducGenerator((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2));
			branches->push_back(com);
			nColumns=nColumns+5;
			break;
		case 20:
			TimedAcVoltageSourceCount=TimedAcVoltageSourceCount+1;
			com =new TimedAcVoltageSource(TimedAcVoltageSourceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6),caseDfnMatrix(i,7),caseDfnMatrix(i,8),caseDfnMatrix(i,9));
			branches->push_back(com);
			break;
		case 21:
			ImpedanceCount+=1;
			com=new Impedance(ImpedanceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4));
			branches->push_back(com);
			break;
		case 22:
			newIGBTCount += 1;
			com =new newIGBT(newIGBTCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5));
			branches->push_back(com);
			break;
		case 23:
			ControlledvoltageSourceCount+=1;
			com = new ControlledVoltageSource(ControlledvoltageSourceCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),(int)caseDfnMatrix(i,4));
			branches->push_back(com);
			break;
		case 24:
			PWMConverterCount += 1;
			com = new PWMConverter((int)caseDfnMatrix(i,1),caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6),caseDfnMatrix(i,7));
			branches->push_back(com);
			nColumns = nColumns + 4;
			break;
		case 25:
			ThreePhaseTransformerCount += 1;
			com = new ThreePhaseTransformer(ThreePhaseTransformerCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6),caseDfnMatrix(i,7),caseDfnMatrix(i,8),caseDfnMatrix(i,9));
			branches->push_back(com);
			nColumns = nColumns + 6;
			break;
		case 77:
			diodeWithSCCount += 1;
			com =new DiodeWithSC(diodeCount,(int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2),caseDfnMatrix(i,3),caseDfnMatrix(i,4),caseDfnMatrix(i,5),caseDfnMatrix(i,6),caseDfnMatrix(i,7));
			branches->push_back(com);
			break;
		case 555:
			Asych5Count += 1;
			com = new Asych5((int)caseDfnMatrix(i,1),(int)caseDfnMatrix(i,2));
			branches->push_back(com);
			nColumns = nColumns + 4;
			break;

		default:
			cout<<"Warning,There is no component with type"<<type<<endl;
		}
	}

	nSolveG=0;
	timeArray.ResizeTo(1,rows);//ʱ������
	nodeVoltageVec.ResizeTo(1,nNodes);//�ڵ��ѹ����
	nodeVoltageVec_1.ResizeTo(1,nNodes);
	nodeVoltageVec_2.ResizeTo(1,nNodes);
	nodeVoltageMatrix.ResizeTo(1,rows,1,nNodes);//�ڵ��ѹ����

	nodeVoltageMatrix_1 = new double*[rows];
	for (int k=0;k<rows;k++)
	{
		nodeVoltageMatrix_1[k] = new double[nNodes];
	}

	for(int i=0;i<rows;i++)
	{
		for (int j=0;j<nNodes;j++)
		{
			nodeVoltageMatrix_1[i][j] = 0;
		}
	}

	conductanceMatrix.ResizeTo(1,nNodes,1,nNodes);//�ڵ㵼�ɾ���
	resistanceMatrix.ResizeTo(1,nNodes,1,nNodes);;//�ڵ��迹����
	nodeNortonEquivalentCurrentArray.ResizeTo(1,nNodes);//�ڵ�ŵ�ٵ�ֵ��������
	lu_conductanceMatrix= new TDecompLU(conductanceMatrix);

	branchCurrentMatrix.ResizeTo(1,rows,1,this->nColumns);

	branchCurrentMatrix_1 = new double*[rows];
	for (int k=0;k<rows;k++)
	{
		branchCurrentMatrix_1[k] = new double[this->nColumns];
	}

	for(int i=0;i<rows;i++)
	{
		for (int j=0;j<this->nColumns;j++)
		{
			branchCurrentMatrix_1[i][j] = 0;
		}
	}

	/************************************************************************/
	/*                     ���ش�����ر���������Ĵ���                     */
	/************************************************************************/
	Component* tempCom;//��ʱ������ɨ��֧·ʱ��

	//���������
	nSwitch = diodeCount+newIGBTCount;
	if (nSwitch != 0)
	{
		diodeNumArray = new int[nSwitch];
		switchNumArray = new int[nSwitch];//��Ҫ�����Ŀ��صı��
		switchRatioArray = new double[nSwitch];//������Ҫ�����Ŀ��ض�Ӧ��ratio
		switchModeArray = new int[nSwitch];//�������ض������õ�mode

		int k1=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if (tempCom->type == 7||tempCom->type==22)
			{
				diodeNumArray[k1] = k;
				k1++;
			}
		}
	}
	else
	{
		diodeNumArray = NULL;
		switchNumArray = NULL;
		switchRatioArray = NULL;
	}
	
	//Inverter���
	nInverter = InverterCount;
	if (nInverter != 0)
	{
		inverterNumArray = new int[nInverter];

		int k2=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->type == 11)
			{
				inverterNumArray[k2] = k;
				k2++;
			}
		}
	}
	else
	{
		inverterNumArray = NULL;
	}

	//TimeSwitch���
	nTimeSwitch = timeSwitchCount;
	if (nTimeSwitch != 0)
	{
		timeSwitchNumArray = new int[nTimeSwitch];

		int k2=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->type == 8)
			{
				timeSwitchNumArray[k2] = k;
				k2++;
			}
		}
	}
	else
	{
		timeSwitchNumArray = NULL;
	}

	/************************************************************************/
	/*                    ŵ�ٵ�ֵ���������������Ĵ���                    */
	/************************************************************************/
	nBranch_NEC=0;
	for (int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		if(tempCom->need_NEC==1)
		{
			nBranch_NEC++;
		}
	}
	if (nBranch_NEC != 0)
	{
		branchNumArray_NEC = new int[nBranch_NEC];
		int k3=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->need_NEC==1)
			{
				branchNumArray_NEC[k3] = k;
				k3++;
			}
		}
	}
	else
	{
		branchNumArray_NEC = NULL;
	}

	/************************************************************************/
	/*                        �Զ���Ԫ������Ĵ���                          */
	/************************************************************************/
	nBranch_userDef=0;
	for (int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		if(tempCom->isUserDef==1)
		{
			nBranch_userDef++;
		}
	}	
	if (nBranch_userDef != 0)
	{
		branchNumArray_userDef = new int[nBranch_userDef];
		int k4=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->isUserDef==1)
			{
				branchNumArray_userDef[k4] = k;
				k4++;
			}
		}
	}
	else
	{
		branchNumArray_userDef = NULL;
	}

	/*************************************************************/
	/*                 �ܿ�Ԫ��������鴴��                      */
	/*************************************************************/

	nControlledBranches = ControlledvoltageSourceCount+newIGBTCount;
	if (nControlledBranches != 0)
	{
		controlledBranches = new int[nControlledBranches];
		int k3=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->type==22||tempCom->type==23)
			{
				controlledBranches[k3] = k;
				k3++;
			}
		}
	}
	else
	{
		controlledBranches = NULL;
	}

	/*************************************************************/
	/*          PWM������ƽ��ģ��������鴴��                    */
	/*************************************************************/
	nPWMConverter = PWMConverterCount;
	if ( nPWMConverter != 0 )
	{
		branchNumArray_PWMConverter = new int[nPWMConverter];
		ctrlNumArray_PWMConverter = new int[4*nPWMConverter];
		int k3=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->type == 24)
			{
				branchNumArray_PWMConverter[k3] = k;
				tempCom->ctrlNodeNumber_PWM(ctrlNumArray_PWMConverter, k3);
				k3++;
			}
		}
	}
	else
	{
		branchNumArray_PWMConverter = NULL;
		ctrlNumArray_PWMConverter = NULL;
	}

	// xuyin, 20130224
	if (nPWMConverter == 0)
	{
		tau_s = NULL;
	} 
	else
	{
		tau_s = new double* [Np];
		for (int k=0; k<Np; k++)
		{
			tau_s[k] = new double [3*nPWMConverter];
		}
	}

	/*************************************************************/
	/*                   �첽���������鴴��                    */
	/*************************************************************/
	nIndGen = WoundInducGeneratorCount;
	if ( nIndGen != 0 )
	{
		branchNumArray_IndGen = new int[nIndGen];
		int k3=0;
		for (int k=0;k<nBranch;k++)
		{
			tempCom = branches->at(k);
			if(tempCom->type == 18)
			{
				branchNumArray_IndGen[k3] = k;
				k3++;
			}
		}
	}
	else
	{
		branchNumArray_IndGen = NULL;
	}

	machineWrMatrix = new double*[rows];
	for (int k=0;k<rows;k++)
	{
		machineWrMatrix[k] = new double[nIndGen];
	}

	for(int i=0;i<rows;i++)
	{
		for (int j=0;j<nIndGen;j++)
		{
			machineWrMatrix[i][j] = 0;
		}
	}
	/************************************************************/
	/*                  ���÷��沽��                            */
	/************************************************************/
	for (int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		tempCom->setDetalT(deltaT);
		tempCom->calculateNortonEquivalentResistance();
	}
}


void EMTP::solveNodeVoltageEquation()//���ڵ��ѹ����
{
	bool ok;
	nodeVoltageVec_1 = nodeVoltageVec;
	nodeVoltageVec=lu_conductanceMatrix->Solve(nodeNortonEquivalentCurrentArray,ok);
}

void EMTP::saveNodeVoltage()//����ڵ��ѹ
{
	if (counter<finishTime/deltaT+2 && (counter%Nsave==1 || Nsave==1))
	{
		//for(int i=1;i<=nNodes;i++)
		//{
		//	nodeVoltageMatrix(counter,i)=nodeVoltageVec(i);
		//}
		int savecounter = (int)(counter/Nsave)+1;
		if (Nsave==1)
		{
			savecounter-=1;
		}
		for(int i=0;i<nNodes;i++)
		{
			nodeVoltageMatrix_1[savecounter-1][i]=nodeVoltageVec(i+1);
		}
	}
}
void EMTP::saveResult(void)
{
	int row=0;            
	row=timeArray.GetNrows();
	//resultMatrix.ResizeTo(1,row,1,nNodes+nColumns);
	char ch1[100];
	char ch2[100];
	string dir1;     
	string dir2;     
#ifdef  WIN32
	sprintf(ch1,".\\result\\result_latest\\result%d.dat",caseNumber);    
	sprintf(ch2,".\\result\\result_latest\\result%dName.dat",caseNumber);  
#else
	sprintf(ch1,"..//result//result_latest//result%d.dat",caseNumber);    
	sprintf(ch2,".//result//result_latest//result%dName.dat",caseNumber);  
#endif
	dir1=ch1;
	dir2=ch2;

	ofstream outfile1(dir1.c_str(),ios::out);
	for(int i=1;i<=row;i++)                                    
	{                        
		outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<timeArray(i)<<" ";
		for(int j=1;j<=nNodes;j++)                                               
		{                                                                   
			//resultMatrix(i,j)=nodeVoltageMatrix(i,j);
			//resultMatrix(i,j) = nodeVoltageMatrix_1[i-1][j-1];
			//outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<nodeVoltageMatrix(i,j)<<" ";
			outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<nodeVoltageMatrix_1[i-1][j-1]<<" ";
		}  
		for(int k=1;k<=nColumns;k++)                                               
		{       
			//resultMatrix(i,nNodes+k)=branchCurrentMatrix(i,k);
			//resultMatrix(i,nNodes+k)=branchCurrentMatrix_1[i-1][k-1];
			//outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<branchCurrentMatrix(i,k)<<" ";
			outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<branchCurrentMatrix_1[i-1][k-1]<<" ";
		}  
		for(int k=1;k<=nIndGen;k++)                                               
		{       
			outfile1<<setiosflags(ios::left)<<setiosflags(ios::showpoint)<<setiosflags(ios::scientific)<<setprecision(15)<<setw(22)<<setfill(' ')<<machineWrMatrix[i-1][k-1]<<" ";
		}
		outfile1<<endl;
	}     

	outfile1.close();                               

	ofstream outfile2(dir2.c_str(),ios::out);
	//ofstream outfile2("varName.dat",ios::out);

	for(int i=1;i<=nNodes;i++)                                               
	{                                                                   
		outfile2<<"Voltage_of_node"<<i<<endl;
	}  
	for(int j=1;j<=nColumns;j++)                                               
	{                                                                   
		outfile2<<"Current_of_branch"<<j<<endl;
	} 
	for (int j=1;j<=nIndGen;j++)
	{
		outfile2<<"Wr_of_InducMachine"<<j<<endl;
	}
	//�ڶ��ַ��������֧·�����֣���Ӧ֧·����		
	/*
	for(int j=1;j<=nBranch;j++)                                               
	{                
	Component * tempCom = branches->at(j-1);

	outfile2<<"Current_of_branch_"<<tempCom->name<<endl;
	}	
	*/  
	outfile2.close();                               
}


void EMTP::getPSCADResult(void)//return result0 = PSCAD result 
{

	int nCol=0;//columns  of result except
	nCol=nNodes+nBranch;

	int units=0;//columns of the last fie
	units=nCol%10;

	int nFiles=0;//count of  *.out files
	nFiles=(nCol-units)/10;		//	step1:cal number of *.out files 
	if(units!=0)//�����λ��Ϊ0
	{
		nFiles=nFiles+1;
	}

	//	step2:scan all *.out files and read init data
	int row=0;            
	row=timeArray.GetNrows();//row of result matrix

	int IndxCol=1;//column index of resultP,fill resultP from column 1

	for (int fileIndex=1;fileIndex<=nFiles;fileIndex++)
	{	
		int num2;
		num2=fileIndex%10;

		int num1;
		num1=(fileIndex-num2)/10;

		char chDir[100];   
#ifdef WIN32
		sprintf(chDir,".\\result\\pscadData\\case%d_%d%d.out",caseNumber,num1,num2);   
#else		
		sprintf(chDir,"..//result//pscadData//case%d_%d%d.out",caseNumber,num1,num2);   
#endif
		string caseDir;
		caseDir=chDir;
		//	cout<<"chDir = "<<chDir<<endl;
		//	cout<<"caseDir = "<<caseDir<<endl;
		ifstream infile(caseDir.c_str(),ios::in);//open cur file stream

		if (fileIndex<nFiles)//not the last file ,read 10 data except var t
		{
			for(int i=1;i<=row;i++)
			{
				double t;
				infile>>t;
				for(int k=0;k<10;k++)//read 10 initdata
				{
					infile>>resultMatrix0(i,IndxCol+k);
				}

			}
			IndxCol=IndxCol+10;
		}
		else// this is the last file
		{
			for(int i=1;i<=row;i++)
			{
				double t;
				infile>>t;
				for(int k=0;k<units;k++)//read units  initdata
				{
					infile>>resultMatrix0(i,IndxCol+k);
				}
			}
			IndxCol=IndxCol+units;

		}
		infile.close();    //close file stream
	}
	resultMatrix0*=1000.0;
}

void EMTP::getMATLABResult(int option)//return result0 = MATLAB result 
{
	//option=0,init0 result
	//option!=0,initP result
	char ch[100];
	string dir1;   

	if(option==0)
	{
#ifdef WIN32
		sprintf(ch,".\\result\\MatlabData_init0\\Mresult%d.dat",caseNumber);   
#else
		sprintf(ch,"..//result//MatlabData_init0//Mresult%d.dat",caseNumber);   
#endif
	}
	else
	{
#ifdef WIN32
		sprintf(ch,".\\result\\MatlabData_initP\\Mresult%d.dat",caseNumber); 
#else
		sprintf(ch,"..//result//MatlabData_initP//Mresult%d.dat",caseNumber); 
#endif
	}

	dir1=ch;
	ifstream infile(dir1.c_str(),ios::in);

	int row=0;            
	row=timeArray.GetNrows();//row of result matrix
	for(int i=1;i<=row;i++)                                    
	{                                                                     
		double t;
		infile>>t;

		for(int j=1;j<=nNodes+this->nColumns;j++)                                               
		{                                                                   
			infile>>resultMatrix0(i,j);

		}  
	}     
	infile.close();    

}
void EMTP::getCPPResult(int option)//return result0 = CPP(Windows) result
{
	//option=0,init0 result
	//option!=0,initP result
	char ch[100];
	string dir1;   

	if(option==0)
	{
#ifdef WIN32
		sprintf(ch,".\\result\\result_CPP_init0\\result%d.dat",caseNumber);   
#else
		sprintf(ch,"..//result//result_CPP_init0//result%d.dat",caseNumber);   
#endif
	}
	else
	{
#ifdef WIN32
		sprintf(ch,".\\result\\result_CPP_initP\\result%d.dat",caseNumber); 
#else
		sprintf(ch,"..//result//result_CPP_initP//result%d.dat",caseNumber); 
#endif
	}

	dir1=ch;
	ifstream infile(dir1.c_str(),ios::in);

	int row=0;            
	row=timeArray.GetNrows();//row of result matrix
	for(int i=1;i<=row;i++)                                    
	{                                                                     
		double t;
		infile>>t;

		for(int j=1;j<=nNodes+nBranch;j++)                                               
		{                                                                   
			infile>>resultMatrix0(i,j);

		}  
	}     
	infile.close();    

}


void EMTP::checkResult(void)
{
	//��ʼ���ο��������
	int row=0;            
	row=timeArray.GetNrows();
	resultMatrix0.ResizeTo(1,row,1,nNodes+this->nColumns);
	resultMatrix0.Zero();
	getPSCADResult();//----��PSCAD�Ľ���Ƚ�ʱ��---//
	//getMATLABResult(1);//-----��MATLAB�Ľ���Ƚ�ʱ��---//0Ϊ0�������������ΪP�������
	//getCPPResult(1);//----��C++(Windows)�Ƚ�ʱ��---//0Ϊ0�������������ΪP�������

	double err;

	TMatrixD absErrMatrix;
	TMatrixD relErrMatrix;
	TMatrixD errMatrix;
	TMatrixD A;
	TMatrixD B;
	absErrMatrix.ResizeTo(1,row,1,nNodes+nColumns);
	relErrMatrix.ResizeTo(1,row,1,nNodes+nColumns);
	errMatrix.ResizeTo(1,row,1,nNodes+nColumns);
	A.ResizeTo(1,row,1,nNodes+nColumns);
	B.ResizeTo(1,row,1,nNodes+nColumns);

	absErrMatrix=resultMatrix-resultMatrix0;
	absErrMatrix.Abs();

	A=absErrMatrix;
	B=resultMatrix0;
	B.Abs();
	B=B+0.000001;
	ElementDiv(A,B);
	relErrMatrix=A;

	double tol=1.0e-4;
	for(int i=1;i<=row;i++)
	{
		for(int j=1;j<=nNodes+nColumns;j++)
		{
			errMatrix(i,j)=min(absErrMatrix(i,j),relErrMatrix(i,j));

		}
	}
	double maxErr=absErrMatrix.Max();
	cout<<setiosflags(ios::scientific)<<"MaxError"<<" = "<<errMatrix.Max()<<"--";
	
	if(maxErr<1.0e-2)
	{
		cout<<"case"<<setiosflags(ios::left)<<setw(3)<<caseNumber<<" is ok  !"<<'\t';
	}
	else
	{
		cout<<"case"<<setiosflags(ios::left)<<setw(3)<<caseNumber<<" is fail!"<<'\t';
	}

}

void EMTP::display()//��ʾ����
{
	cout<<deltaT<<endl;
	cout<<finishTime<<endl;
	cout<<counter<<endl;
	cout<<curTime<<endl;
	cout<<nNodes<<endl;
	cout<<nBranch;//֧·����
	timeArray.Print();//ʱ������
	nodeVoltageVec.Print();//�ڵ��ѹ����
	nodeVoltageMatrix.Print();//�ڵ��ѹ����
	branchCurrentMatrix.Print();//֧·��������
	conductanceMatrix.Print();//�ڵ㵼�ɾ���
	resistanceMatrix.Print();//�ڵ��迹����
	nodeNortonEquivalentCurrentArray.Print();//�ڵ�ŵ�ٵ�ֵ��������
}
//��ò���ϵͳ�������
void EMTP::getCaseMsrDfnMatrix(int caseNumber)
{
	caseMsrDfnMatrix.ResizeTo(2000,3);
	char ch[100];     
	string dir1;
#ifdef WIN32
	sprintf(ch,".\\caseDfn\\case%dmsr.txt",caseNumber);    
#else
	sprintf(ch,"..//caseDfn//case%dmsr.txt",caseNumber);    
#endif
	dir1=ch;
	ifstream infile(dir1.c_str(),ios::in);
	nMsrComps=0;
	while(!infile.eof())
	{
		infile>>caseMsrDfnMatrix(nMsrComps,0);
		if(infile.fail()) //ֻ���� >> ֮������жϣ�������Ч
		{				
			caseExistFlag=0;
			cout<<"Measure file for case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//ֱ�ӷ���
			break;
		}

		infile>>caseMsrDfnMatrix(nMsrComps,1);
		infile>>caseMsrDfnMatrix(nMsrComps,2);
		nMsrComps++;
	}
	caseMsrDfnMatrix.ResizeTo(nMsrComps,3);
	infile.close();      
}
//����ϵͳ���뺯��
void EMTP::	specifyMsrSystem()
{
	getCaseMsrDfnMatrix(caseNumber);
	int type=0;

	for(int i=0;i<nMsrComps;i++)
	{   
		type=(int)caseMsrDfnMatrix(i,0);
		MeasureComponent * com = NULL;

		switch (type)
		{
		case 1: 
			com = new NodeVoltageMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//�ڵ��ѹ����Ԫ��
			msrComps->push_back(com);
			break;
		case 2:
			com = new BranchVoltageMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//֧·��ѹ����Ԫ��
			msrComps->push_back(com);
			break;
		case 3:
			com = new BranchCurrentMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//֧·��������Ԫ��
			msrComps->push_back(com);
			break;
		case 4:
			com = new InducGenWrMsr(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//֧·��������Ԫ��
			msrComps->push_back(com);
			break;
		case 5:
			com = new InducGenAngleMsr(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//֧·��������Ԫ��
			msrComps->push_back(com);
			break;

		default:
			cout<<"Warning,There is no measure component with type"<<type<<endl;
		}
	}

}
//��ÿ���ϵͳ�������
void EMTP::getCaseCtrlDfnMatrix(int caseNumber)
{
	caseCtrlDfnMatrix.ResizeTo(2000,10);
	char ch[100];     
	string dir1;
#ifdef WIN32
	sprintf(ch,".\\caseDfn\\case%dctrl.txt",caseNumber);    
#else
	sprintf(ch,"..//caseDfn//case%dctrl.txt",caseNumber);    
#endif
	dir1=ch;
	ifstream infile(dir1.c_str(),ios::in);
	nCtrlBranches=0;
	while(!infile.eof())
	{
		infile>>caseCtrlDfnMatrix(nCtrlBranches,0);
		if(infile.fail()) //ֻ���� >> ֮������жϣ�������Ч
		{				
			caseExistFlag=0;
			cout<<"Control file for case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//ֱ�ӷ���
			break;
		}

		for(int j=1;j<10;j++)
		{
			infile>>caseCtrlDfnMatrix(nCtrlBranches,j);
		}

		nCtrlBranches++;
	}
	caseCtrlDfnMatrix.ResizeTo(nCtrlBranches,10);
	infile.close();      
}
//����ϵͳ���뺯��
void EMTP::	specifyCtrlSystem()
{
	getCaseCtrlDfnMatrix(caseNumber);
	int type=0;
	nCtrlNodes=0;//����ϵͳ�ܽڵ����
	int sigmaCtrlCompCount=0;
	int PICtrlCompCount=0;
	int PCtrlCompCount=0;
	int constantCount=0;
	int pulseGenCount=0;
	int triangleGenCount=0;
	int sinGenCount=0;
	int T2DTransCount = 0;
	int D2TTransCount = 0;
	int S2RTransCount=0;
	int R2STransCount=0;
	int PRCoordinateCount=0;
	int ComparatorCount=0;
	int LimiterCount=0;

	for(int i=0;i<nCtrlBranches;i++)
	{   
		type=(int)caseCtrlDfnMatrix(i,0);
		CtrlComponent * com = NULL;

		switch (type)
		{
		case 1: 
			com = new SigmaCtrlComp(++sigmaCtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,3));
			break;
		case 2:
			com = new PICtrlComp(++PICtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//PIԪ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			break;
		case 3:
			com = new ConstantCtrlComp(++constantCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2));//��������Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 4:
			com = new PulseGenComp(++pulseGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//��������Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 5:
			com = new TriangleGenComp(++triangleGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//��������Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 6:
			com = new SinGenComp(++sinGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//��������Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 7: 
			com = new PCtrlComp(++PCtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			break;
		case 8:
			com = new T2DTrans(++T2DTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 9:
			com = new D2TTrans(++D2TTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 10:
			com = new S2RTrans(++S2RTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 11:
			com = new R2STrans(++R2STransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 12:
			com = new PRCoordinate(++PRCoordinateCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<5;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 13:
			com = new Comparator(++ComparatorCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<5;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 14:
			com = new Limiter(++LimiterCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//�Ӻ���Ԫ��
			ctrlBranches->push_back(com);
			for (int j=1;j<3;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;


		default:
			cout<<"Warning,There is no control component with type"<<type<<endl;
		}
	}

	nodeCalMark = new int[nCtrlNodes];
	ctrlNodeValue = new double[nCtrlNodes];
	branchCalMark = new int[nCtrlBranches];

	for (int i=0;i<nCtrlBranches;i++)
	{
		branchCalMark[i] = 0;
	}
	for (int i=0;i<nCtrlNodes;i++)
	{
		nodeCalMark[i] = 0;
	}

	CtrlComponent* tempCom;
	for (int k=0;k<nCtrlBranches;k++)
	{
		tempCom = ctrlBranches->at(k);
		tempCom->setDeltaT(deltaT);
	}

	//// xuyin, 20121013
	//int r = (int)((finishTime+deltaT/2)/deltaT) + 1;
	//ctrlResultMatrix = new double* [r];
	//for (int k=0; k<r; k++)
	//	ctrlResultMatrix[k] = new double [nCtrlNodes];
	// xuyin, 20130224
	if ( nPWMConverter == 0 )
	{
		ctrlResultMatrix_2 = NULL;
	}
	else
	{
		int r2 = (int)((len_subinterval+deltaT/2)/deltaT)+1;
		ctrlResultMatrix_2 = new double* [r2];
		for (int k=0; k<r2; k++)
			ctrlResultMatrix_2[k] = new double [4*nPWMConverter];
	}

	// xuyin, 20121110
	// PI Controller
	nPIController = PICtrlCompCount;
	if (nPIController != 0)
	{
		PIControllerBranches = new int[nPIController];
		int k3=0;
		for (int k=0;k<nCtrlBranches;k++)
		{
			CtrlComponent* tempCtrlCom = ctrlBranches->at(k);
			if( tempCtrlCom->type == 2 )
			{
				PIControllerBranches[k3] = k;
				k3++;
			}
		}

		PIControllerInitialValue = new double[2*nPIController];
	}
	else
	{
		PIControllerBranches = NULL;
		PIControllerInitialValue = NULL;
	}
}


//������ϵͳ���������뵽����ϵͳ��
void EMTP::transferMeasurands()
{
	//// debug
	//cout << "����ǰ��ctrlNodeValue:" << endl;
	//for (int k=0; k<nCtrlNodes; k++)
	//	cout << ctrlNodeValue[k] << "  ";
	//cout << endl << endl;

	MeasureComponent* tempMsrComp;
	for (int i=0;i<nMsrComps;i++)
	{
		tempMsrComp = msrComps->at(i);
		ctrlNodeValue[tempMsrComp->ctrlNode-1] = tempMsrComp ->getMeasurands(nodeVoltageVec,branches,branchCurrentMatrix_1[counter-1]);
		nodeCalMark[tempMsrComp->ctrlNode-1] = 1;
	}

	//// debug
	//cout << "�����ctrlNodeValue:" << endl;
	//for (int k=0; k<nCtrlNodes; k++)
	//	cout << ctrlNodeValue[k] << "  ";
	//cout << endl << endl;
}

//������ϵͳ����
void EMTP::solveCtrlSystem()
{
	CtrlComponent* tempCtrlComp;
	int branchCaled=0;
	while (!(branchCaled==nCtrlBranches))
	{
		branchCaled=0;
		for (int i=0;i<nCtrlBranches;i++)
		{
			if (branchCalMark[i]==1)
			{
				branchCaled++;
			}
			else
			{
				tempCtrlComp = ctrlBranches->at(i);
				if (tempCtrlComp->checkCalCondition(nodeCalMark)==1)
				{
					tempCtrlComp->saveInNodeValue(ctrlNodeValue);
					tempCtrlComp->calculateOutputValue(curTime);
					tempCtrlComp->saveOutNodeValue(ctrlNodeValue);
					tempCtrlComp->markOutputNode(nodeCalMark);
					branchCalMark[i] = 1;
				}
			}
		}
	}
	for (int i=0;i<nCtrlBranches;i++)
	{
		branchCalMark[i]=0;
	}
	for (int i=0;i<nCtrlNodes;i++)
	{
		nodeCalMark[i]=0;
	}

	//// xuyin, 20121013
	//// �洢����ϵͳ�����
	//for (int k=0; k<nCtrlNodes; k++)
	//	ctrlResultMatrix[counter-1][k] = ctrlNodeValue[k];
	// xuyin, 20130224
	for (int k=0; k<4*nPWMConverter; k++)
	{
		ctrlResultMatrix_2[counter2][k] = ctrlNodeValue[ctrlNumArray_PWMConverter[k]-1];
	}

	// debug
	//cout << ctrlNodeValue[64] << "  " << ctrlNodeValue[65] << "  "
	//	<< ctrlNodeValue[66] << "  " << ctrlNodeValue[67] << endl; 
	//cout << ctrlNodeValue[40] << "  " << ctrlNodeValue[42] << "  "
	//	<< ctrlNodeValue[48] << endl; 
	//cout << ctrlNodeValue[145] << "  " << ctrlNodeValue[148] << endl; 
	//cout << ctrlNodeValue[96] << "  " << ctrlNodeValue[97] << "  "
	//	<< ctrlNodeValue[98] << endl; 
	//cout << ctrlNodeValue[123] << "  " << ctrlNodeValue[124] << "  "
	//	<< ctrlNodeValue[125] << endl; 
}


//�������źŴ���������ϵͳ
void EMTP::transferControlVariables()
{
	for(int k=0;k<nControlledBranches;k++)
	{
		Component * tempCom = branches->at(controlledBranches[k]);
		tempCom->setControlledVariable(ctrlNodeValue);
	}
}

void EMTP::transferControlVariablesForSwitch()
{
	for(int k=0;k<nControlledBranches;k++)
	{
		Component * tempCom = branches->at(controlledBranches[k]);
		tempCom->setControlledVariableForSwitch(ctrlNodeValue);
	}
}


//����ʼ��ʱ�Ŀ���ϵͳ����
void EMTP::solveInitCtrlSystem()
{
	CtrlComponent* tempCtrlComp;
	int branchCaled=0;
	while (!(branchCaled==nCtrlBranches))
	{
		//// debug
		//for (int i=0; i<nCtrlBranches; i++)
		//{
		//	cout << branchCalMark[i] << "  ";
		//}
		//cout << endl;

		branchCaled=0;
		for (int i=0;i<nCtrlBranches;i++)
		{
			if (branchCalMark[i]==1)
			{
				branchCaled++;
			}
			else
			{
				tempCtrlComp = ctrlBranches->at(i);
				if (tempCtrlComp->checkCalCondition(nodeCalMark)==1)
				{
					tempCtrlComp->saveInNodeValue(ctrlNodeValue);
					tempCtrlComp->calculateInitOutputValue(curTime);
					tempCtrlComp->saveOutNodeValue(ctrlNodeValue);
					tempCtrlComp->markOutputNode(nodeCalMark);
					branchCalMark[i] = 1;
				}
			}		
		}
	}
	for (int i=0;i<nCtrlBranches;i++)
	{
		branchCalMark[i]=0;
	}
	for (int i=0;i<nCtrlNodes;i++)
	{
		nodeCalMark[i]=0;
	}
	
	//// xuyin, 20121013
	//// �洢���
	//for (int k=0; k<nCtrlNodes; k++)
	//	ctrlResultMatrix[0][k] = ctrlNodeValue[k];
	// xuyin, 20130224
	for (int k=0; k<4*nPWMConverter; k++)
		ctrlResultMatrix_2[0][k] = ctrlNodeValue[ctrlNumArray_PWMConverter[k]-1];

	// debug
	//cout << ctrlNodeValue[64] << "  " << ctrlNodeValue[65] << "  "
	//	<< ctrlNodeValue[66] << "  " << ctrlNodeValue[67] << endl; 
	//cout << ctrlNodeValue[40] << "  " << ctrlNodeValue[42] << "  "
	//	<< ctrlNodeValue[48] << endl;
	// cout << ctrlNodeValue[145] << "  " << ctrlNodeValue[148] << endl;
	//cout << ctrlNodeValue[96] << "  " << ctrlNodeValue[97] << "  "
	//	<< ctrlNodeValue[98] << endl; 
	//cout << ctrlNodeValue[123] << "  " << ctrlNodeValue[124] << "  "
	//	<< ctrlNodeValue[125] << endl; 
}

/* PWM������ƽ��ģ����� */

// Ԥ�⿪�ض���ʱ��
// ver 0.1����ʱ��Ϊ����PWM���������ز�Ƶ����ͬ�����ز�����t=0ʱ��ȡ���ֵ
void EMTP::predictSwitchingInstants()
{
	// ����һ��������巨
	// Ԥ�����PWM�������Ŀ���ʱ��
	//int ptr = 0;
	//for (int k=0; k<nPWMConverter; k++) {
	//	Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
	//	tempCom->predictSwithcingInstants();
	//}
		
	// �����γɽڵ㵼�ɾ���LU�ֽ�
	//formConductanceMatrix();

	/* ����2��������ϵͳ */
	// double fre = 1000; // Hz, ������֪�ز�Ƶ��, case271
	double fre = 5000; // Hz, ������֪�ز�Ƶ��, case272
	int Nc = 1.0/(fre*deltaT); // ����һ���ز����ڰ�������ʱ��
	if ( counter % Nc == 1 ) // ������ز����ڿ�ʼʱ�̣�����Ԥ��
	{
		// ������ʱ�������ڴ洢Ԥ��ʱ������ϵͳ�õ��ĵ��Ʋ�˲ʱֵ���ز�˲ʱֵ
		double** PWMPredictMatrix = new double* [Nc+1];
		for (int k=0; k<=Nc; k++) {
			PWMPredictMatrix[k] = new double [4*nPWMConverter];
		}

		// ������ϵͳ��һ���������ڵĽ�
		solveCtrlSystemforPrediction(Nc, PWMPredictMatrix);

		// Ԥ�����PWM�������Ŀ���ʱ��
		int ptr = 0;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			tempCom->predictSwithcingInstants(PWMPredictMatrix, Nc, ptr);
		}
		
		// �����γɽڵ㵼�ɾ���LU�ֽ�
		formConductanceMatrix();
	}
}

// Ԥ�⿪��ʱ��ʱ��������ϵͳ����
void EMTP::solveCtrlSystemforPrediction(int Nc, double** PWMPredictMatrix)
{
	CtrlComponent* tempCtrlComp;

	// �洢����Ԫ��״̬
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->storeInternalVariables_pre();
	}

	// �洢��ǰʱ�����
	//// ���case271,��ʱ����ֻ��1��PWM������������֪��Ӧ����ϵͳ�Ľڵ���
	//for (int j=0; j<4; j++) {
	//	PWMPredictMatrix[0][j] = ctrlNodeValue[j];
	//}
	// ���case272
	for (int j=0; j<4; j++) {
		PWMPredictMatrix[0][j] = ctrlNodeValue[j+64];
	}

	// ������ϵͳ��һ�����������еĽ�
	for (int k=1; k<=Nc; k++) {
		int branchCaled=0;
		transferMeasurands();
		// ����ϵͳ���
		while (!(branchCaled==nCtrlBranches))
		{
			branchCaled=0;
			for (int i=0;i<nCtrlBranches;i++)
			{
				if (branchCalMark[i]==1)
				{
					branchCaled++;
				}
				else
				{
					tempCtrlComp = ctrlBranches->at(i);
					if (tempCtrlComp->checkCalCondition(nodeCalMark)==1)
					{
						tempCtrlComp->saveInNodeValue(ctrlNodeValue); // ��ȡ����ڵ��ֵ
						tempCtrlComp->calculateOutputValue(curTime+k*deltaT); // ����ģ��
						tempCtrlComp->saveOutNodeValue(ctrlNodeValue); // �������д��ctrlNodeValue������
						tempCtrlComp->markOutputNode(nodeCalMark); // ���ñ�־λ����ʾ����ź�ֵ��ȷ��
						branchCalMark[i] = 1; // ���ñ�־λ����ʾ��ģ����ɼ���
					}
				}		
			}
		}

		// �����־λ
		for (int i=0;i<nCtrlBranches;i++)
		{
			branchCalMark[i]=0;
		}
		for (int i=0;i<nCtrlNodes;i++)
		{
			nodeCalMark[i]=0;
		}

		// �洢�����
		//// ���case271,��ʱ����ֻ��1��PWM������������֪��Ӧ����ϵͳ�Ľڵ���
		//for (int j=0; j<4; j++) {
		//	PWMPredictMatrix[k][j] = ctrlNodeValue[j];
		//}
		// ���case272
		for (int j=0; j<4; j++) {
			PWMPredictMatrix[k][j] = ctrlNodeValue[j+64];
		}
	}

	// �ָ�����Ԫ��״̬
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->restoreInternalVariables_pre();
	}
}

// У�����ض���ʱ��
// ver 0.1����ʱ��Ϊ����PWM���������ز�Ƶ����ͬ�����ز�����t=0ʱ��ȡ���ֵ
void EMTP::correctSwitchingInstants(double tol)
{
	int Nc = (len_subinterval+deltaT/2)/deltaT; // ����һ������У���������������ʱ����

	if ( counter % Nc == 1 ) // ��������������ʱ�̣�����У��
	{
		//// ������ʱ�������ڴ洢У��ʱ����ĵ��Ʋ�˲ʱֵ���ز�˲ʱֵ
		//double** PWMCorrectMatrix = new double* [Nc+1];
		//for (int k=0; k<=Nc; k++) {
		//	PWMCorrectMatrix[k] = new double [4*nPWMConverter];
		//}

		//// �ӿ���ϵͳ������л�ȡ��PWM��������صĲ���
		//for (int i=0; i<=Nc; i++)
		//	for (int j=0; j<4*nPWMConverter; j++)
		//		PWMCorrectMatrix[i][j] = ctrlResultMatrix[counter-Nc-1+i][ctrlNumArray_PWMConverter[j]-1];

		////debug
		//cout << "PWMCorrectMatrix:" << endl;
		//for (int i=0; i<=1; i++) {
		//	for (int j=0; j<4*nPWMConverter; j++)
		//		cout << PWMCorrectMatrix[i][j] << "\t";
		//	cout << endl;
		//}
		//cout << endl;
		//cout << "ctrlResultMatrix_2:" << endl;
		//for (int i=0; i<=1; i++) {
		//	for (int j=0; j<4*nPWMConverter; j++)
		//		cout << ctrlResultMatrix_2[i][j] << "\t";
		//	cout << endl;
		//}
		//cout << endl;

		// �������PWM�������Ŀ���ʱ��
		int ptr = 0;
		int changeOrNot = 0;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			//if ( tempCom->correctSwithcingInstants(PWMCorrectMatrix, tol, Nc, ptr) == 1 )
			if ( tempCom->correctSwithcingInstants(ctrlResultMatrix_2, tol, Nc, ptr) == 1 )
				changeOrNot = 1;
		}

		// У��
		if ( changeOrNot == 1 ) // ����ҪУ������
		{
			// �����γɽڵ㵼�ɾ���
			formConductanceMatrix();

			// ����ʱ��ص����������ʼʱ��
			counter = counter - Nc;
			curTime = curTime - Nc * deltaT;
			counter2 = 0;

			// �ָ�����Ԫ��״̬
			for (int k=0; k<nBranch; k++) {
				Component* tempCom = branches->at(k);
				tempCom->restoreInternalVariables();
			}

			// �ָ�����Ԫ��״̬
			for (int k=0; k<nCtrlBranches; k++) {
				CtrlComponent* tempCtrl = ctrlBranches->at(k);
				tempCtrl->restoreInternalVariables();
			}
		}
		else // ������У������
		{
			// �洢��ǰʱ�̵���Ԫ��״̬
			for (int k=0; k<nBranch; k++) {
				Component* tempCom = branches->at(k);
				tempCom->storeInternalVariables();
			}

			// �洢����Ԫ��״̬
			for (int k=0; k<nCtrlBranches; k++) {
				CtrlComponent* tempCtrl = ctrlBranches->at(k);
				tempCtrl->storeInternalVariables();
			}

			for (int k=0; k<4*nPWMConverter;k++)
			{
				ctrlResultMatrix_2[0][k]=ctrlResultMatrix_2[counter2][k];
			}
			counter2 = 0;

			// ������һ�������ڵ�Ԥ��
			// predictSwitchingInstants();

			// debug
			// cout << "\r" << curTime - len_subinterval << endl;
		}
	}
}

// ===============================================================
// xuyin, 20130224
// ===============================================================

// ���ϵͳ��һ���ֶ��������ϵĽ⣬������tau_new=G(tau)
void EMTP::solveG(double * tau_new)
{
	nSolveG++;
	// ���һ���ֶ��������ϵĽ�
	for (int k=0; k<Ns; k++)
	{
		curTime += deltaT;
		counter += 1;
		counter2 += 1;
		
		if (counter%Nsave==1 || Nsave==1)
		{
			int savecounter = (int) (counter/Nsave)+1;
			if (Nsave==1)
			{
				savecounter-=1;
			}
			timeArray(savecounter) = curTime;
		}

		transferControlVariables(); //�������ź�ֵ����������ϵͳ		
		calculateBranchNortonEquivalentCurrent(curTime); //����֧·ŵ�ٵ�Ч����
		formNodeNortonEquivalentCurrentArray(); //�γɽڵ�ŵ�ٵ�Ч��������
		solveNodeVoltageEquation(); //���ڵ��ѹ����
		saveNodeVoltage(); //����ڵ��ѹ����,�����������ȫ�ֵ�ѹ�����У����֧·�޹ء�
		calculateBranchVoltage(); //����֧·��ѹ������
		calculateBranchCurrent(); //����֧·����
		saveBranchCurrent(); //����֧·����			
		saveMahineWr(); //������ת��
		transferMeasurands(); //������ϵͳ����ֵ���䵽����ϵͳ
		solveCtrlSystem(); //������ϵͳ
	}

	// ���㿪�ض���ʱ��
	int ptr = 0;
    double tao_s_new[3];
	for (int kk=0; kk<nPWMConverter; kk++)
	{
		for (int k=0; k<3; k++)
		{
			tao_s_new[k] = 0.0;
			for (int i=1; i<=Ns; i++)
			{
				if (ctrlResultMatrix_2[i-1][ptr+k] >= ctrlResultMatrix_2[i-1][ptr+3]
				&& ctrlResultMatrix_2[i][ptr+k] >= ctrlResultMatrix_2[i][ptr+3]) 
				{
					tao_s_new[k] += 1;
				}
				else if (ctrlResultMatrix_2[i-1][ptr+k] < ctrlResultMatrix_2[i-1][ptr+3]
				&& ctrlResultMatrix_2[i][ptr+k] < ctrlResultMatrix_2[i][ptr+3]) 
				{
					// tao_s_new[k] += 0;
				}
				else
				{
					double temp1 = ctrlResultMatrix_2[i][ptr+k]-ctrlResultMatrix_2[i][ptr+3];
					double temp2 = ctrlResultMatrix_2[i-1][ptr+k]-ctrlResultMatrix_2[i-1][ptr+3];
					tao_s_new[k] += 0.5*(1+(temp1+temp2)/(fabs(temp1)+fabs(temp2)));
				}
			}
			tao_s_new[k] /= Ns;
			if (tao_s_new[k]>1)
			{
				tao_s_new[k]=1;
			}
			if (tao_s_new[k]<0)
			{
				tao_s_new[k]=0;
			}
			tau_new[3*kk+k] = tao_s_new[k];
		}
		ptr += 4;
	}
}

// �ж��Ƿ���Ҫ����
int EMTP::check_tau(double * tau_new, double tol)
{
	// ����������
	double maxError = 0.0;
	for (int k=0; k<3*nPWMConverter; k++)
		if ( fabs(tau_new[k]-tau_s[cnt_p][k]) > maxError )
		{
			//cout<<k<<":"<<tau_new[k]-tau_s[cnt_p][k]<<endl;
			maxError = fabs(tau_s[cnt_p][k]-tau_new[k]);
		}
	if (maxError > tol)
	{
		return 1;
	} 
	else
	{	
		for (int k=0; k<3*nPWMConverter; k++)
		{
			tau_s[cnt_p][k] = tau_new[k];
		}
		return 0;
	}
}

// Ԥ��
void EMTP::predict_tau()
{
	//// �ͽ�Ԥ�ⷨ
	//for (int k=0; k<3*nPWMConverter; k++)
	//{
	//	tau_s[cnt_p][k] = tau_s[cnt_p-1][k];
	//}

	// �������Ʒ�
	if (cnt_p == 1)
	{
		for (int k=0; k<3*nPWMConverter; k++)
			tau_s[cnt_p][k] = tau_s[cnt_p-1][k];
	} 
	else
	{
		for (int k=0; k<3*nPWMConverter; k++)
			tau_s[cnt_p][k] = 2*tau_s[cnt_p-1][k] - tau_s[cnt_p-2][k];
	}

	//����Ԥ���ռ�ձ���0��1֮��
	for (int k=0; k<3*nPWMConverter; k++)
	{
		if (tau_s[cnt_p][k]>1)
		{
			tau_s[cnt_p][k]=1;
		}
		if (tau_s[cnt_p][k]<0)
		{
			tau_s[cnt_p][k]=0;
		}
	}

	// ����PWM��������ŵ�ٵ�ֵ���ɾ���
	double tao_s[3];
	for (int k=0; k<nPWMConverter; k++)
	{
		for (int j=0; j<3; j++)
		{
			tao_s[j] = tau_s[cnt_p][3*k+j];
		}
		Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
		tempCom->correntYne(tao_s);
	}

	// ���½ڵ㵼�ɾ���
	formConductanceMatrix();
}

// У��
void EMTP::correct_tau(double * tau_new)
{
	// ���Ե�����
	for (int k=0; k<3*nPWMConverter; k++)
	{
		tau_s[cnt_p][k] = tau_s[cnt_p][k] + 0.6 * (tau_new[k] - tau_s[cnt_p][k]);
	}

	// ����PWM��������ŵ�ٵ�ֵ���ɾ���
	double tao_s[3];
	for (int k=0; k<nPWMConverter; k++)
	{
		for (int j=0; j<3; j++)
		{
			tao_s[j] = tau_s[cnt_p][3*k+j];
		}
		Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
		tempCom->correntYne(tao_s);
	}

	// ���½ڵ㵼�ɾ���
	formConductanceMatrix();
}

// ����״̬
void EMTP::store_state()
{
	// �洢��ǰʱ�̵���Ԫ��״̬
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->storeInternalVariables();
	}

	// �洢����Ԫ��״̬
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->storeInternalVariables();
	}

	for (int k=0; k<4*nPWMConverter;k++)
	{
		ctrlResultMatrix_2[0][k]=ctrlResultMatrix_2[counter2][k];
	}
	counter2 = 0;
}

// �ָ�״̬
void EMTP::restore_state()
{
	// ����ʱ��ص����������ʼʱ��
	counter = counter - Ns;
	curTime = curTime - Ns * deltaT;
	counter2 = 0;

	// �ָ�����Ԫ��״̬
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->restoreInternalVariables();
	}

	// �ָ�����Ԫ��״̬
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->restoreInternalVariables();
	}
}
// ===============================================================