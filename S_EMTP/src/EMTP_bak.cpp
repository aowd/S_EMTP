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

void EMTP::initializeSystem() // 初始化函数，0启动或从PSCAD启动
{
/////////////////////////////////////////
// edited by xuyin, 20130119
// 定义INITIAL_VI时采用处理过的PSCAD数据文件进行初始化
// 处理后的数据文件有2个，分别存储电压电流，命名为case***_V.dat和case***_I.dat
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

		// 读取电压初值文件
		double t; // 第一列为时间
		char chDir[100];   
		string caseDir;
#ifdef WIN32
		sprintf(chDir,".\\result\\pscadData\\case%d_V.dat",caseNumber);   
#else
		sprintf(chDir,"..//result//pscadData//case%d_V.dat",caseNumber);   
#endif
		caseDir=chDir;
		ifstream infile(caseDir.c_str(),ios::in);

		// 过滤掉startTime以前的行
		infile >> t;
		double tmp = 0;
		while ( fabs(t-startTime) > initDeltaT/2 ) {
			for(int k=0;k<nNodes;k++)
				infile >> tmp;
			infile >> t;
		}

		// 形成电压初值向量
		for(int k=1;k<=nNodes;k++)
			infile>>initialVoltageVec(k);

		infile.close();
		infile.clear();

		// 读取电流初值文件
#ifdef WIN32
		sprintf(chDir,".\\result\\pscadData\\case%d_I.dat",caseNumber);   
#else
		sprintf(chDir,"..//result//pscadData//case%d_I.dat",caseNumber);   
#endif
		caseDir=chDir;
		infile.open(caseDir.c_str(),ios::in);

		// 过滤掉startTime以前的行
		infile >> t;
		while ( fabs(t-startTime) > initDeltaT/2 ) {
			for(int k=0;k<nColumns;k++)
				infile >> tmp;
			infile >> t;
		}

		// 形成电流初值向量
		for(int k=1;k<=nColumns;k++)
			infile>>initialCurrentVec(k);

		infile.close();
		infile.clear();

		initialVoltageVec *= 1000.0;
		initialCurrentVec *= 1000.0;
	}

	// 初始化节点电压
	nodeVoltageVec = initialVoltageVec;

	// 初始化支路
	int ptr = 1;//指针
	for (int i = 1; i <= nBranch; i ++)
	{
		Component * tempCom = branches->at(i-1);
		tempCom->initializeBranch(initialVoltageVec,initialCurrentVec,ptr);
	}

	// 存储初始化结果
	for (int i=1; i<=nNodes; i++)
		nodeVoltageMatrix_1[0][i-1] = initialVoltageVec(i);

	for (int i=1; i<=nColumns; i++)
		branchCurrentMatrix_1[0][i-1] = initialCurrentVec(i);

#else
	int nCol=0;//初值向量的维数
	nCol = nNodes + nColumns;

	int units=0;//初值向量维数的个位数
	units=nCol%10;

	int nFiles=0;//PSCAD.out文件的个数，每个,out文件中存有11列数，第一列是时间，后面10列是10个变量的数据
	nFiles=(nCol-units)/10;

	int nVecIndex=1;
	TVectorD initiaValVec,branchCurrentVec;
	initiaValVec.ResizeTo(1,nCol);
	branchCurrentVec.ResizeTo(1,nColumns);

	if(initializeOption==0)//如果option=0
	{
		initiaValVec*=0.0;
		cout<<"S_EMTP begin with 0 state"<<endl;
	}
	else
	{
		cout<<"S_EMTP begin with PSCAD state"<<endl;
		//	step1:cal number of *.out files 
		if(units!=0)//如果个位不为0
		{
			nFiles=nFiles+1;
		}

		//	step2:scan all *.out files and read init data
		for (int fileIndex=1;fileIndex<=nFiles;fileIndex++)
		{
			int num1;
			int num2;
			double t;//第一列为时间
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
			// 过滤掉startTime以前的行
			if (fileIndex<nFiles)//不是最后一个文件，直接读走第一行的10个数据（除去时间变量）
			{
				infile >> t;
				// 过滤掉startTime以前的数据
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
			else//到了最后一个文件，只需要读走units个变量即可
			{
				if (units==0)
				{
					units=10;
				}

				infile>>t;
				// 过滤掉startTime以前的数据
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

	int ptr = 1;//指针
	for (int i = 1; i <= nBranch; i ++)
	{
		Component * tempCom = branches->at(i-1);
		tempCom->initializeBranch(nodeVoltageVec,branchCurrentVec,ptr);
	}
#endif

	////////////////////////////////////////////////////////////////////
	// PWM变流器平均模型开关动作时刻初始化,xuyin,20121011
	// edit by xuyin, 20130119, 修改使之适应包含多个初始化文件的情况
	if (nPWMConverter != 0) {
		int nCol_PWM = 3*nPWMConverter;
		int units_PWM = nCol_PWM%10;
		int nFiles_PWM = (nCol_PWM-units_PWM)/10;
		if(units_PWM!=0)
		{
			nFiles_PWM = nFiles_PWM + 1;
		}

		// 计算列数和行数
		int c = nCol_PWM + 1;
		int r = (len_subinterval+initDeltaT/2)/initDeltaT+1;

		// 形成PWM变流器平均模型初始化矩阵
		double** PWMInitialMatrix = new double* [r];
		for (int k=0; k<r; k++) {
			PWMInitialMatrix[k] = new double [c];
		}

		// 读入数据
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

			double t; // 时间
			if (fileIndex < nFiles_PWM)//不是最后一个文件，直接读走第一行的10个数据（除去时间变量）
			{
				// 过滤掉startTime以前的数据
				double tmp = 0;
				infile >> t;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给PWM初始化数组赋值
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

				// 过滤掉startTime以前的数据
				double tmp = 0;
				infile >> t;
				while ( fabs(t-startTime) > initDeltaT/2 ) {
					for(int k=0; k<units_PWM; k++)
						infile >> tmp;
					infile >> t;
				}

				// 给PWM初始化数组赋值
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

		// 初始化开关动作时刻
		ptr = 1;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			tempCom->initializeSwitchingInstants(PWMInitialMatrix, r-1, ptr);
		}

		// 初始化tau_s, xuyin, 20130224
		for (int k=0; k<nPWMConverter; k++)
		{
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			for (int j=0; j<3; j++)
			{
				tau_s[0][3*k+j] = tempCom->tao_s[j];
			}
		}

		// 删除临时数组
		for (int k=0; k<r; k++) {
			delete[] PWMInitialMatrix[k];
		}
		delete[] PWMInitialMatrix;
	}
	//////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////
	// 异步电机初始化,xuyin,20121226
	if (nIndGen != 0)
	{
		int nCol_Gen = 15*nIndGen;
		int units_Gen = nCol_Gen%10;
		int nFiles_Gen = (nCol_Gen-units_Gen)/10;
		if(units_Gen!=0)
		{
			nFiles_Gen = nFiles_Gen + 1;
		}

		// 电机初始化矩阵
		double* GenInitialMatrix[3];
		for (int k=0; k<3; k++)
		{
			GenInitialMatrix[k] = new double [nCol_Gen+1];
		}

		// 读入数据
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

			// 过滤掉startTime-2*deltaT以前的行
			double t; // 时间
			if (fileIndex < nFiles_Gen)//不是最后一个文件，直接读走第一行的10个数据（除去时间变量）
			{
				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给电机初始化数组赋值 
				for(int i=0;i<3;i++)//read 10 initdata
				{
					GenInitialMatrix[i][0] = t;
					for(int j=1; j<=10; j++)
						infile >> GenInitialMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			else//到了最后一个文件，只需要读走units个变量即可
			{
				if (units_Gen == 0)
				{
					units_Gen = 10;
				}

				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<units_Gen;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给电机初始化数组赋值 
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

		// 初始化电机
		ptr = 1; // 第一列为时间，故从1开始
		for (int k=0; k<nIndGen; k++) {
			Component* tempCom = branches->at(branchNumArray_IndGen[k]);
			tempCom->initializeGen(GenInitialMatrix, ptr);
		}

		// 删除临时数组
		for (int k=0; k<3; k++) {
			delete[] GenInitialMatrix[k];
		}

#ifdef WIND_VELOCITY_DATA_INPUT
		// 读入风速文件
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

			// 过滤掉startTime-2*deltaT以前的行
			double t; // 时间
			if (fileIndex < nFiles_Gen)//不是最后一个文件，直接读走第一行的10个数据（除去时间变量）
			{
				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给电机初始化数组赋值 
				for(int i=0;i<rows;i++)//read VwData
				{
					for(int j=0; j<10; j++)
						infile >> VwMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			else//到了最后一个文件，只需要读走units个变量即可
			{
				if (units_GenVw == 0)
				{
					units_GenVw = 10;
				}

				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime+2*initDeltaT) > initDeltaT/2 )
				{
					for(int k=0;k<units_Gen;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给电机初始化数组赋值 
				for(int i=0;i<rows;i++)//read 10 initdata
				{
					for(int j=0; j<units_GenVw; j++)
						infile >> VwMatrix[i][j+10*(fileIndex-1)];
					infile >> t;
				}
			}
			infile.close();
		}

		// 初始化电机
		ptr = 1; // 第一列为时间，故从1开始
		for (int k=0; k<nIndGen; k++) {
			Component* tempCom = branches->at(branchNumArray_IndGen[k]);
			tempCom->getWindVelocityData(VwMatrix,rows,ptr);
		}

		// 删除临时数组
		for (int k=0; k<rows; k++) {
			delete[] VwMatrix[k];
		}
#endif
	}	
	//////////////////////////////////////////////////////////////////////

	//测量系统初始化
	MeasureComponent* tempMsrCom;
	for(int i=0;i<nMsrComps;i++)
	{
		tempMsrCom = msrComps->at(i);
		tempMsrCom->initializeBranch();
	}

	//控制系统初始化
	///////////////////////////////////////////////////////////////////////
	// edited by xuyin, 20121109
	// 读取PI控制器初值文件
	if ( nPIController != 0 )
	{
		int nCol_PI = 2*nPIController;
		int units_PI = nCol_PI%10;
		int nFiles_PI = (nCol_PI-units_PI)/10;
		if(units_PI!=0)
		{
			nFiles_PI = nFiles_PI + 1;
		}

		// PI控制器初始化向量
		PIControllerInitialValue = new double [nCol_PI+1];

		// 读入数据
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

			// 过滤掉startTime以前的行
			double t; // 时间
			if (fileIndex < nFiles_PI)//不是最后一个文件，直接读走第一行的10个数据（除去时间变量）
			{
				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 )
				{
					for(int k=0;k<10;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给PI控制器初始化数组赋值 
				PIControllerInitialValue[0] = t;
				for(int j=1; j<=10; j++)
					infile >> PIControllerInitialValue[j+10*(fileIndex-1)];
			}
			else//到了最后一个文件，只需要读走units个变量即可
			{
				if (units_PI == 0)
				{
					units_PI = 10;
				}

				infile >> t;
				// 过滤掉startTime以前的数据
				double tmp = 0;
				while ( fabs(t-startTime) > initDeltaT/2 )
				{
					for(int k=0;k<units_PI;k++)
						infile >> tmp;
					infile >> t;
				}

				// 给PI控制器初始化数组赋值 
				PIControllerInitialValue[0] = t;
				for(int j=1; j<=units_PI; j++)
					infile >> PIControllerInitialValue[j+10*(fileIndex-1)];
			}
			infile.close();
		}
	}
	
	// 初始化
	CtrlComponent* tempCtrlCom;
	for(int i=0;i<nCtrlBranches;i++)
	{
		tempCtrlCom = ctrlBranches->at(i);
		tempCtrlCom->initializeCtrlBranch();
		tempCtrlCom->calculateCtrlEquivalentParameter();
	}
	
	// 给PI控制器赋初始值 
	ptr = 1;
	for (int i=0; i<nPIController; i++)
	{
		// 赋初值
		tempCtrlCom = ctrlBranches->at(PIControllerBranches[i]);
		tempCtrlCom->initializePICtrl(PIControllerInitialValue, ptr);
	}
	////////////////////////////////////////////////////////////////////////

	formConductanceMatrix();//形成节点导纳矩阵
	transferMeasurands();
	solveInitCtrlSystem();//求解初始化方程，得到控制信号的初始量值
	transferControlVariables();//将控制信号值传输到电气系统

	if (!initializeOption)
	{
		calculateBranchNortonEquivalentCurrent(curTime);//计算支路诺顿等效电流
		formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
		solveNodeVoltageEquation();//求解节点电压方程
		saveNodeVoltage();//保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。
		calculateBranchVoltage();//计算支路电压并保存
		calculateBranchCurrent();//计算支路电流
		saveBranchCurrent();//保存支路电流
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// edit by xuyin, 20121123
	// 迭代校正专用
	// 存储当前时刻电气元件状态
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->storeInternalVariables();
	}
	// 存储控制元件状态
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->storeInternalVariables();
	}
	/////////////////////////////////////////////////////////////////////////////////////////
}

void EMTP::saveBranchCurrent()//保存支路电流
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

void EMTP::saveMahineWr()//保存电机转速
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

void EMTP::formConductanceMatrix()//形成节点导纳阵
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

void EMTP::calculateBranchNortonEquivalentCurrent(double time)//计算支路诺顿等效电流
{
	for(int k=0;k<nBranch_NEC;k++)
	{
		Component * tempCom = branches->at(branchNumArray_NEC[k]);
		tempCom->calculateNortonEquivalentCurrent(time);
	}
}

void EMTP::formNodeNortonEquivalentCurrentArray()
{
	//形成节点诺顿等效电流向量
	for(int i=1;i<=nNodes;i++)
		nodeNortonEquivalentCurrentArray(i) = 0;

	for(int k=0;k<nBranch_NEC;k++)
	{
		Component * tempCom = branches->at(branchNumArray_NEC[k]);
		tempCom->formNodeNortonEquivalentCurrentArray(nodeNortonEquivalentCurrentArray);
	}
}

void EMTP::calculateBranchVoltage(){
	//计算支路电压
	for(int k=0;k<nBranch;k++)
	{
		Component * tempCom = branches->at(k);
		tempCom->readNodeVoltage(nodeVoltageVec);
		tempCom->calculateBranchVoltage();
	}
}

void EMTP::calculateBranchCurrent(){
	//计算支路电流
	for(int k=0;k<nBranch;k++)
	{
		Component * tempCom = branches->at(k);
		tempCom->calculateBranchCurrent();
	}
}

//////////////////////////////////////////////////////////////////////////
//                      开关处理相关函数起始点                          //
//////////////////////////////////////////////////////////////////////////

//开关处理主函数
void EMTP::switchTreatment()
{
	//确定开关动作时刻
	bool switchOrNot = 0;//标志位，是否有开关动作
	bool checkSwitchOrNot = 1;//标志位，是否需要进行开关过程的判断和处理，初始值为1（每个时步开始时总是要检测）
	double ratio;//记录最早动作的开关的插值比
	double timeToSwitch = curTime-2*deltaT;//记录开关动作时间
	int switchID;//记录最早动作的开关的编号
	int switchCounter = 0;//统计需要动作的开关数
	int switchMode = 0;//0为自然开断，1为控制开断
	double maxRatio = 1;//为确保开关动作时刻都在curTime-deltaT时刻之前而设置的阈值

	while (checkSwitchOrNot)
	{
		//每次开过过程处理的开始，初始化变量
		switchID = -1;checkSwitchOrNot=0;switchMode=0;

		//确定需要动作的开关，确定首个动作开关的插值比以及开关动作时间
		checkSwitch(ratio,switchID,maxRatio,switchCounter);
		
		//如果有开关动作，进行相应处理
		if (switchID!=-1)
		{
			timeToSwitch = timeToSwitch + ratio*deltaT;//如果有开关动作，更新timeToSwitch
			maxRatio = (curTime-deltaT-timeToSwitch)/deltaT;

			switchOrNot = 1;
			checkSwitchOrNot = 1;//仍然需要进一步判断开关动作
			doSwitch(ratio,switchCounter, switchMode);//动作需要动作的开关，重新形成导纳矩阵
	
			linearInt(ratio);//(1)线性插值，求解timeToSwitch时刻的解
			if (switchMode!=0)//控制开断过程，保持动作点诺顿等效电流值不变，求解开关动作t+时刻状态
			{
				networkSolution();//求解t+时刻状态
				//cout<<"内层循环开始"<<endl;
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
				//cout<<"内层循环结束"<<endl;
			}
			calculateBranchNortonEquivalentCurrent_forBasicComp(timeToSwitch+deltaT);//计算支路诺顿等效电流
			networkSolution();//(2)求解timeToSwitch+deltaT时刻的解
		}		
	}

	//消除Chatter和存储结果
	if (switchOrNot)
	{
		//removeChatter(timeToSwitch);

		linearInt((curTime-deltaT-timeToSwitch)/deltaT);

		transferMeasurands();//将电气系统测量值传输到控制系统

		solveCtrlSystem();

		transferControlVariablesForSwitch();//将控制信号值传输至电气系统

		saveNewResult(counter-1);
	}
	else
	{
		/************************************************************************/
		/*    检查是否有Inverter状态发生变化，若有变化则重新形成导纳矩阵        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/
	}
} 

//判断需要动作的开关，以及首先动作的开关对应的ratio
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

//动作需要动作的开关，重新形成导纳矩阵
void EMTP::doSwitch(double ratio,int switchCounter, int& switchMode)
{
	
	//变换需要动作的开关
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

	//重新对导纳矩阵进行lu分解
	lu_conductanceMatrix->SetMatrix(conductanceMatrix);
	lu_conductanceMatrix->Decompose();
}

void EMTP::removeChatter(double timeToSwitch)
{
	if(timeToSwitch<curTime-3.0/2.0*deltaT)//如果最后的开关动作点在前半个步长内
	{
		//调用半步长插值算法，得到curTime-deltaT时刻的解
		linearInt(0.5);
		calculateBranchNortonEquivalentCurrent_forBasicComp(timeToSwitch+3.0/2.0*deltaT);//计算支路诺顿等效电流
		networkSolution();
		linearInt((curTime-timeToSwitch-deltaT*3.0/2.0)/deltaT);//
		//更新输出矩阵中curTime-deltaT时刻的存储结果
		saveNewResult(counter-1);
		updateResult(3);

		/************************************************************************/
		/*    检查是否有Inverter状态发生变化，若有变化则重新形成导纳矩阵        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/
	}
	else//如果最后的开关动作点在后半个步长内
	{
		//（1）首先插值回到整时步点curTime-deltaT，并更新输出矩阵中的存储结果
		linearInt((curTime-deltaT-timeToSwitch)/deltaT);
		saveNewResult(counter-1);

		/************************************************************************/
		/*    检查是否有Inverter状态发生变化，若有变化则重新形成导纳矩阵        */
		/************************************************************************/
		inverterTreatment();
		/************************************************************************/

		//（2）向前积分一次，得到curTime时刻的解
		updateResult(3);
		calculateBranchNortonEquivalentCurrent(curTime);//计算支路诺顿等效电流
		networkSolution();
		//消除Chatter前将_1变量的值存入_2变量中进行保存
		updateResult(1);
		//（3）调用半步长插值算法，去掉curTime时刻的解中的Chatter
		linearInt(0.5);
		calculateBranchNortonEquivalentCurrent_forBasicComp(curTime+deltaT/2.0);//计算支路诺顿等效电流
		networkSolution();
		//cout<<"第2次半步长插值："<<endl;
		linearInt(0.5);
		//更将curTime时刻的解存入输出矩阵中
		saveNewResult(counter);
		//消除Chatter后用_2变量的值还原_1变量的值
		updateResult(2);
		//（4）时间向前推进一个时步，进入新一轮的开关判断
		advanceTime();
		switchTreatment();
	}
}

//一次线性插值
void EMTP::linearInt(double ratio)
{
	Component* tempCom;

	//节点电压的插值
	for (int k=1;k<=nNodes;k++)
	{
		nodeVoltageVec(k) = nodeVoltageVec_1(k) + ratio*(nodeVoltageVec(k)-nodeVoltageVec_1(k));
	}

	//各支路电压电流的插值
	for(int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		tempCom->interpolate(ratio);
	}
}

//计算支路诺顿等效电流
void EMTP::calculateBranchNortonEquivalentCurrent_forBasicComp(double time)
{
	Component * tempCom;
	for(int k=0;k<nBranch_NEC;k++)
	{
		tempCom = branches->at(branchNumArray_NEC[k]);
		if (!tempCom->isUserDef)//若不是自定义元件
		{
			tempCom->calculateNortonEquivalentCurrent(time);
		}
	}
}

//一次网络求解
void EMTP::networkSolution()
{
	formNodeNortonEquivalentCurrentArray();//形成节点诺顿等效电流向量
	solveNodeVoltageEquation();//求解节点电压方程
	calculateBranchVoltage();//计算支路电压
	calculateBranchCurrent();//计算支路电流
}

//更新开关处理过程中用于存储结果的变量
void EMTP::updateResult(int updateTypeNum)
{
	Component * tempCom;//指向元件的指针
	double temp;//临时变量，用于数据交换

	switch (updateTypeNum)
	{
	case 1://将_1变量的数值存入_2变量中
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
	case 2://将_2变量的数值存入_1变量中
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
		cerr<<"请输入正确的更新类型编号！"<<endl;
		exit(1);
	}
}

//将插值得到的结果存入结果矩阵中
void EMTP::saveNewResult(int counter)
{
	if (counter<finishTime/deltaT+2 && (counter%Nsave==1 || Nsave==1))
	{
		Component* tempCom;//临时变量，扫描支路时用
		int savecounter = (int)(counter/Nsave)+1;
		if (Nsave==1)
		{
			savecounter-=1;
		}

		//更新节点电压计算结果
		for(int i=1;i<=nNodes;i++)
		{
			//nodeVoltageMatrix(counter,i)=nodeVoltageVec(i);
			nodeVoltageMatrix_1[savecounter-1][i-1] = nodeVoltageVec(i);
		}
	
		//更新支路电流计算结果
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

//时间前进一个步长
void EMTP::advanceTime()
{
	curTime += deltaT;
	counter++;
	if (counter<finishTime/deltaT+2)
	{
		timeArray(counter) = curTime;
	}
}

//检查是否有Inverter状态发生变化，若有变化则重新形成导纳矩阵
void EMTP::inverterTreatment()
{
	//判断是否有Inverter状态变化
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

	//如果有Inverter状态发生变化，重新形成导纳矩阵
	if (changeOrNot)
	{
		lu_conductanceMatrix->SetMatrix(conductanceMatrix);
		lu_conductanceMatrix->Decompose();
	}
}

//////////////////////////////////////////////////////////////////////////
//              开关处理相关函数终止点                                  //
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
		if(infile.fail()) //只有在 >> 之后才能判断，否则无效
		{				
			caseExistFlag=0;
			cout<<"case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//直接返回
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
	nNodes=(int)caseDfnMatrix.GetSub(0,caseDfnMatrix.GetRowUpb(),1,2).Max();//节点总数

	rows=1;//时间，电压电流等全局数组的长度
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
			//case 14:	//浙大的控制方法
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
			//case 16://根据节点电压判断的算例
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
	timeArray.ResizeTo(1,rows);//时间数组
	nodeVoltageVec.ResizeTo(1,nNodes);//节点电压向量
	nodeVoltageVec_1.ResizeTo(1,nNodes);
	nodeVoltageVec_2.ResizeTo(1,nNodes);
	nodeVoltageMatrix.ResizeTo(1,rows,1,nNodes);//节点电压矩阵

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

	conductanceMatrix.ResizeTo(1,nNodes,1,nNodes);//节点导纳矩阵
	resistanceMatrix.ResizeTo(1,nNodes,1,nNodes);;//节点阻抗矩阵
	nodeNortonEquivalentCurrentArray.ResizeTo(1,nNodes);//节点诺顿等值电流向量
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
	/*                     开关处理相关变量和数组的创建                     */
	/************************************************************************/
	Component* tempCom;//临时变量，扫描支路时用

	//二极管相关
	nSwitch = diodeCount+newIGBTCount;
	if (nSwitch != 0)
	{
		diodeNumArray = new int[nSwitch];
		switchNumArray = new int[nSwitch];//需要动作的开关的编号
		switchRatioArray = new double[nSwitch];//各个需要动作的开关对应的ratio
		switchModeArray = new int[nSwitch];//各个开关动作采用的mode

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
	
	//Inverter相关
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

	//TimeSwitch相关
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
	/*                    诺顿等值电流计算相关数组的创建                    */
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
	/*                        自定义元件数组的创建                          */
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
	/*                 受控元件相关数组创建                      */
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
	/*          PWM变流器平均模型相关数组创建                    */
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
	/*                   异步电机相关数组创建                    */
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
	/*                  设置仿真步长                            */
	/************************************************************/
	for (int k=0;k<nBranch;k++)
	{
		tempCom = branches->at(k);
		tempCom->setDetalT(deltaT);
		tempCom->calculateNortonEquivalentResistance();
	}
}


void EMTP::solveNodeVoltageEquation()//求解节点电压方程
{
	bool ok;
	nodeVoltageVec_1 = nodeVoltageVec;
	nodeVoltageVec=lu_conductanceMatrix->Solve(nodeNortonEquivalentCurrentArray,ok);
}

void EMTP::saveNodeVoltage()//保存节点电压
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
	//第二种方案，输出支路的名字，对应支路电流		
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
	if(units!=0)//如果个位不为0
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
	//初始化参考结果数组
	int row=0;            
	row=timeArray.GetNrows();
	resultMatrix0.ResizeTo(1,row,1,nNodes+this->nColumns);
	resultMatrix0.Zero();
	getPSCADResult();//----和PSCAD的结果比较时用---//
	//getMATLABResult(1);//-----和MATLAB的结果比较时用---//0为0启动结果，其他为P启动结果
	//getCPPResult(1);//----和C++(Windows)比较时用---//0为0启动结果，其他为P启动结果

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

void EMTP::display()//显示函数
{
	cout<<deltaT<<endl;
	cout<<finishTime<<endl;
	cout<<counter<<endl;
	cout<<curTime<<endl;
	cout<<nNodes<<endl;
	cout<<nBranch;//支路总数
	timeArray.Print();//时间数组
	nodeVoltageVec.Print();//节点电压向量
	nodeVoltageMatrix.Print();//节点电压矩阵
	branchCurrentMatrix.Print();//支路电流矩阵
	conductanceMatrix.Print();//节点导纳矩阵
	resistanceMatrix.Print();//节点阻抗矩阵
	nodeNortonEquivalentCurrentArray.Print();//节点诺顿等值电流向量
}
//获得测量系统定义矩阵
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
		if(infile.fail()) //只有在 >> 之后才能判断，否则无效
		{				
			caseExistFlag=0;
			cout<<"Measure file for case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//直接返回
			break;
		}

		infile>>caseMsrDfnMatrix(nMsrComps,1);
		infile>>caseMsrDfnMatrix(nMsrComps,2);
		nMsrComps++;
	}
	caseMsrDfnMatrix.ResizeTo(nMsrComps,3);
	infile.close();      
}
//测量系统输入函数
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
			com = new NodeVoltageMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//节点电压测量元件
			msrComps->push_back(com);
			break;
		case 2:
			com = new BranchVoltageMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//支路电压测量元件
			msrComps->push_back(com);
			break;
		case 3:
			com = new BranchCurrentMsrComp(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//支路电流测量元件
			msrComps->push_back(com);
			break;
		case 4:
			com = new InducGenWrMsr(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//支路电流测量元件
			msrComps->push_back(com);
			break;
		case 5:
			com = new InducGenAngleMsr(i+1,(int)caseMsrDfnMatrix(i,1),(int)caseMsrDfnMatrix(i,2));//支路电流测量元件
			msrComps->push_back(com);
			break;

		default:
			cout<<"Warning,There is no measure component with type"<<type<<endl;
		}
	}

}
//获得控制系统定义矩阵
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
		if(infile.fail()) //只有在 >> 之后才能判断，否则无效
		{				
			caseExistFlag=0;
			cout<<"Control file for case "<<caseNumber<<" is not exist,to be continue ..."<<endl;

			return;//直接返回
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
//控制系统输入函数
void EMTP::	specifyCtrlSystem()
{
	getCaseCtrlDfnMatrix(caseNumber);
	int type=0;
	nCtrlNodes=0;//控制系统总节点个数
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
			com = new SigmaCtrlComp(++sigmaCtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4));//加和器元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,3));
			break;
		case 2:
			com = new PICtrlComp(++PICtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//PI元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			break;
		case 3:
			com = new ConstantCtrlComp(++constantCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2));//常量输入元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 4:
			com = new PulseGenComp(++pulseGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//常量输入元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 5:
			com = new TriangleGenComp(++triangleGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//常量输入元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 6:
			com = new SinGenComp(++sinGenCount,(int)caseCtrlDfnMatrix(i,1),caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//常量输入元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			break;
		case 7: 
			com = new PCtrlComp(++PCtrlCompCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3));//加和器元件
			ctrlBranches->push_back(com);
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,1));
			nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,2));
			break;
		case 8:
			com = new T2DTrans(++T2DTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 9:
			com = new D2TTrans(++D2TTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 10:
			com = new S2RTrans(++S2RTransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 11:
			com = new R2STrans(++R2STransCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<6;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 12:
			com = new PRCoordinate(++PRCoordinateCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),(int)caseCtrlDfnMatrix(i,5));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<5;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 13:
			com = new Comparator(++ComparatorCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),(int)caseCtrlDfnMatrix(i,3),(int)caseCtrlDfnMatrix(i,4),caseCtrlDfnMatrix(i,5),caseCtrlDfnMatrix(i,6));//加和器元件
			ctrlBranches->push_back(com);
			for (int j=1;j<5;j++)
			{
				nCtrlNodes = max(nCtrlNodes, (int)caseCtrlDfnMatrix(i,j));
			}
			break;
		case 14:
			com = new Limiter(++LimiterCount,(int)caseCtrlDfnMatrix(i,1),(int)caseCtrlDfnMatrix(i,2),caseCtrlDfnMatrix(i,3),caseCtrlDfnMatrix(i,4));//加和器元件
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


//将电气系统测量量输入到控制系统中
void EMTP::transferMeasurands()
{
	//// debug
	//cout << "量测前，ctrlNodeValue:" << endl;
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
	//cout << "量测后，ctrlNodeValue:" << endl;
	//for (int k=0; k<nCtrlNodes; k++)
	//	cout << ctrlNodeValue[k] << "  ";
	//cout << endl << endl;
}

//求解控制系统方程
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
	//// 存储控制系统求解结果
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


//将控制信号传输至电气系统
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


//求解初始化时的控制系统方程
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
	//// 存储结果
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

/* PWM变流器平均模型相关 */

// 预测开关动作时刻
// ver 0.1，暂时认为所有PWM变流器的载波频率相同，且载波均在t=0时刻取最大值
void EMTP::predictSwitchingInstants()
{
	// 方法一：线性外插法
	// 预测各个PWM变流器的开关时刻
	//int ptr = 0;
	//for (int k=0; k<nPWMConverter; k++) {
	//	Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
	//	tempCom->predictSwithcingInstants();
	//}
		
	// 重新形成节点导纳矩阵并LU分解
	//formConductanceMatrix();

	/* 方法2：求解控制系统 */
	// double fre = 1000; // Hz, 假设已知载波频率, case271
	double fre = 5000; // Hz, 假设已知载波频率, case272
	int Nc = 1.0/(fre*deltaT); // 计算一个载波周期包含几个时步
	if ( counter % Nc == 1 ) // 如果是载波周期开始时刻，进行预测
	{
		// 创建临时数组用于存储预测时求解控制系统得到的调制波瞬时值和载波瞬时值
		double** PWMPredictMatrix = new double* [Nc+1];
		for (int k=0; k<=Nc; k++) {
			PWMPredictMatrix[k] = new double [4*nPWMConverter];
		}

		// 求解控制系统在一个开关周期的解
		solveCtrlSystemforPrediction(Nc, PWMPredictMatrix);

		// 预测各个PWM变流器的开关时刻
		int ptr = 0;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			tempCom->predictSwithcingInstants(PWMPredictMatrix, Nc, ptr);
		}
		
		// 重新形成节点导纳矩阵并LU分解
		formConductanceMatrix();
	}
}

// 预测开关时刻时，求解控制系统方程
void EMTP::solveCtrlSystemforPrediction(int Nc, double** PWMPredictMatrix)
{
	CtrlComponent* tempCtrlComp;

	// 存储控制元件状态
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->storeInternalVariables_pre();
	}

	// 存储当前时步结果
	//// 针对case271,暂时假设只有1个PWM变流器，且已知对应控制系统的节点编号
	//for (int j=0; j<4; j++) {
	//	PWMPredictMatrix[0][j] = ctrlNodeValue[j];
	//}
	// 针对case272
	for (int j=0; j<4; j++) {
		PWMPredictMatrix[0][j] = ctrlNodeValue[j+64];
	}

	// 求解控制系统在一个开关周期中的解
	for (int k=1; k<=Nc; k++) {
		int branchCaled=0;
		transferMeasurands();
		// 控制系统求解
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
						tempCtrlComp->saveInNodeValue(ctrlNodeValue); // 读取输入节点的值
						tempCtrlComp->calculateOutputValue(curTime+k*deltaT); // 求解各模块
						tempCtrlComp->saveOutNodeValue(ctrlNodeValue); // 将求解结果写入ctrlNodeValue向量中
						tempCtrlComp->markOutputNode(nodeCalMark); // 设置标志位，表示输出信号值已确定
						branchCalMark[i] = 1; // 设置标志位，表示该模块完成计算
					}
				}		
			}
		}

		// 清除标志位
		for (int i=0;i<nCtrlBranches;i++)
		{
			branchCalMark[i]=0;
		}
		for (int i=0;i<nCtrlNodes;i++)
		{
			nodeCalMark[i]=0;
		}

		// 存储求解结果
		//// 针对case271,暂时假设只有1个PWM变流器，且已知对应控制系统的节点编号
		//for (int j=0; j<4; j++) {
		//	PWMPredictMatrix[k][j] = ctrlNodeValue[j];
		//}
		// 针对case272
		for (int j=0; j<4; j++) {
			PWMPredictMatrix[k][j] = ctrlNodeValue[j+64];
		}
	}

	// 恢复控制元件状态
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->restoreInternalVariables_pre();
	}
}

// 校正开关动作时刻
// ver 0.1，暂时认为所有PWM变流器的载波频率相同，且载波均在t=0时刻取最大值
void EMTP::correctSwitchingInstants(double tol)
{
	int Nc = (len_subinterval+deltaT/2)/deltaT; // 计算一个迭代校正子区间包含几个时步数

	if ( counter % Nc == 1 ) // 如果是子区间结束时刻，进行校正
	{
		//// 创建临时数组用于存储校正时所需的调制波瞬时值和载波瞬时值
		//double** PWMCorrectMatrix = new double* [Nc+1];
		//for (int k=0; k<=Nc; k++) {
		//	PWMCorrectMatrix[k] = new double [4*nPWMConverter];
		//}

		//// 从控制系统求解结果中获取与PWM变流器相关的部分
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

		// 计算各个PWM变流器的开关时刻
		int ptr = 0;
		int changeOrNot = 0;
		for (int k=0; k<nPWMConverter; k++) {
			Component* tempCom = branches->at(branchNumArray_PWMConverter[k]);
			//if ( tempCom->correctSwithcingInstants(PWMCorrectMatrix, tol, Nc, ptr) == 1 )
			if ( tempCom->correctSwithcingInstants(ctrlResultMatrix_2, tol, Nc, ptr) == 1 )
				changeOrNot = 1;
		}

		// 校正
		if ( changeOrNot == 1 ) // 若需要校正，则
		{
			// 重新形成节点导纳矩阵
			formConductanceMatrix();

			// 仿真时间回到子区间的起始时刻
			counter = counter - Nc;
			curTime = curTime - Nc * deltaT;
			counter2 = 0;

			// 恢复电气元件状态
			for (int k=0; k<nBranch; k++) {
				Component* tempCom = branches->at(k);
				tempCom->restoreInternalVariables();
			}

			// 恢复控制元件状态
			for (int k=0; k<nCtrlBranches; k++) {
				CtrlComponent* tempCtrl = ctrlBranches->at(k);
				tempCtrl->restoreInternalVariables();
			}
		}
		else // 若无需校正，则
		{
			// 存储当前时刻电气元件状态
			for (int k=0; k<nBranch; k++) {
				Component* tempCom = branches->at(k);
				tempCom->storeInternalVariables();
			}

			// 存储控制元件状态
			for (int k=0; k<nCtrlBranches; k++) {
				CtrlComponent* tempCtrl = ctrlBranches->at(k);
				tempCtrl->storeInternalVariables();
			}

			for (int k=0; k<4*nPWMConverter;k++)
			{
				ctrlResultMatrix_2[0][k]=ctrlResultMatrix_2[counter2][k];
			}
			counter2 = 0;

			// 启动下一开关周期的预测
			// predictSwitchingInstants();

			// debug
			// cout << "\r" << curTime - len_subinterval << endl;
		}
	}
}

// ===============================================================
// xuyin, 20130224
// ===============================================================

// 求解系统在一个分段子区间上的解，即计算tau_new=G(tau)
void EMTP::solveG(double * tau_new)
{
	nSolveG++;
	// 求解一个分段子区间上的解
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

		transferControlVariables(); //将控制信号值传输至电气系统		
		calculateBranchNortonEquivalentCurrent(curTime); //计算支路诺顿等效电流
		formNodeNortonEquivalentCurrentArray(); //形成节点诺顿等效电流向量
		solveNodeVoltageEquation(); //求解节点电压方程
		saveNodeVoltage(); //保存节点电压向量,将求解结果存在全局电压矩阵中，与各支路无关。
		calculateBranchVoltage(); //计算支路电压并保存
		calculateBranchCurrent(); //计算支路电流
		saveBranchCurrent(); //保存支路电流			
		saveMahineWr(); //保存电机转速
		transferMeasurands(); //将电气系统测量值传输到控制系统
		solveCtrlSystem(); //求解控制系统
	}

	// 计算开关动作时刻
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

// 判断是否需要迭代
int EMTP::check_tau(double * tau_new, double tol)
{
	// 计算最大误差
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

// 预测
void EMTP::predict_tau()
{
	//// 就近预测法
	//for (int k=0; k<3*nPWMConverter; k++)
	//{
	//	tau_s[cnt_p][k] = tau_s[cnt_p-1][k];
	//}

	// 线性外推法
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

	//限制预测的占空比在0到1之间
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

	// 更新PWM变流器的诺顿等值导纳矩阵
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

	// 更新节点导纳矩阵
	formConductanceMatrix();
}

// 校正
void EMTP::correct_tau(double * tau_new)
{
	// 线性迭代法
	for (int k=0; k<3*nPWMConverter; k++)
	{
		tau_s[cnt_p][k] = tau_s[cnt_p][k] + 0.6 * (tau_new[k] - tau_s[cnt_p][k]);
	}

	// 更新PWM变流器的诺顿等值导纳矩阵
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

	// 更新节点导纳矩阵
	formConductanceMatrix();
}

// 保存状态
void EMTP::store_state()
{
	// 存储当前时刻电气元件状态
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->storeInternalVariables();
	}

	// 存储控制元件状态
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

// 恢复状态
void EMTP::restore_state()
{
	// 仿真时间回到子区间的起始时刻
	counter = counter - Ns;
	curTime = curTime - Ns * deltaT;
	counter2 = 0;

	// 恢复电气元件状态
	for (int k=0; k<nBranch; k++) {
		Component* tempCom = branches->at(k);
		tempCom->restoreInternalVariables();
	}

	// 恢复控制元件状态
	for (int k=0; k<nCtrlBranches; k++) {
		CtrlComponent* tempCtrl = ctrlBranches->at(k);
		tempCtrl->restoreInternalVariables();
	}
}
// ===============================================================