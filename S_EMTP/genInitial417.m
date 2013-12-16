% =========================================================================
%                    由4171的结果文件生成多风机算例的初始化文件
%                          风机共同接到同一个电源上
% =========================================================================
clear all
clc

% 读入文件
load ./result/pscadData/case4171_I.dat
load ./result/pscadData/case4171_V.dat
load ./result/pscadData/case4171_Gen_01.out
load ./result/pscadData/case4171_Gen_02.out
case4171_Gen = [case4171_Gen_01(:,2:end) case4171_Gen_02(:,2:end)];
load ./result/pscadData/case4171_PI_01.out
load ./result/pscadData/case4171_PI_02.out
case4171_PI = [case4171_PI_01(:,2:end) case4171_PI_02(:,2:end)];
load ./result/pscadData/case4171_PWM_01.out
case4171_PWM = case4171_PWM_01(:,2:end);

% 定义DFIG数量
nDFIG = 18;

nElecNode = 3;
nElecBranch = 3;
nDFIGNode = 14;
nDFIGBranch = 25;

%% 生成电气系统初始化矩阵
VMatrix = case4171_V;
IMatrix = case4171_I;
for k = 2:nDFIG
    VMatrix = [VMatrix case4171_V(:,nElecNode+1+1:end)];
    IMatrix = [IMatrix case4171_I(:,nElecBranch+1+1:end)];
end
Matrix=[VMatrix IMatrix];

stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_V.dat']);
save(stg,'-ascii','VMatrix');
stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_I.dat']);
save(stg,'-ascii','IMatrix');

%% 生成电机初始化矩阵
t_Gen = case4171_Gen_01(:,1);
GenMatrix = [];
for k = 1:nDFIG
    GenMatrix = [GenMatrix case4171_Gen];
end

counter = 1;
stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_Gen_',num2str(counter,'%02d')],'.out');
while size(GenMatrix,2)>10
    temp = [t_Gen GenMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    GenMatrix = GenMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_Gen_',num2str(counter,'%02d')],'.out');
end
temp = [t_Gen GenMatrix];
save(stg,'-ascii','temp');

%% 生成PI初始化矩阵
t_PI = case4171_PI_01(:,1);
PIMatrix = [];
for k = 1:nDFIG
    PIMatrix = [PIMatrix case4171_PI];
end

counter = 1;
stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_PI_',num2str(counter,'%02d')],'.out');
while size(PIMatrix,2)>10
    temp = [t_PI PIMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PIMatrix = PIMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_PI_',num2str(counter,'%02d')],'.out');
end
temp = [t_PI PIMatrix];
save(stg,'-ascii','temp');

%% 生成PWM初始化矩阵
t_PWM = case4171_PWM_01(:,1);
PWMMatrix = [];
for k = 1:nDFIG
    PWMMatrix = [PWMMatrix case4171_PWM];
end

counter = 1;
stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_PWM_',num2str(counter,'%02d')],'.out');
while size(PWMMatrix,2)>10
    temp = [t_PWM PWMMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PWMMatrix = PWMMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case417',num2str(nDFIG),'_PWM_',num2str(counter,'%02d')],'.out');
end
temp = [t_PWM PWMMatrix];
save(stg,'-ascii','temp');