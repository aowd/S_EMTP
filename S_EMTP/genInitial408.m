% =========================================================================
%                    由408的结果文件生成多风机单独算例的初始化文件
% =========================================================================
clear all
clc

% 读入文件
load ./result/pscadData/case408_01_s.out
load ./result/pscadData/case408_02_s.out
load ./result/pscadData/case408_03_s.out
load ./result/pscadData/case408_04_s.out
load ./result/pscadData/case408_05_s.out
case408 = [case408_01_s(:,2:end) case408_02_s(:,2:end) case408_03_s(:,2:end) case408_04_s(:,2:end) case408_05_s(:,2:end)];
load ./result/pscadData/case408_Gen_01_s.out
load ./result/pscadData/case408_Gen_02_s.out
case408_Gen = [case408_Gen_01_s(:,2:end) case408_Gen_02_s(:,2:end)];
load ./result/pscadData/case408_PI_01_s.out
load ./result/pscadData/case408_PI_02_s.out
case408_PI = [case408_PI_01_s(:,2:end) case408_PI_02_s(:,2:end)];
load ./result/pscadData/case408_PWM_01_s.out
case408_PWM = case408_PWM_01_s(:,2:end);

% 定义DFIG数量
nDFIG = 72;

nElecNode = 17;

%% 生成电气系统初始化矩阵
t = case408_01_s(:,1);
VMatrix = [];
IMatrix = [];
for k = 1:nDFIG
    VMatrix = [VMatrix case408(:,1:nElecNode)];
    IMatrix = [IMatrix case408(:,nElecNode+1:end)];
end
Matrix=[VMatrix IMatrix];

counter = 1;
stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_',num2str(counter,'%02d')],'.out');
while size(Matrix,2)>10
    temp = [t Matrix(:,1:10)];
    save(stg,'-ascii','temp');
    Matrix = Matrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_',num2str(counter,'%02d')],'.out');
end
temp = [t Matrix];
save(stg,'-ascii','temp');

%% 生成电机初始化矩阵
t_Gen = case408_Gen_01_s(:,1);
GenMatrix = [];
for k = 1:nDFIG
    GenMatrix = [GenMatrix case408_Gen];
end

counter = 1;
stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_Gen_',num2str(counter,'%02d')],'.out');
while size(GenMatrix,2)>10
    temp = [t_Gen GenMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    GenMatrix = GenMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_Gen_',num2str(counter,'%02d')],'.out');
end
temp = [t_Gen GenMatrix];
save(stg,'-ascii','temp');

%% 生成PI初始化矩阵
t_PI = case408_PI_01_s(:,1);
PIMatrix = [];
for k = 1:nDFIG
    PIMatrix = [PIMatrix case408_PI];
end

counter = 1;
stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_PI_',num2str(counter,'%02d')],'.out');
while size(PIMatrix,2)>10
    temp = [t_PI PIMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PIMatrix = PIMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_PI_',num2str(counter,'%02d')],'.out');
end
temp = [t_PI PIMatrix];
save(stg,'-ascii','temp');

%% 生成PWM初始化矩阵
t_PWM = case408_PWM_01_s(:,1);
PWMMatrix = [];
for k = 1:nDFIG
    PWMMatrix = [PWMMatrix case408_PWM];
end

counter = 1;
stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_PWM_',num2str(counter,'%02d')],'.out');
while size(PWMMatrix,2)>10
    temp = [t_PWM PWMMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PWMMatrix = PWMMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case408',num2str(nDFIG),'_PWM_',num2str(counter,'%02d')],'.out');
end
temp = [t_PWM PWMMatrix];
save(stg,'-ascii','temp');
