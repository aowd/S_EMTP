clear all
clc

nElecNode=5;
% 读入文件
load ./result/pscadData/case272_01.out
load ./result/pscadData/case272_02.out
case272 = [case272_01(:,2:end) case272_02(:,2:end)];
load ./result/pscadData/case272_PI_01.out
load ./result/pscadData/case272_PWM_01.out

% 定义PWM变流器数量
nPWMconverter = 180;

% 生成电气系统初始化矩阵
t = case272_01(:,1);
VMatrix = [];
IMatrix = [];
for k = 1:nPWMconverter
    VMatrix = [VMatrix case272(:,1:nElecNode)];
    IMatrix = [IMatrix case272(:,nElecNode+1:end)];
end
Matrix=[VMatrix IMatrix];

counter = 1;
stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case2725_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case27210_',num2str(counter,'%02d')],'.out');
while size(Matrix,2)>10
    temp = [t Matrix(:,1:10)];
    save(stg,'-ascii','temp');
    Matrix = Matrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case2725_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case27210_',num2str(counter,'%02d')],'.out');
end
temp = [t Matrix];
save(stg,'-ascii','temp');

% 生成PI初始化矩阵
t_PI = case272_PI_01(:,1);
PIMatrix = [];
for k = 1:nPWMconverter
    PIMatrix = [PIMatrix case272_PI_01(:,2:end)];
end

counter = 1;
stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_PI_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case2725_PI_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case27210_PI_',num2str(counter,'%02d')],'.out');
while size(PIMatrix,2)>10
    temp = [t_PI PIMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PIMatrix = PIMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_PI_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case2725_PI_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case27210_PI_',num2str(counter,'%02d')],'.out');
end
temp = [t_PI PIMatrix];
save(stg,'-ascii','temp');

% 生成PWM初始化矩阵
t_PWM = case272_PWM_01(:,1);
PWMMatrix = [];
for k = 1:nPWMconverter
    PWMMatrix = [PWMMatrix case272_PWM_01(:,2:end)];
end

counter = 1;
stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_PWM_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case2725_PWM_',num2str(counter,'%02d')],'.out');
% stg = strcat(['case27210_PWM_',num2str(counter,'%02d')],'.out');
while size(PWMMatrix,2)>10
    temp = [t_PWM PWMMatrix(:,1:10)];
    save(stg,'-ascii','temp');
    PWMMatrix = PWMMatrix(:,11:end);
    counter = counter + 1;
    stg = strcat(['./result/pscadData/case272',num2str(nPWMconverter),'_PWM_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case2725_PWM_',num2str(counter,'%02d')],'.out');
%     stg = strcat(['case27210_PWM_',num2str(counter,'%02d')],'.out');
end
temp = [t_PWM PWMMatrix];
save(stg,'-ascii','temp');
