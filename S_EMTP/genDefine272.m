clear all
clc

% 读入初始文件
load ./caseDfn/case272.txt
load ./caseDfn/case272ctrl.txt
load ./caseDfn/case272msr.txt

% 定义节点数
nElecNode = 5;
nCtrlNode = 85;
nElecBranch = size(case272,1);

% 定义PWM变流器数量
nPWMconverter = 180;

% 生成电气系统定义矩阵
elecMatrix = case272;
temp = case272;
for k = 2:nPWMconverter
    for i = 1:size(temp,1)
        if temp(i,2) ~= 0
            temp(i,2) = temp(i,2) + nElecNode;
        end
        if temp(i,3) ~= 0
            temp(i,3) = temp(i,3) + nElecNode;
        end
        if temp(i,1)==24
            temp(i,7)=temp(i,7) + nCtrlNode;
        end
    end
    elecMatrix = [elecMatrix; temp];
end

% 生成量测系统定义矩阵
msrMatrix = case272msr;
temp = case272msr;
for k = 2:nPWMconverter
    for i = 1:size(temp,1)
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nElecNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 3
            temp(i,2) = temp(i,2) + nElecBranch;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
    end
    msrMatrix = [msrMatrix; temp];
end

% 生成控制系统定义矩阵
ctrlMatrix = case272ctrl;
temp = case272ctrl;
for k = 2:nPWMconverter
    for i = 1:size(temp,1)
        if temp(i,1) == 3 || temp(i,1) == 5
            temp(i,2) = temp(i,2) + nCtrlNode;
        end
        if temp(i,1) == 2 || temp(i,1) == 7 || temp(i,1) == 14
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
        end
        if temp(i,1) == 12 || temp(i,1) == 13
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
        end
        if temp(i,1) == 8 || temp(i,1) == 9 || temp(i,1) == 10 || temp(i,1) == 11
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
            temp(i,6) = temp(i,6) + nCtrlNode;
        end
    end
    ctrlMatrix = [ctrlMatrix; temp];
end

stg = strcat(['./caseDfn/case272',num2str(nPWMconverter)],'.txt');
save(stg,'-ascii','elecMatrix');
stg = strcat(['./caseDfn/case272',num2str(nPWMconverter)],'ctrl.txt');
save(stg,'-ascii','ctrlMatrix');
stg = strcat(['./caseDfn/case272',num2str(nPWMconverter)],'msr.txt');
save(stg,'-ascii','msrMatrix');
% save case2722.txt elecMatrix -Ascii
% save case2722ctrl.txt ctrlMatrix -Ascii
% save case2722msr.txt msrMatrix -Ascii
% save case2725.txt elecMatrix -Ascii
% save case2725ctrl.txt ctrlMatrix -Ascii
% save case2725msr.txt msrMatrix -Ascii
% save case27210.txt elecMatrix -Ascii
% save case27210ctrl.txt ctrlMatrix -Ascii
% save case27210msr.txt msrMatrix -Ascii