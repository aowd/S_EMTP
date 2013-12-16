% =============================================================
%         以case 4171为基础，生成一系列风机算例
%                 它们都与同一个电源相连
% =============================================================
clear all
clc

% 读入初始文件
load ./caseDfn/case4171.txt
load ./caseDfn/case4171ctrl.txt
load ./caseDfn/case4171msr.txt

% 定义节点数
nElecNode=3;
nElecBranch=3;
nDFIGNode = 14;
nDFIGBranch = 12;
nCtrlNode = 193;

% 定义DFIG数量
nDFIG = 33;

%% 生成电气系统定义矩阵
elecMatrix = case4171;
temp = case4171(nElecBranch+1:end,:);
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,2) ~= 0 && temp(i,2) ~= 1 && temp(i,2) ~= 2 && temp(i,2) ~= 3
            temp(i,2) = temp(i,2) + nDFIGNode;
        end
        if temp(i,3) ~= 0 && temp(i,3) ~= 1 && temp(i,3) ~= 2 && temp(i,3) ~= 3
            temp(i,3) = temp(i,3) + nDFIGNode;
        end
        if temp(i,1)==24
            temp(i,7)=temp(i,7) + nCtrlNode;
        end
    end
    elecMatrix = [elecMatrix; temp];
end

%% 生成量测系统定义矩阵
msrMatrix = case4171msr;
temp = case4171msr;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nDFIGNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 3 || temp(i,1) == 4 ||temp(i,1) == 5 
            temp(i,2) = temp(i,2) + nDFIGBranch;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
    end
    msrMatrix = [msrMatrix; temp];
end

%% 生成控制系统定义矩阵
ctrlMatrix = case4171ctrl;
temp = case4171ctrl;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 3 || temp(i,1) == 5
            temp(i,2) = temp(i,2) + nCtrlNode;
        elseif temp(i,1) == 2 || temp(i,1) == 7 || temp(i,1) == 14
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        elseif temp(i,1) == 1
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
        elseif temp(i,1) == 12 || temp(i,1) == 13
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
        elseif temp(i,1) == 8 || temp(i,1) == 9 || temp(i,1) == 10 || temp(i,1) == 11
            temp(i,2) = temp(i,2) + nCtrlNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
            temp(i,4) = temp(i,4) + nCtrlNode;
            temp(i,5) = temp(i,5) + nCtrlNode;
            temp(i,6) = temp(i,6) + nCtrlNode;
        end
    end
    ctrlMatrix = [ctrlMatrix; temp];
end

%% 存储结果
stg = strcat(['./caseDfn/case417',num2str(nDFIG)],'.txt');
save(stg,'-ascii','elecMatrix');
stg = strcat(['./caseDfn/case417',num2str(nDFIG)],'ctrl.txt');
save(stg,'-ascii','ctrlMatrix');
stg = strcat(['./caseDfn/case417',num2str(nDFIG)],'msr.txt');
save(stg,'-ascii','msrMatrix');