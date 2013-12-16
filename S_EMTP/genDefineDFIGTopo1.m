clear all
clc

% 读入初始文件
load('./caseDfn/WF_CTRL.txt');
load('./caseDfn/WF_Topo1_MSR.txt');
load('./caseDfn/WF_Topo1_CS_Part1.txt');
load('./caseDfn/WF_Topo1_CS_Part2.txt');

% 定义节点数
nNode_CS_P1 =16;
nNode_CS_P2 = 24;
nBranch_CS_P1 = size(WF_Topo1_CS_Part1,1);
nBranch_CS_P2 = size(WF_Topo1_CS_Part2,1);
nCtrlNode = 193;

nDFIG = 1;


% 生成电气系统定义矩阵
elecMatrix = [WF_Topo1_CS_Part1;WF_Topo1_CS_Part2];
temp = WF_Topo1_CS_Part2;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,2) ~= 0
            temp(i,2) = temp(i,2) + nNode_CS_P2;
        end
        if temp(i,3) ~= 0 || temp(i,3) ~= 11 || temp(i,3) ~= 13 || temp(i,3) ~= 15
            temp(i,3) = temp(i,3) + nNode_CS_P2;
        end
        if temp(i,1)==24
            temp(i,7)=temp(i,7) + nCtrlNode;
        end
    end
    elecMatrix = [elecMatrix; temp];
end

% 生成量测系统定义矩阵
temp = WF_Topo1_MSR;
for i = 1:size(temp,1)
    if temp(i,1) == 3 || temp(i,1) == 4 || temp(i,1)==5
        temp(i,2) = temp(i,2) + nBranch_CS_P1;
    end
end
msrMatrix = temp;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nNode_CS_P2;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 3 || temp(i,1) == 4 || temp(i,1)==5
            temp(i,2) = temp(i,2) + nBranch_CS_P2;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
    end
    msrMatrix = [msrMatrix; temp];
end

% 生成控制系统定义矩阵
ctrlMatrix = WF_CTRL;
temp = WF_CTRL;
for k = 2:nDFIG
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

save ./caseDfn/case413.txt elecMatrix -Ascii
save ./caseDfn/case413ctrl.txt ctrlMatrix -Ascii
save ./caseDfn/case413msr.txt msrMatrix -Ascii