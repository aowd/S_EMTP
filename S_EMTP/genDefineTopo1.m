clear all
clc

% 读入初始文件
load('./caseDfn/case407.txt')
case409ctrl=load('./caseDfn/case407ctrl.txt');
case409msr=load('./caseDfn/case407msr.txt');

% 定义节点数
nElecNode =66;
nDFIGNode = 14;
nElecBranch = 48;
nDFIGBranch = 12;
nCtrlNode = 193;
nElecBranch_1=51;

nDFIG = 5;


% 生成量测系统定义矩阵
temp = case409msr;
for i = 1:size(temp,1)
    if temp(i,1) == 1
        temp(i,2) = temp(i,2) + nElecNode;
    end
    if temp(i,1) == 3
        temp(i,2) = temp(i,2) + nElecBranch;
    end
    if temp(i,1) == 4 || temp(i,1)==5
        temp(i,2) = temp(i,2) + nElecBranch_1;
    end
end
msrMatrix = temp;
for k = 2:nDFIG
    for i = 1:size(temp,1)
        if temp(i,1) == 1
            temp(i,2) = temp(i,2) + nDFIGNode;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
        if temp(i,1) == 3 || temp(i,1) == 4 || temp(i,1)==5
            temp(i,2) = temp(i,2) + nDFIGBranch;
            temp(i,3) = temp(i,3) + nCtrlNode;
        end
    end
    msrMatrix = [msrMatrix; temp];
end

% 生成控制系统定义矩阵
ctrlMatrix = case409ctrl;
temp = case409ctrl;
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

save ./caseDfn/case409ctrl.txt ctrlMatrix -Ascii
save ./caseDfn/case409msr.txt msrMatrix -Ascii