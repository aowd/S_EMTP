clear all;
clc;

%% 读取文件
% 算例信息
iCase = 295;
nFile = 3;
deltaT =5;

% 读取PSCAD数据文件
result_PSC = [];
for k = 1:nFile
    if k<9.5
        stg = sprintf('../result/pscadData/case%d_0%d.out',iCase,k);
    else
        stg = sprintf('../result/pscadData/case%d_%d.out',iCase,k);
    end
    result_temp = load(stg);
    time_PSC = result_temp(:,1);
    result_PSC = [result_PSC result_temp(:,2:end)*1000];
end

% 读取CPP数据文件
stg = sprintf('../result/result_latest/result%d.dat',iCase);
result_temp=load(stg);
t=result_temp(:,1);%时间序列
result_CPP=result_temp(:,2:end);%CPP数据

% 数据处理
idx = 1;
t_tmp = time_PSC(idx);
while abs(t_tmp - t(1))>1e-6
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
startIdx = idx;
while abs(t_tmp - t(end))>1e-6
    idx = idx + 1;
    t_tmp = time_PSC(idx);
end
endIdx = idx;

if deltaT >= 5
    steps = deltaT/5;
    result_PSC = result_PSC(startIdx:steps:endIdx,:);
else
    t = time_PSC(startIdx:endIdx);
    steps = 5/deltaT;
    result_CPP = result_CPP(1:steps:end,:);
    result_PSC = result_PSC(startIdx:endIdx,:);
end

%% 选择需要观察的波形
Iabc_PSC = result_PSC(:,16:18);
Iabc_CPP = result_CPP(:,16:18);

%% 滤波
% 滤波器分子系数
load Num_ac_5us.mat
Num_ac = Num_ac_5us;

% testStart = size(Num_ac,2)+1;
testStart = int32(size(t,1)*0.8);

% 三相电感电流abc瞬时值
Iabc_PSC_filter = filter(Num_ac,1,Iabc_PSC);
Iabc_CPP_filter = filter(Num_ac,1,Iabc_CPP);

% 幅值相角计算
for k = 1:size(Iabc_PSC_filter,1)
    [Mag_PSC_filter(k) phase_PSC_filter(k)] = cal_mag_theta(Iabc_PSC_filter(k,:)');
    [Mag_CPP_filter(k) phase_CPP_filter(k)] = cal_mag_theta(Iabc_CPP_filter(k,:)');
end

% delta_Mag = abs( ( Mag_CPP_filter(10000:end) - Mag_PSC_filter(10000:end) ) ./ Mag_PSC_filter(10000:end) );
delta_Mag = abs( ( Mag_CPP_filter(testStart:end) - Mag_PSC_filter(testStart:end) ) ./ Mag_PSC_filter(testStart:end) );
Err_Mag = max(delta_Mag)*100

tStart = int32(size(t,1)*0.1);
RE_mag = abs( Mag_CPP_filter(tStart:end) - Mag_PSC_filter(tStart:end) ) ./ abs(Mag_CPP_filter(tStart:end)) * 100;

% 相角差计算
for k = 1:size(Iabc_PSC_filter,1)
    delta_phase(k) = phase_CPP_filter(k) - phase_PSC_filter(k);
    if delta_phase(k) > pi
        delta_phase(k) = delta_phase(k) - 2*pi;
    end
    if delta_phase(k) < -pi
        delta_phase(k) = delta_phase(k) + 2*pi;
    end
    delta_phase(k) = delta_phase(k)/pi*180;
end

% Err_phase = max(abs(delta_phase(10000:end)));
Err_phase = max(abs(delta_phase(testStart:end)))

% 实部虚部计算
Ir_PSC_filter = (2/3)*Iabc_PSC_filter(:,1) - (1/3)*Iabc_PSC_filter(:,2) ...
     - (1/3)*Iabc_PSC_filter(:,3);
Ii_PSC_filter = (1/sqrt(3))*(Iabc_PSC_filter(:,2)-Iabc_PSC_filter(:,3));

Ir_CPP_filter = (2/3)*Iabc_CPP_filter(:,1) - (1/3)*Iabc_CPP_filter(:,2) ...
     - (1/3)*Iabc_CPP_filter(:,3);
Ii_CPP_filter = (1/sqrt(3))*(Iabc_CPP_filter(:,2)-Iabc_CPP_filter(:,3));

TVE = sqrt( ( (Ir_PSC_filter-Ir_CPP_filter).^2+(Ii_PSC_filter-Ii_CPP_filter).^2 ) ./ ( Ir_PSC_filter.^2 + Ii_PSC_filter.^2 ) );
% Err_iabc = max(TVE(10000:end))*100
Err_iabc = max(TVE(testStart:end))*100

%% 作图比较
figure(2)
plot(t,Iabc_PSC_filter,t,Iabc_CPP_filter);
% 
figure(3)
plot(t,Mag_PSC_filter,t,Mag_CPP_filter);

figure(4)
plot(t(tStart:end),RE_mag);
% 
% figure(4)
% plot(t,phase_PSC_filter,t,phase_CPP_filter);
% 
figure(5)
plot(t,delta_phase);
