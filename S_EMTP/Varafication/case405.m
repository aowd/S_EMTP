clear all;
clc;

%% 读取文件
% 算例信息
iCase = 405;
nFile = 5;
deltaT = 5;

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

if deltaT >= 1
    steps = deltaT/1;
    result_PSC = result_PSC(startIdx:steps:endIdx,:);
else
    t = time_PSC(startIdx:endIdx);
    steps = 1/deltaT;
    result_CPP = result_CPP(1:steps:end,:);
    result_PSC = result_PSC(startIdx:endIdx,:);
end

% %% 选择需要观察的波形
% Udc_PSC = result_PSC(:,10)-result_PSC(:,11);
% Udc_CPP = result_CPP(:,10)-result_CPP(:,11);
% 
% % Iabc_PSC = result_PSC(:,15:17);
% % Iabc_CPP = result_CPP(:,15:17);
% 
% % Iabc_PSC = result_PSC(:,21:23);
% % Iabc_CPP = result_CPP(:,21:23);
% 
Iabc_PSC = result_PSC(:,35:37);
Iabc_CPP = result_CPP(:,35:37);
% 
% Iabc_PSC = result_PSC(:,30:32);
% Iabc_CPP = result_CPP(:,30:32);
% 
% %% 滤波
% % 滤波器分子系数
% load Num_dc.mat
load Num_ac_5us.mat
Num_ac = Num_ac_5us;
% 
% % testStart = size(Num_ac,2)+1;
testStart = int32(size(t,1)*0.8);
% 
% % 直流电压
% Udc_PSC_filter = filter(Num_dc,1,Udc_PSC);
% Udc_CPP_filter = filter(Num_dc,1,Udc_CPP);
% 
% % delta_Udc = abs( ( Udc_CPP_filter(10000:end) - Udc_PSC_filter(10000:end) ) ./ Udc_PSC_filter(10000:end) );
% delta_Udc = abs( ( Udc_CPP_filter(testStart:end) - Udc_PSC_filter(testStart:end) ) ./ Udc_PSC_filter(testStart:end) );
% Err_Udc = max(delta_Udc)*100
% 
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
% 
% %% 作图比较
% figure(1)
% plot(t,Udc_PSC_filter,t,Udc_CPP_filter);
% % 
figure(2)
plot(t,Iabc_PSC_filter,t,Iabc_CPP_filter);
% % 
figure(3)
plot(t,abs(Mag_CPP_filter-Mag_PSC_filter)./Mag_PSC_filter*100);
% % plot(t,Mag_PSC_filter,t,Mag_CPP_filter);
% % 
% % figure(4)
% % plot(t,phase_PSC_filter,t,phase_CPP_filter);
% % 
figure(5)
plot(t,delta_phase);
% 
% 
% Idc_PSC_filter = filter(Num_dc,1,result_PSC(:,38));
% Idc_CPP_filter = filter(Num_dc,1,result_CPP(:,38));
% figure(10)
% plot(t,Idc_PSC_filter,t,Idc_CPP_filter);
% 
% % figure(11)
% % plot(t,(Iabc_PSC_filter(:,1)+Iabc_PSC_filter(:,2)+Iabc_PSC_filter(:,3))/3,...
% %     t,(Iabc_CPP_filter(:,1)+Iabc_CPP_filter(:,2)+Iabc_CPP_filter(:,3))/3);
% 
% % Vr_PSC = result_PSC(:,4:6);
% % Vr_CPP = result_CPP(:,4:6);
% % Vr_PSC_filter = filter(Num_ac,1,Vr_PSC);
% % Vr_CPP_filter = filter(Num_ac,1,Vr_CPP);
% % figure(12)
% % plot(t,Vr_PSC_filter,t,Vr_CPP_filter)
% 
% Vr2_PSC = result_PSC(:,7:9);
% Vr2_CPP = result_CPP(:,7:9);
% Vr2_PSC_filter = filter(Num_ac,1,Vr2_PSC);
% Vr2_CPP_filter = filter(Num_ac,1,Vr2_CPP);
% figure(13)
% plot(t,Vr2_PSC_filter,t,Vr2_CPP_filter)
% 
% % Igsc_PSC = result_PSC(:,35:37);
% % Igsc_CPP = result_CPP(:,35:37);
% % Igsc_PSC_filter = filter(Num_ac,1,Igsc_PSC);
% % Igsc_CPP_filter = filter(Num_ac,1,Igsc_CPP);
% % figure(14)
% % plot(t,Igsc_PSC_filter,t,Igsc_CPP_filter)
% 
% % Pr_PSC = Vr2_PSC_filter(:,1).*Iabc_PSC_filter(:,1) + ...
% %     Vr2_PSC_filter(:,2).*Iabc_PSC_filter(:,2) + ...
% %     Vr2_PSC_filter(:,3).*Iabc_PSC_filter(:,3);
% % Pr_CPP = Vr2_CPP_filter(:,1).*Iabc_CPP_filter(:,1) + ...
% %     Vr2_CPP_filter(:,2).*Iabc_CPP_filter(:,2) + ...
% %     Vr2_CPP_filter(:,3).*Iabc_CPP_filter(:,3);
% % figure(15)
% % plot(t,Pr_PSC,t,Pr_CPP)
% 
% Pdc_PSC = Udc_PSC_filter.*Idc_PSC_filter;
% Pdc_CPP = Udc_CPP_filter.*Idc_CPP_filter;
% figure(16)
% plot(t,Pdc_PSC,t,Pdc_CPP)
% 
% % Ploss_CPP = 0.001*(Iabc_PSC_filter(:,1).^2 + ...
% %     Iabc_PSC_filter(:,2).^2 + Iabc_PSC_filter(:,3).^2);
% % figure(17)
% % plot(t,Ploss_CPP)
% 
% 
% Pr2_PSC = Vr2_PSC(:,1).*Iabc_PSC(:,1) + ...
%     Vr2_PSC(:,2).*Iabc_PSC(:,2) + ...
%     Vr2_PSC(:,3).*Iabc_PSC(:,3);
% Pr2_CPP = Vr2_CPP(:,1).*Iabc_CPP(:,1) + ...
%     Vr2_CPP(:,2).*Iabc_CPP(:,2) + ...
%     Vr2_CPP(:,3).*Iabc_CPP(:,3);
% Pr2_PSC_filter = filter(Num_dc,1,Pr2_PSC);
% Pr2_CPP_filter = filter(Num_dc,1,Pr2_CPP);
% figure(18)
% plot(t,Pr2_PSC_filter,t,Pr2_CPP_filter)

% t = t(int32(end/2):end);
% 
% % pscad结果
% vra_psc = result_PSC(int32(end/2):end,7);
% vrb_psc = result_PSC(int32(end/2):end,8);
% vrc_psc = result_PSC(int32(end/2):end,9);
% ira_psc = result_PSC(int32(end/2):end,30);
% irb_psc = result_PSC(int32(end/2):end,31);
% irc_psc = result_PSC(int32(end/2):end,32);
% udc1_psc = result_PSC(int32(end/2):end,10);
% udc2_psc = result_PSC(int32(end/2):end,11);
% idc1_psc = result_PSC(int32(end/2):end,38);
% idc2_psc = result_PSC(int32(end/2):end,39);
% 
% Pac_psc = vra_psc.*ira_psc + vrb_psc.*irb_psc + vrc_psc.*irc_psc;
% Pdc_psc = udc1_psc.*idc1_psc + udc2_psc.*idc2_psc;
% Ploss_psc = 0.001*(ira_psc.^2+irb_psc.^2+irc_psc.^2);
% 
% Pac_mean_psc = mean(Pac_psc)
% Pdc_mean_psc = mean(Pdc_psc)
% Ploss_mean_psc = mean(Ploss_psc)
% 
% 
% % CPP结果
% vra_cpp = result_CPP(int32(end/2):end,7);
% vrb_cpp = result_CPP(int32(end/2):end,8);
% vrc_cpp = result_CPP(int32(end/2):end,9);
% ira_cpp = result_CPP(int32(end/2):end,30);
% irb_cpp = result_CPP(int32(end/2):end,31);
% irc_cpp = result_CPP(int32(end/2):end,32);
% udc1_cpp = result_CPP(int32(end/2):end,10);
% udc2_cpp = result_CPP(int32(end/2):end,11);
% idc1_cpp = result_CPP(int32(end/2):end,38);
% idc2_cpp = result_CPP(int32(end/2):end,39);
% 
% Pac_cpp = vra_cpp.*ira_cpp + vrb_cpp.*irb_cpp + vrc_cpp.*irc_cpp;
% Pdc_cpp = udc1_cpp.*idc1_cpp + udc2_cpp.*idc2_cpp;
% Ploss_cpp = 0.001*(ira_cpp.^2+irb_cpp.^2+irc_cpp.^2);
% 
% Pac_mean_cpp = mean(Pac_cpp)
% Pdc_mean_cpp = mean(Pdc_cpp)
% Ploss_mean_cpp = mean(Ploss_cpp)
% 
% plot(t,Pac_psc,t,Pac_cpp)
% 
% idc1_mean_cpp = mean(idc1_cpp)
% idc2_mean_cpp = mean(idc2_cpp)
% 
% udc1_mean_cpp = mean(udc1_cpp)
% udc2_mean_cpp = mean(udc2_cpp)
