function [ A, theta ] = cal_mag_theta( xabc )
%% ���������ڼ�������Գ��źŵķ�ֵ�����
% xabc: ��������Գ��źţ�abc�����µ�˲ʱֵ��������
% A: �����źŵķ�ֵ
% theta: �����źŵ����

%% ���źű任��alpha_beta����ϵ��
% alpha_beta�任����
Tabc2ab = [
    2/3 -1/3 -1/3
    0 1/sqrt(3) -1/sqrt(3)
    ];

% ab�任�� abc ==> alpha_beta
xab = Tabc2ab * xabc;

%% �����ֵ�����
% ��ֵ
A = sqrt( xab(1)^2 + xab(2)^2 );

% ���
if abs(xab(1)) < 1e-10 && abs(xab(2)) < 1e-10
    theta = 0;
elseif abs(xab(1)) < 1e-10
    if xab(2) > 0
        theta = pi/2;
    else
        theta = -pi/2;
    end
elseif abs(xab(2)) < 1e-10
    if xab(1) > 0
        theta = 0;
    else
        theta = pi;
    end
else
    theta = atan( xab(2) / xab(1) );
    if xab(1) < 0
        if theta < 0
            theta = theta + pi;
        else
            theta = theta - pi;
        end
    end
end

end

