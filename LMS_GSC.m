%% MMSE criterion using LMS algorithm
%% 初始化参数  initial parameter
close all;clear all;clc;
N=7;                       % sensor阵元数
theta=[10 -60];      % 到达角 10为信号 其他为干扰 -20 40
theta_s = 10;
theta_i = -60;
% direction of arrival. the first one is expect signal; others are interference signals
ss=4400;                    % snapshot  快拍数  32
snr=[10 30];                     % SNR 信噪比 40 30
j=sqrt(-1);
%% 信号复包络  signal
w=[pi/5 pi/6]';             % pi/4 pi/3
for m = 1:length(theta)
    S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*(0:ss-1));   
end
% d=10.^(snr(1)/10)*exp(-j*w(1)*(0:ss-1));   %期望信号  pilot signal
% for m=1:length(theta)
%     S(m,:) = 10.^(snr(m)/10)*(randn(1,ss)+j*randn(1,ss));         %Signal and interference
% end
%% 阵列流形  steering vector
A=exp(-j*(0:N-1)'*pi*sin(theta/180*pi));                %8*4

%% 噪声   noise
n=randn(N,ss)+j*randn(N,ss);

%% 观测信号   received signal
X=A*S+n;
d=A(:,1)'*X;
%% LMS
% B0=[-1 exp(-j*pi*sin(theta_s)) 0 0 0 0 0;...        % 阻塞矩阵
%     0 -1 exp(-j*pi*sin(theta_s)) 0 0 0 0;...
%     0 0 -1 exp(-j*pi*sin(theta_s)) 0 0 0;...
%     0 0 0 -1 exp(-j*pi*sin(theta_s)) 0 0;...
%     0 0 0 0 -1 exp(-j*pi*sin(theta_s)) 0;...
%     0 0 0 0 0 -1 exp(-j*pi*sin(theta_s));...
%     1 1 1 1 1 1 1];
% B0H_B0 = B0'*B0;
% B0H_B0=[0 0 0 0 0 0 0;...           % 阻塞矩阵
%         0 1 0 0 0 0 0;...
%         0 0 1 0 0 0 0 ;...
%         0 0 0 1 0 0 0;...
%         0 0 0 0 1 0 0;...
%         0 0 0 0 0 1 0;...
%         0 0 0 0 0 0 1];

Wx = zeros(N,1);                       % ;A(:,1)
dataout=zeros(1,ss);               % 输出数据初始化  
dataout(1,1)=d(1,1);            % 输出数据初值
u=2^(-45);                          % 步长因子
for i=1:length(dataout)-1           % 迭代 求权值
      Wx=Wx - u*X(:,i)*dataout(1,i);        % 步长因子的影响很大  conj
      dataout(1,i+1) = d(1,i+1) - Wx'*X(:,i+1);
end

%% 阵列方向图

phi=-89:1:90;
a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
F=Wx'*a;
G=abs(F).^2./max(abs(F).^2);
G_dB=10*log10(G);
figure();
plot(phi,G_dB,'linewidth',2);hold on;grid on;
% legend('SMI','DL-SMI');                        % *******
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
%axis([-90 90 -50 0]);legend('N=16,d=lamda/2');