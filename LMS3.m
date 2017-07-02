%% MVDR criterion   LMS ALGORITHM
%% 初始化参数
close all;clear all;clc;
source=1;           % 信源  signal number
interference=1;     % 干扰  interference number
N=7;                % array number     阵元数
theta_s=20;         % DOA of signal
theta_i=-60;        % DOA of interference
ss=1024;            % snapshot  快拍数
snr=[-10 40];    %  SNR  信噪比 20 30 50
j=sqrt(-1);
%% 信号复包络
w=[pi/6 pi/4]';
for m=1:(source+interference)
    S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*(0:ss-1));        %3*1024
end
% for m=1:(source+interference)
%     S(m,:) = 10.^(snr(m)/10)*(randn(1,ss)+j*randn(1,ss));         %Signal and interference
% end
%% 阵列流形
A_i=exp(-j*pi*(0:N-1)'*sin(theta_i/180*pi));%8*4
A_s=exp(-j*pi*(0:N-1)'*sin(theta_s/180*pi));%8*1信号方向40
A = [A_s A_i(:,1:interference)];
%% 噪声
n=randn(N,ss)+j*randn(N,ss);
%% 观测信号
X=A*S+n;
%% LMS
Wx = A_s';                          % 静态权矢量
u=2^(-35);                          % 步长因子
% B0H_B0=[0 0 0 0 0 0 0;...           % 阻塞矩阵
%         0 1 0 0 0 0 0;...
%         0 0 1 0 0 0 0 ;...
%         0 0 0 1 0 0 0;...
%         0 0 0 0 1 0 0;...
%         0 0 0 0 0 1 0;...
%         0 0 0 0 0 0 1];
B0=[-1 exp(-j*pi*sin(theta_s)) 0 0 0 0 0;...        % 阻塞矩阵
    0 -1 exp(-j*pi*sin(theta_s)) 0 0 0 0;...
    0 0 -1 exp(-j*pi*sin(theta_s)) 0 0 0;...
    0 0 0 -1 exp(-j*pi*sin(theta_s)) 0 0;...
    0 0 0 0 -1 exp(-j*pi*sin(theta_s)) 0;...
    0 0 0 0 0 -1 exp(-j*pi*sin(theta_s))];
B0H_B0 = B0'*B0;
dataout=zeros(1,ss);               % 输出数据初始化  
dataout(1,1)=Wx*X(:,1);            % 输出数据初值
for i=1:length(dataout)-1           % 迭代 求权值
      Wx=Wx-u*(X(:,i)')*B0H_B0*dataout(1,i);        % 步长因子的影响很大
      dataout(1,i+1)=Wx*X(:,i+1);
end
%% 阵列方向图
phi=-89:1:90;
a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
F=Wx*a;
G=abs(F).^2./max(abs(F).^2);
G_dB=10*log10(G);
figure();
plot(phi,G_dB,'linewidth',2);legend('N=7,d=lamda/2');
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
grid on;
