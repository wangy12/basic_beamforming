%% MMSE criterion using SMI algorithm
%% 初始化参数  initial parameter
close all;clear all;clc;
N=16;                       % sensor阵元数
theta=[0 -60 -20 40];      % 到达角 10为信号 其他为干扰
% direction of arrival. the first one is expect signal; others are interference signals
ss=1024;                    % snapshot  快拍数  32
snr=[50 80 40 30];                     % SNR 信噪比
j=sqrt(-1);
%% 信号复包络  signal
w=[pi/5 pi/6 pi/4 pi/3]';
for m = 1:length(theta)
    S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*[0:ss-1]);   
end
d=10.^(snr(1)/10)*exp(-j*w(1)*[0:ss-1]);   %期望信号  pilot signal
%% 阵列流形  steering vector
A=exp(j*(0:N-1)'*pi*sin(theta/180*pi));                %8*4
%% 噪声   noise
n=randn(N,ss)+j*randn(N,ss);
%% 观测信号   received signal
X=A*S+n;
%% 阵列协方差矩阵   covariance matrix of pilot signal and received signal
R=X*X'/ss;
%% 期望信号与接收信号复包络协方差矩阵  covariance matrix
rxd=X*d'/ss;
%% 阵列方向图
W=inv(R)*rxd;           % weighting vector
phi=-89:1:90;
a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));
F=W'*a;
G=abs(F).^2./max(abs(F).^2);
G_dB=10*log10(G);
figure();
plot(phi,G_dB,'linewidth',2);hold on;grid on;     
legend('MMSE');                       
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
