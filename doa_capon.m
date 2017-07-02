%% DOA CAPON
%% 初始化参数 initial parameter
close all;
clear all;clc;
source=2;       %信源  signal number 期望信号
N=4;            %array number
theta_s=[20 60];      %DOA of signal
ss=1024;        %snapshot  快拍数
snr=[1 1];    %  SNR  信噪比
j=sqrt(-1);
%% 信号复包络 SIGNAL
w=[pi/6 pi/5]';
for m=1:source
    S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*[0:ss-1]);
end
%% 阵列流形  STEERING VECTOR
A=exp(-j*(0:N-1)'*pi*sin(theta_s/180*pi));
%% 噪声  NOISE
n=randn(N,ss)+j*randn(N,ss);
%% 观测信号  SIGNAL RECEIVED
X=A*S+n;
%% 阵列协方差矩阵  COVIARIANCE MATRIX
R=X*X'/ss;
inv_R = inv(R);
%% Capon DOA
for phi=1:1:90;
    a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
    Pcapon(phi)=1/(a'*inv_R*a);
end
phi_scan=1:1:90;
Pcapon=10*log10(Pcapon);
figure();
plot(phi_scan,Pcapon,'linewidth',2);legend('N=7,d=lamda/2');
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
grid on;

