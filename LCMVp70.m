%% MNV criterion  MVDR beamforming
%% 初始化参数 initial parameter
close all;clear all;clc;
source=1;           %信源  signal number
interference=2;     %干扰  interference number
N=16;               %array number     阵元数
theta_s=0;          %DOA of signal
theta_i=[-60 40];  %DOA of interference
ss=1024;            %snapshot  快拍数
snr=[-10 40 20];    %  SNR  信噪比
j=sqrt(-1);
%% 信号复包络 SIGNAL
% w=[pi/6 pi/6 pi/3 pi/4]';
% for m=1:4
%     S(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*[0:ss-1]);        %3*1024
% end

for m=1:(source+interference)
    S(m,:) = 10.^(snr(m)/10)*(randn(1,ss)+j*randn(1,ss));         %Signal and interference
end
%% 阵列流形  STEERING VECTOR

A_i=exp(-j*(0:N-1)'*pi*sin(theta_i/180*pi));%8*4
A_s=exp(-j*pi*(0:N-1)'*sin(theta_s*pi/180));%8*1信号方向40
A = [A_s A_i(:,1:interference)];
%% 噪声  NOISE
n=randn(N,ss)+j*randn(N,ss);

%% 观测信号  SIGNAL RECEIVED
X=A*S+n;

%% 阵列协方差矩阵  COVIARIANCE MATRIX
R=X*X'/ss;
Inv_Rx = inv(R);
As = [A_s A_i(:,1:interference)];     %****
f = [1 0 0]';       %****
  W_opt=Inv_Rx*As*inv(As'*Inv_Rx*As)*f;%计算权值   
  W_opt=W_opt/sqrt(W_opt'*W_opt);%归一化

%% 阵列方向图  pattern

phi=-89:1:90;
a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
F=W_opt'*a;
%figure();
%plot(phi,F);
G=abs(F).^2./max(abs(F).^2);
% G=abs(F).^2;
G_dB=10*log10(G);
%figure();
%plot(phi,G);legend('N=8,d=lamda/2');
figure();
plot(phi,G_dB,'linewidth',2);legend('N=16,d=lamda/2');
xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
grid on;
%axis([-90 90 -50 0]);