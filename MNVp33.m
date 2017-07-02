%% beamforming anti-jamming
%% 初始化参数 initial parameter
close all;
clear all;clc;
source=1;           %信源  signal number 期望信号
interference=3;     %干扰  interference number   6
N=7;               %array number     阵元数   7
theta_s=0;          %DOA of signal   -60  0
theta_i=[-60 -20 40 20 60 -40];  %DOA of interference
ss=2048;            %snapshot  快拍数
snr=[-10 40 20 30 50 30 30];    %  SNR  信噪比
j=sqrt(-1);
%% 信号复包络 SIGNAL
w=[pi/6 pi/6 pi/3 pi/4 pi/4 pi/4 pi/4]';
for m=1:(source+interference)
    SS(m,:)=10.^(snr(m)/10)*exp(-j*w(m)*[0:ss-1]);        %3*1024
    S(m,:) = awgn(SS(m,:),4,'measured');
end

% for m=1:(source+interference)
%     SS(m,:) = 10.^(snr(m)/10)*(randn(1,ss)+j*randn(1,ss));         %Signal and interference
%     S(m,:) = awgn(SS(m,:),10,'measured');
% end
%% 阵列流形  STEERING VECTOR
A_i=exp(-j*(0:N-1)'*pi*sin(theta_i/180*pi));%8*4
A_s=exp(-j*(0:N-1)'*pi*sin(theta_s/180*pi));%8*1信号方向40
A = [A_s A_i(:,1:interference)];
%% 噪声  NOISE
n=randn(N,ss)+j*randn(N,ss);
%% 观测信号  SIGNAL RECEIVED
X=A*S+n;
% X=A*S;
%% 阵列协方差矩阵  COVIARIANCE MATRIX
R=X*X'/ss;
inv_R = inv(R);
%% 最小噪声方差MNV   CAPON   MVDR  LVDS
P=inv(A_s'*inv_R*A_s);% 最小输出功率
W_mnv = P*inv_R*A_s;
% W_gui1 = W_mnv/W_mnv(1,1);
% W_guiyi = W_gui1/max(abs(real(W_gui1)));
% W_mnv = W_guiyi;
% % W_mnv = P*inv_R*A_s*exp(j*pi/9)*1.5;
% W_mnv2 = inv_R*A_s;
% % W_mnv=W_mnv/max(W_mnv);
% amp_w = abs(W_mnv);
% pha_w = angle(W_mnv);

% yy = W_mnv'*X;
% figure();plot(1:ss,yy);
% fft_yy = fft(yy);
% p = unwrap(angle(fft_yy)); 
% figure();plot(1:ss,p);
% % 最大信干噪比MSINR
% Ai = A_i(:,1:interference);
% X_i_n = Ai*S(source+1:source+interference,:);
% R_i_n = X_i_n*X_i_n'/ss;
% W_msinr = inv(R_i_n)*A_s;
% % W_msinr=W_msinr/sqrt(W_msinr'*W_msinr);%归一化

% % 最小均方误差 MMSE
% snr_d = -10;
% % dd = 10.^(snr_d/10)*(randn(1,ss)+j*randn(1,ss));
% % d = awgn(dd,10,'measured');
% d = 10.^(snr_d/10)*exp(-j*w(1)*[0:ss-1]);
% r_xd = X*d'/ss;
% W_mmse = inv_R*r_xd;

%% 阵列方向图  pattern
phi=-89:1:90;
a=exp(-j*pi*(0:N-1)'*sin(phi*pi/180));
F=W_mnv'*a;     % W_mnv  W_msinr  W_mmse
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
% axis equal;
% axis([-90 90 -310 0]);