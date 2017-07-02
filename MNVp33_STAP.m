clc;clear all;close all;
 
sou_num=1;  %信源数**
int_num=1;
int_in=1;   %int_in=1 加干扰信号，int_in=0 不加干扰信号
sou_in=0;   %sou_in=1 加信号，sou_in=0 不加信号
noi_in=1;   %noi_in=1 加噪声，noi_in=0 不加噪声

% BD_F=2.49175e9;
BD_F=1.26852e9;
c=3e8;
Fs=62e6;  %采样频率
lamda=c/BD_F;

%BD_lamda=c/BD_F;
%position=[0 0;-82.4456 47.6;-82.4456 47.6;0 -95.2]*0.001;
%position=[0,0;0.5,0;0.25,sqrt(3)/4;-0.25,sqrt(3)/4;-0.5,0;-0.25,-sqrt(3)/4;0.25,-sqrt(3)/4]*lamda;%阵元位置
% position=[1,1;-1,1;-1,-1;1,-1]*lamda/4;%阵元位置   4*2
% position=[1,1;-1,1;-1,-1;1,-1]*lamda/8;%阵元位置

position=[-3;-1;1;3]*lamda/4;
array_num=4;

theta_s=[0 10 20 30]*pi/180;            %信号源到达角
theta_s=theta_s(1:sou_num);
% phi_s=[200 60 180 300]*pi/180;
% phi_s=phi_s(1:sou_num);
% theta_i=[50 20 50 70]*pi/180;         %干扰到达角
theta_i=[50 20 50 70]*pi/180;
theta_i=theta_i(1:int_num);     %1*2
% phi_i=[120 50 210 300]*pi/180;
% phi_i=[120 50 210 300]*pi/180;
% phi_i=phi_i(1:int_num);

Fsi=46.52*1e6;
Fi=[46.52;40;30;50]*1e6;%36.75
Fi=Fi(1:int_num);
%Fi=[36.75 25 20 10]*1e6;   %干扰频率
% lambda_i=[c/Fi(1) c/Fi(2) c/Fi(3) c/Fi(4)];
% lambda_i=lambda_i(1:int_num);
% ux_s=sin(theta_s).*cos(phi_s);
% uy_s=sin(theta_s).*sin(phi_s);
% uz_s=cos(theta_s);%1*1

% ux_i=sin(theta_i).*cos(phi_i);
% uy_i=sin(theta_i).*sin(phi_i);
% uz_i=cos(theta_i);%1*2
%% 阵列流形  STEERING VECTOR

A_i=exp(j*2*pi*position*sin(theta_i)/lamda);%8*4  -j
A_s=exp(j*2*pi*position*sin(theta_s)/lamda);%8*1信号方向40
% As = [A_s A_i];  %(:,1:interference)

% As=exp(j*2*pi*position*[ux_s;uy_s]/lamda);                                 %导向矩阵，远场平面波  4*1
% for i=1:int_num
% Ai(:,i)=exp(j*2*pi*position*[ux_i(i);uy_i(i)]/lamda); %4*2
% end
%Ai=exp(j*2*pi*position*[ux_i;uy_i]/BD_lamda); 

%产生信号及干扰
N=1024; %采样快拍数
ts=1/Fs;
t=(0:N-1)*ts;
% S=100*exp(j*2*pi*Fi*t);
S=exp(j*2*pi*Fsi*t);
% figure();
% plot(1:N,S.');
Si=exp(j*2*pi*Fi*t);%2*N
%产生噪声
  Noise=zeros(array_num,N);
 Noise=sqrt(0.5)*(randn(array_num,N)+j*randn(array_num,N));
 Noise=diag(1./std(Noise,0,2))*Noise;   
 
SNR=[-140 -130 -120]+114;                                                 %信噪比,dB,[]号内为dbw
SNR=SNR(1:sou_num);
INR=[-45 -60 -60 -65]+85;                                                                                                                         
INR=INR(1:int_num);                                                           %干噪比
%阵列接收信号及干扰
SigEnp=zeros(array_num,N);
IntEnp=zeros(array_num,N);
X=zeros(array_num,N);
SigEnp=A_s*diag(10.^(SNR(1:sou_num)/20))*S(1:sou_num,:);
IntEnp=A_i*diag(10.^(INR(1:int_num)/20))*Si(1:int_num,:);% Ai  4*2  Si 2*N      IntEnp 4*N
% SS=As*10.^(SNR(sou_num)/20)*S(1:sou_num,:);
X=SigEnp*sou_in+IntEnp*int_in+Noise*noi_in;
% %--------------------------------------------------------------------------
% 延时为1个时钟
tap_num=3;
NN=array_num*tap_num;
L=N-tap_num+1;
data_st=zeros(NN,L);
% data_st(1:tap_num:NN-tap_num+1,:)=X(:,1:L);
% data_st(2:tap_num:NN,:)=X(:,2:N);
data_st(1:array_num,:)=X(:,1:L);
data_st(array_num+1:2*array_num,:)=X(:,2:L+1);
data_st(2*array_num+1:NN,:)=X(:,3:N);
% %--------------------------------------------------------------------------
% % 延时为2个时钟
% tap_num=3;
% NN=array_num*tap_num;
% L=N-(tap_num-1)*2;
% data_st=zeros(NN,L);
% % data_st(1:tap_num:NN-tap_num+1,:)=X(:,1:L);
% % data_st(2:tap_num:NN,:)=X(:,2:N);
% data_st(1:array_num,:)=X(:,1:L);
% data_st(array_num+1:2*array_num,:)=X(:,3:L+2);
% data_st(2*array_num+1:NN,:)=X(:,5:N);
%--------------------------------------------------------------------------
% % 延时为4个时钟
% tap_num=3;
% NN=array_num*tap_num;
% L=N-(tap_num-1)*4;
% data_st=zeros(NN,L);
% % data_st(1:tap_num:NN-tap_num+1,:)=X(:,1:L);
% % data_st(2:tap_num:NN,:)=X(:,2:N);
% data_st(1:array_num,:)=X(:,1:L);
% data_st(array_num+1:2*array_num,:)=X(:,5:L+4);
% data_st(2*array_num+1:NN,:)=X(:,9:N);
%--------------------------------------------------------------------------
IntF=abs(fft(data_st.'));
fseq=[0:L-1]*Fs/L;
figure();
plot(fseq,20*log10(IntF/max(max(IntF))));
title('原始数据图');

rx=data_st*data_st'/L;
tap_no=[0:1:(tap_num-1)]';
delay_num=1;%延时时钟数
% ast_t=[exp(j*2*pi*Fi*tap_no/Fs)];
ast_t=[exp(j*2*pi*Fsi*tap_no*delay_num/Fs)];
%ast_t=[1;zeros(tap_num-1,1)];
% ast_s=As;
ast_s=[1;1;1;1];
ast=kron(ast_t,ast_s);
%ast=[ones(array_num,1);zeros(NN-array_num,1)];
w=inv(rx)*ast;

y_lcmv=w'*data_st;
IntF_lcmv=abs(fft(y_lcmv.'));
figure();
plot(fseq,20*log10(IntF_lcmv/max(max(IntF_lcmv))));
title('空时算法加权后数据图');

% A_s=exp(j*2*pi*position*sin(theta_s)/lamda);%8*1信号方向40

k=0;
% for phi_se=1:360
%     i=i+1;
%     k=0;
     for the_se=0:1:180
         k=k+1;
         ux_s=sin(the_se*pi/180);uy_s=sin(the_se*pi/180);
         As_se1=exp(j*2*pi*position*sin(theta_s)/lamda); 
         As_se=kron(ast_t,As_se1);
         modu(k)=w'*As_se;
         
     end
% end
modu=20*log10(abs(modu)/max(max(abs(modu))));
figure();
% mesh([0:90],[1:360],modu);
phi=1:1:181;
plot(phi,modu,'linewidth',2);
%  mesh([0:90],[1:360],modu);
% xlabel('俯仰角(\circ)');ylabel('方位角(\circ)');zlabel('阵列增益(dB)');
% title('空时算法抗干扰方向图');
% figure();
% plot(modu(120,:),'g');grid on;zoom on;
% xlabel('俯仰角(\circ)');ylabel('阵列增益(dB)');
% figure();
% plot(modu(:,51),'r');grid on;zoom on;
% xlabel('方位角(\circ)');ylabel('阵列增益(dB)');



