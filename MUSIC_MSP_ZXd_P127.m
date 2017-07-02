%% 七阵元Capon MUSIC 波束形成
clc;clear all;close all;
%% 参数初始化
sensor=7;       %传感器七阵元
%sensor=4;      %四阵元
source=2;       %信源
j=sqrt(-1);
GPS_F=1.57542e9;
c=3e8;
Fs=62e6;  %采样频率
lamda=c/Fs;  %采样的波长？？lambda
ss=1024;%快拍
ts=1/Fs;
t=(0:ss-1)*ts;
%snr=10;
snr=[10 10 10];%信噪比不同，波束形成上的影响？
position=[0,0;0.5,0;0.25,sqrt(3)/4;-0.25,sqrt(3)/4;-0.5,0;-0.25,-sqrt(3)/4;0.25,-sqrt(3)/4]*lamda;%阵元位置 七

%position=[0,0;-sqrt(3)/2,1/2;sqrt(3)/2,1/2;0,-1]*lamda/2;%四阵元

theta=[60 40 80]*pi/180;   %信号源到达角 俯仰角
phi=[160 250 160]*pi/180;    %水平角

Fi=[1.575 1 2.5]*1e9;   %什么频率?频率有何影响？
%lamda_i=[c/Fi(1) c/Fi(2) c/Fi(3)];%干扰波长
%% 坐标换算
ux=sin(theta).*cos(phi);
uy=sin(theta).*sin(phi);
uz=cos(theta);
%% 导向矩阵/方向向量/响应向量
%A_sig=exp(-j*2*pi*position*[ux_s;uy_s]/lamda);  %导向矩阵，远场平面波7*1
for i=1:source;
    A(:,i)=exp(-j*2*pi*position*[ux(i);uy(i)]/lamda);%7*3
end
% 期望方向
% A_desired=[1;0;0;0;0;0;0];%7*3
theta_desired=[20]*pi/180;   %信号源到达角 俯仰角
phi_desired=[300]*pi/180;    %水平角
ux_desired=sin(theta_desired).*cos(phi_desired);
uy_desired=sin(theta_desired).*sin(phi_desired);
A_desired = exp(-j*2*pi*position*[ux_desired;uy_desired]/lamda);
%% 噪声
%Noise=zeros(array_num,N);
Noise=randn(sensor,ss)+j*randn(sensor,ss);%7*512
%% 信号复包络
for i=1:1:source
    S(i,:)=10.^(snr(i)/10)*exp(-j*2*pi*Fi(i)'*t);%3*512
end
%d=10.^(snr/10)*exp(-j*2*pi*GPS_F*t);%1*512
%% 协方差矩阵
X=A*S+Noise;%7*512
Rx=X*X'/ss;%7*7
%rxd=X*d'/ss;%7*1
Rx_inv = inv(Rx);

W_opt=Rx_inv*A_desired/(A_desired'*Rx_inv*A_desired);%计算权值
W=W_opt/sqrt(W_opt'*W_opt);%归一化
anti_jam_X = [W(1,1)*X(1,:);W(2,1)*X(2,:);W(3,1)*X(3,:);W(4,1)*X(4,:);W(5,1)*X(5,:);W(6,1)*X(6,:);W(7,1)*X(7,:)];
Rx_anti=anti_jam_X*anti_jam_X'/ss;%7*7
Rx_inv_anti = inv(Rx_anti);

%% CAPON 波束形成
angle_step=1;
theta_pos=0;
phi_pos=0;
for phi_s=1:angle_step:360
    phi_pos=phi_pos+1;
    theta_pos=0;
    for theta_s=1:angle_step:90
      theta_pos=theta_pos+1;  
      ux_s=sin(theta_s*pi/180).*cos(phi_s*pi/180);
      uy_s=sin(theta_s*pi/180).*sin(phi_s*pi/180);
      A_temp=exp(-j*2*pi*position*[ux_s;uy_s]/lamda);
      Pcapon(phi_pos,theta_pos,:)=1/(A_temp'*Rx_inv*A_temp);
%       P(phi_pos,theta_pos,:)=W'*A_temp;
        P(phi_pos,theta_pos)=A_temp'*Rx_inv*A_temp;
      end
end

Max2=max(max(abs(Pcapon)));
G2=abs(Pcapon)/Max2;%归一化
% G1_dB=20*log10(G1);
figure();
mesh([1:angle_step:90],[1:angle_step:360],G2);
colorbar;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Pcapon/dB');title('Capon方法测向谱');
figure();
contour([1:angle_step:90],[1:angle_step:360],G2);grid on;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('Capon方法测向谱俯视图');

% figure();
% meshc([1:angle_step:90],[1:angle_step:360],20*log10(abs(P)));
% 
Max1=max(max(abs(Pcapon)));
G1=abs(Pcapon)/Max1;%归一化
G1_dB=20*log10(G1);
figure();
mesh([1:angle_step:90],[1:angle_step:360],G1_dB);
colorbar;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Pcapon/dB');title('Capon方法测向谱');
figure();
contour([1:angle_step:90],[1:angle_step:360],G1_dB);grid on;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('Capon方法测向谱俯视图');

figure();
mesh([1:angle_step:90],[1:angle_step:360],abs(P));
colorbar;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Pcapon/dB');title('Capon方法测向谱');
figure();
contour([1:angle_step:90],[1:angle_step:360],P);grid on;
xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('Capon方法测向谱俯视图');
% 
% figure();meshc([1:angle_step:90],[1:angle_step:360],20*log10(abs(P)));%colorbar;%rmal
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('方向图增益/dB');title('无约束时抗干扰方向图');
% figure();contour([1:angle_step:90],[1:angle_step:360],20*log10(abs(P)));grid on;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('无约束时抗干扰方向图俯视图');
% 
% 
% angle_step=1;
% theta_pos=0;
% phi_pos=0;
% for phi_s=1:angle_step:360
%     phi_pos=phi_pos+1;
%     theta_pos=0;
%     for theta_s=1:angle_step:90
%       theta_pos=theta_pos+1;  
%       ux_s=sin(theta_s*pi/180).*cos(phi_s*pi/180);
%       uy_s=sin(theta_s*pi/180).*sin(phi_s*pi/180);
%       A_temp=exp(-j*2*pi*position*[ux_s;uy_s]/lamda);
%       Pcapon_anti(phi_pos,theta_pos,:)=1/(A_temp'*Rx_inv_anti*A_temp);
%       end
% end
% 
% Max3=max(max(abs(Pcapon_anti)));
% G3=abs(Pcapon_anti)/Max3;%归一化
% figure();
% meshc([1:angle_step:90],[1:angle_step:360],G3);
% colorbar;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Pcapon/dB');title('Capon方法后测向谱');
% figure();
% contour([1:angle_step:90],[1:angle_step:360],G3);grid on;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('Capon方法后测向谱俯视图');
% 
% Max4=max(max(abs(Pcapon_anti)));
% G4=abs(Pcapon_anti)/Max4;%归一化
% G4_dB=20*log10(G4);
% figure();
% mesh([1:angle_step:90],[1:angle_step:360],G4_dB);
% colorbar;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Pcapon/dB');title('Capon方法后测向谱');
% figure();
% contour([1:angle_step:90],[1:angle_step:360],G4_dB);grid on;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('Capon方法后测向谱俯视图');



%% MUSIC
%参考张星doa_music
% 
% [Vec,Val]=eig(Rx);%特征向量、特征值
% [Val Seq]=sort(max(Val));
% %七阵元
% Vec_s=Vec(:,Seq(sensor-source:7));%信源个数估计sensor-GDE_num+1:sensor信号特征向量，明显大;取5到7列
% Vec_n=Vec(:,Seq(1:sensor-source));%1:sensor-GDE_num噪声特征向量  7*4；取1到4列  sensor-source
% %四阵元
% %Vec_s=Vec(:,Seq(2:4));%信源个数估计sensor-GDE_num+1:sensor信号特征向量，明显大;取5到7列
% %Vec_n=Vec(:,Seq(1:1));%1:sensor-GDE_num噪声特征向量  7*4；取1到4列
% 
% theta_pos=0;
% phi_pos=0;
% for theta_s=0:angle_step:360
%     theta_pos=theta_pos+1;
%     phi_pos=0;
%     for phi_s=0:angle_step:90
%       phi_pos=phi_pos+1;  
%       ux_s=cos(theta_s*pi/180).*sin(phi_s*pi/180);
%       uy_s=sin(theta_s*pi/180).*sin(phi_s*pi/180);
%       As_temp=exp(-j*2*pi*position*[ux_s;uy_s]/lamda);%7*1
%       B=(As_temp'*Vec_n*Vec_n'*As_temp);
%       P(theta_pos,phi_pos)=1/B;%real(1/B)
%     end 
% end
% P_guiyi=20*log10(abs(P)/max(max(abs(P))));
% figure();
% mesh([0:angle_step:90],[0:angle_step:360],P_guiyi);colorbar;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');zlabel('Parten/dB');title('MUSIC方法测向谱');
% figure();
% contour([0:angle_step:90],[0:angle_step:360],P_guiyi);grid on;
% xlabel('俯仰角/ \circ');ylabel('方位角/ \circ');title('MUSIC方法测向谱俯视图');

%% 如何确定信源个数？盖尔圆法估计信号源个数?


