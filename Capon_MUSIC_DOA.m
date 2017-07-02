%% 七阵元Capon MUSIC 波束形成
clc;clear all;close all;
%% 参数设定
% % sensor=11;          %传感器11阵元  阵元数
sensor=25;          %传感器7阵元  阵元数
% sensor=6;          %传感器19阵元  阵元数
source=5;           %信源  卫星数
%source=8;           %信源  卫星数
j=sqrt(-1);
Freq=46.5e6;
c=3e8;
lambda=c/Freq;
snapshot=2048;      %快拍
snr=[30 30 30 30 40];  %信噪比不同，波束形成上的影响
theta=[17 30 50 70 50]*pi/180;   %信号源到达角 俯仰角
phi=[243 100 150 200 250]*pi/180;    %水平角
%% 阵列排布
% 线阵
% 7ULA
% position=[-1.5,0;-1,0;-0.5,0;0,0;0.5,0;1,0;1.5,0]*lambda;% ;2,0   ;2,0;0,0.5;0,1
% 2 ULA
% position=[0,0;0.5,0]*lambda;
% 4 ULA
% position=[0,0;0.5,0;1,0;1.5,0]*lambda;
% 7 HEX ARRAY
% position=[0,0;0.5,0;0.25,sqrt(3)/4;-0.25,sqrt(3)/4;-0.5,0;-0.25,-sqrt(3)/4;0.25,-sqrt(3)/4]*lambda;%阵元位置 七

% 5+4 L ARRAY
% position=[0,0;0.5,0;1,0;1.5,0;2,0;0,0.5;0,1;0,1.5;0,2]*lambda;% ;2,0   ;2,0;0,0.5;0,1

% 5+4 十 ARRAY
% position=[-1,0;-0.5,0;0,0;0.5,0;1,0;0,0.5;0,1;0,-0.5;0,-1]*lambda;% ;2,0
% ;2,0;0,0.5;0,1

% 5X5方针
position=[-1,0;-0.5,0;0,0;0.5,0;1,0;....
          -1,0.5;-0.5,0.5;0,0.5;0.5,0.5;1,0.5;....
          -1,1;-0.5,1;0,1;0.5,1;1,1;....
          -1,-0.5;-0.5,-0.5;0,-0.5;0.5,-0.5;1,-0.5;....
          -1,-1;-0.5,-1;0,-1;0.5,-1;1,-1]*lambda;
      
% 19阵元  圆阵两圈结合正六边形
% position=[0,0;1/2,0;1,0;0.25,sqrt(3)/4;3/4,sqrt(3)/4;0.5,sqrt(3)/2;....
%     0,1;-1/2,0;-1,0;-0.25,sqrt(3)/4;-3/4,sqrt(3)/4;-0.5,sqrt(3)/2;....
%     -0.25,-sqrt(3)/4;-3/4,-sqrt(3)/4;-0.5,-sqrt(3)/2;....
%     0,-1;0.25,-sqrt(3)/4;3/4,-sqrt(3)/4;0.5,-sqrt(3)/2]*lambda;

% position=[0,0;0.5,0;0.25,1/(sqrt(3)*4);0.25,3/(sqrt(3)*4);-0.5,0;-0.25,1/(sqrt(3)*4);-0.25,3/(sqrt(3)*4);....
% 0.25,-1/(sqrt(3)*4);0.25,-3/(sqrt(3)*4);-0.25,-1/(sqrt(3)*4);-0.25,-3/(sqrt(3)*4)]*lambda;
                     % 阵元位置 半径为1/2波长的七阵子圆阵加正三角形中心位置
% position=[0,0;sqrt(3)/2,0;sqrt(3)/4,0.25;sqrt(3)/4,3/4;-sqrt(3)/2,0;-sqrt(3)/4,0.25;-sqrt(3)/4,3/4;sqrt(3)/4,-0.25;....
%     sqrt(3)/4,-3/4;-sqrt(3)/4,-0.25;-sqrt(3)/4,-3/4]*lambda;
                    % 阵元位置 半径为1/2波长的五阵子圆阵加正三角形定点位置
% 19阵元  六边形阵  
% position = [0,0; 1/2,0; 1,0; 0.25,sqrt(3)/4; 3/4,sqrt(3)/4; 0.5,sqrt(3)/2;...
%     0,sqrt(3)/2; -1/2,0; -1,0; -0.25,sqrt(3)/4; -3/4,sqrt(3)/4; -0.5,sqrt(3)/2;...
%     -0.25,-sqrt(3)/4; -3/4,-sqrt(3)/4; -0.5,-sqrt(3)/2;...
%     0,-sqrt(3)/2; 0.25,-sqrt(3)/4; 3/4,-sqrt(3)/4;
%     0.5,-sqrt(3)/2]*lambda;

figure();plot(position(:,1),position(:,2),'s');
xlabel('X');ylabel('Y');title('阵列排布方式');
% axis tight;
grid on;
%% 坐标换算
ux=cos(theta).*sin(phi);
uy=cos(theta).*cos(phi);
%uz=cos(theta);
%% 导向矩阵/方向向量/响应向量
%A_sig=exp(-j*2*pi*position*[ux_s;uy_s]/lamda);  %导向矩阵，远场平面波7*1
for i=1:source
    A(:,i)=exp(-j*2*pi*position*[ux(i);uy(i)]/lambda);  %sensor*source
end
%% 噪声
%Noise=zeros(array_num,N);
Noise=(randn(sensor,snapshot)+j*randn(sensor,snapshot));
%% 信号复包络
% S=2*randn(source,snapshot);
f = [200 100 500 300 465];
w=2*pi./f;
% w=[pi/6 pi/5 pi/4 0 pi/2]';
for i=1:1:source-1
    S(i,:)=10.^(snr(i)/10)*exp(-j*w(i)*(0:snapshot-1));  %source*snapshot
end
S(i+1,:)=10.^(snr(i+1)/10)*(randn(1,snapshot)+j*randn(1,snapshot));
%d=10.^(snr/10)*exp(-j*2*pi*GPS_F*t);%1*512
%% 协方差矩阵
X=A*S;%+Noise;%7*512
Rx=X*X'/snapshot;%7*7
%rxd=X*d'/snapshot;%7*1
inv_Rx = inv(Rx);
%% CAPON 波束形成
angle_step=1;
% theta_pos=0;
phi_pos=0;
for phi_s=0:angle_step:360
    phi_pos=phi_pos+1;
    theta_pos=0;
    for theta_s=0:angle_step:90
      theta_pos=theta_pos+1;  
%       ux_s=sin(theta_s*pi/180).*cos(phi_s*pi/180);
%       uy_s=sin(theta_s*pi/180).*sin(phi_s*pi/180);
        ux_s=cos(theta_s*pi/180).*sin(phi_s*pi/180);
      uy_s=cos(theta_s*pi/180).*cos(phi_s*pi/180);
      A_temp=exp(-j*2*pi*position*[ux_s;uy_s]/lambda);
      Pcapon(phi_pos,theta_pos)=1/(A_temp'*inv_Rx*A_temp);
%       P(phi_pos,theta_pos)=abs(A_temp'*inv_Rx*A_temp);
    end
end

figure();
mesh([0:angle_step:90],[0:angle_step:360],abs(Pcapon));
xlabel('pitch/ \circ');ylabel('azimuth/ \circ');zlabel('magnitude/dB');title('Capon method DOA estimation');

Max1=max(max(abs(Pcapon)));
G1=abs(Pcapon/Max1);%归一化
G1_dB=10*log10(G1);
figure();
mesh([0:angle_step:90],[0:angle_step:360],G1_dB);% colorbar;
xlabel('pitch/ \circ');ylabel('azimuth/ \circ');zlabel('magnitude/dB');title('Capon method DOA estimation');
% PP=G1_dB(:,20);
% figure();plot([0:angle_step:360],G1_dB(:,60),'r');grid on;
% figure();plot([0:angle_step:90],G1_dB(200,:),'k');grid on;
% figure();plot([0:angle_step:360],G1_dB(200,:),'r');% grid on;
%% MUSIC
[Vec,Val]=eig(Rx);                              %特征向量、特征值
[Val Seq]=sort(max(Val));                       % 从小到大排序
Vec_s=Vec(:,Seq(sensor-source+1:sensor));      %信源个数估计sensor-GDE_num+1:sensor 信号特征向量，明显大;取5到7列
Vec_n=Vec(:,Seq(1:sensor-source));             %1:sensor-GDE_num 噪声特征向量  7*4；取1到4列

%四阵元
%Vec_s=Vec(:,Seq(2:4));%信源个数估计sensor-GDE_num+1:sensor信号特征向量，明显大;取5到7列
%Vec_n=Vec(:,Seq(1:1));%1:sensor-GDE_num噪声特征向量  7*4；取1到4列

theta_pos=0;
phi_pos=0;
for theta_s=0:angle_step:360
    theta_pos=theta_pos+1;
    phi_pos=0;
    for phi_s=0:angle_step:90
      phi_pos=phi_pos+1;  
      ux_s=cos(theta_s*pi/180).*sin(phi_s*pi/180);
      uy_s=sin(theta_s*pi/180).*sin(phi_s*pi/180);
      As_temp=exp(-j*2*pi*position*[ux_s;uy_s]/lambda);
      B=As_temp'*Vec_n*Vec_n'*As_temp;
      Pmusic(theta_pos,phi_pos)=1/B;         %real(1/B)
    end 
end
figure();
mesh([0:angle_step:90],[0:angle_step:360],abs(Pmusic));
xlabel('pitch/ \circ');ylabel('azimuth/ \circ');zlabel('magnitude/dB');title('MUSIC method DOA estimation');
P_guiyi=10*log10(abs(Pmusic)/max(max(abs(Pmusic))));
figure();
mesh([0:angle_step:90],[0:angle_step:360],P_guiyi);% colorbar;
xlabel('pitch/ \circ');ylabel('azimuth/ \circ');zlabel('magnitude/dB');title('MUSIC method DOA estimation');

figure();plot([0:angle_step:360],P_guiyi(:,60),'r');grid on;
figure();plot([0:angle_step:90],P_guiyi(200,:),'k');grid on;