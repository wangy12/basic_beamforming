%% Reference 王永良,丁前军,李荣锋.自适应阵列处理[M]. 北京:清华大学出版社,2009.
%% pp:19-24
%  阵列方向图  
%   1) 阵元数的影响

close all;clear all;clc;
% N=14;                       % 14 array ULA
% phi=[-89:1:90];             % scan the pitch angle
% W=ones(N,1);                % uniformly weighting  均匀加权  [1 0 0 0 0 0 0]';%
% a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));      % steering vector 导向矢量
% F1=W'*a;                    % Patterns of arrays
% % figure();
% % plot(phi,F1);
% G1=abs(F1)/max(abs(F1));    % uniform  归一化
% G_dB1=20*log10(G1);
% figure();
% plot(phi,G_dB1,'linewidth',2);grid on;
% figure();
% polar(phi*pi/180,F1);%xlabel('角度/\circ');ylabel('方向图增益/dB');
% figure();
% plot(phi,G_dB1);%legend('N=8,d=lambda/2');legend('N=15,d=lamda/2');
% axis([-90 90 -50 0]);
% hold on;
% 
% N=7;
% % phi=-89:1:90;
% W=ones(N,1);%均匀加权
% a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));
% F2=W'*a;
% G2=abs(F2)/max(abs(F2));
% G_dB2=20*log10(G2);
% % figure();
% plot(phi,G_dB2,'r--','linewidth',2);
% xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
% % xlabel('角度/\circ');ylabel('方向图增益/dB');    % axis([-90 90 -70 0]);
% % figure();
% % plot(phi,F2,'--');        
% % figure();
% % polar(phi*pi/180,F2,'r--');
% legend('N=14,d=lambda/2','N=7,d=lambda/2');title('influence of array number');
% %  axis([2*pi 15]);
% %  xlabel('角度/\circ');ylabel('方向图增益/dB');
% % grid on;legend('N=7,d=lamda/2');
% 
% figure();
% polar(phi*pi/180,G1);hold on;
% polar(phi*pi/180,G2,'r--');title('influence of array number');% title('阵元个数对方向图的影响');
% 
% %% 2) 空间模糊  孔径相同，阵元间距对方向图的影响
% %  2) the influence of array spacing
% N=8;
% phi=-89:1:90;
% W=ones(N,1);                                %均匀加权
% a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));      % d=lambda/2
% F3=W'*a;
% G3=abs(F3).^2./max(abs(F3).^2);
% G_dB3=10*log10(G3);
% % figure();
% % plot(phi,F3,':');legend('N=8,d=lamda/4');
% figure();
% plot(phi,G_dB3,'linewidth',2);grid on;hold on; axis([-90 90 -50 0]);
% % axis([-90 90 -50 0]);legend('N=8,d=lamda/4');
% 
% N=4;
% phi=-89:1:90;
% W=ones(N,1);                                %均匀加权
% a=exp(j*pi*2*(0:N-1)'*sin(phi*pi/180));    % d=lambda
% F4=W'*a;
% G4=abs(F4).^2./max(abs(F4).^2);
% G_dB4=10*log10(G4);
% % figure();
% % plot(phi,F4,'r--');legend('N=8,d=lamda/4');
% % figure();
% plot(phi,G_dB4,'r--','linewidth',2);legend('N=8,d=lambda/2','N=4,d=lambda');
% xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
% % axis([-90 90 -50 0]);
% title('influence of array spacing with same array apertures');
% 
% figure();
% polar(phi*pi/180,G3);hold on;
% polar(phi*pi/180,G4,'r--');title('孔径相同，阵元间距对方向图的影响');
% 
% %% 波束指向phi=10
% %   3) 低仰角的影响
% %   3) influrnce of pitch angle
% % N=7;
% % phi=-89:1:90;
% % W=ones(N,1);                                            %均匀加权
% % a=exp(j*pi*(0:N-1)'*(sin(phi*pi/180)-sin(pi/18)));     % phi=10
% % F5=W'*a;
% % G5=abs(F5).^2./max(abs(F5).^2);
% % G_dB5=10*log10(G5);
% % figure();
% % plot(phi,G_dB5,'linewidth',2');xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');grid on;
% % axis([-90 90 -50 0]);
% % % figure();
% % % polar(phi*pi/180,F5);
% % hold on;
% % % 波束指向  phi=80
% % N=7;
% % phi=-89:1:90;
% % W=ones(N,1);%均匀加权
% % a=exp(j*pi*(0:N-1)'*(sin(phi*pi/180)-sin(pi*6/18)));   % phi=60
% % F6=W'*a;
% % G6=abs(F6).^2./max(abs(F6).^2);
% % G_dB6=10*log10(G6);
% % % figure();
% % plot(phi,G_dB6,'r--','linewidth',2);legend('N=7,d=lambda/2');grid on;
% % % axis([-90 90 -50 0]);
% % title('influrnce of pitch angle')
% % 
% % figure();
% % polar(phi*pi/180,G5);hold on;
% % polar(phi*pi/180,G6,'r--');title('influrnce of pitch angle')

N=8;
phi=-179:1:180;
W=exp(j*pi*(0:N-1)'*(sin(10*pi/180))); % CBF % phi=5
a=exp(j*pi*(0:N-1)'*(sin(phi*pi/180)));     
F5=W'*a;
G5=abs(F5).^2./max(abs(F5).^2);
G_dB5=10*log10(G5);
figure();
plot(phi,G_dB5,'linewidth',2');xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');grid on;
axis([-180 180 -50 0]);
% figure();
% polar(phi*pi/180,F5);
hold on;
% 波束指向  phi=80
N=7;
phi=-89:1:90;
W=exp(j*pi*(0:N-1)'*(sin(50*pi/180)));     % CBF % phi=50
a=exp(j*pi*(0:N-1)'*(sin(phi*pi/180)));   
F6=W'*a;
G6=abs(F6).^2./max(abs(F6).^2);
G_dB6=10*log10(G6);
% figure();
plot(phi,G_dB6,'r--','linewidth',2);legend('N=7,d=lambda/2');grid on;
% axis([-90 90 -50 0]);
title('influrnce of pitch angle')

figure();
polar(phi*pi/180,G5);hold on;
% polar(phi*pi/180,G6,'r--');title('influrnce of pitch angle')
%% 
%  4) 加窗影响 低副瓣锥化
%  4) window
N=21;
w_cheb = chebwin(N,30);     % tapering vector -- -30dB chebyshev win
W=ones(N,1).*w_cheb;%均匀加权
a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));
F7=W'*a;
G7=abs(F7)/max(abs(F7));
G_dB7=20*log10(G7);

N=21;
% phi=-89:1:90;
W=ones(N,1);%均匀加权
a=exp(j*pi*(0:N-1)'*sin(phi*pi/180));
F2=W'*a;
G2=abs(F2)/max(abs(F2));
G_dB2=20*log10(G2);
figure();
plot(phi,G_dB2,'r--','linewidth',2);hold on
plot(phi,G_dB7,'linewidth',2);xlabel('Picth Angle (\circ)');ylabel('Magnitude (dB)');
grid on;axis([-90 90 -60 0]);
title('windowing and tapering');
% figure();
% plot(phi,F2,'--');        
% figure();
% polar(phi*pi/180,F2,'r--');
% legend('N=15,d=lambda/2','N=7,d=lambda/2');title('阵元个数对方向图的影响');
%  axis([2*pi 15]);
%  xlabel('角度/\circ');ylabel('方向图增益/dB');
% legend('N=7,d=lamda/2');

% figure();
% polar(phi*pi/180,G7);   % hold on;
% polar(phi*pi/180,G2,'r--');title('阵元个数对方向图的影响');