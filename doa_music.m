clear all;close all;clc;

source_num=1;
sensor_num=8;
B=4;%what is B?
N=1024;
snapshot_num=N;
w=[pi/6 pi/3 pi/2]';
l=((2*pi*3e8)/w(1)+(2*pi*3e8)/w(2)+(2*pi*3e8)/w(3))/3;
d=l/2;
snr=30;
source_doa=[pi/3 pi/4 pi/6];
j=sqrt(-1);

A=[exp(-j*(0:sensor_num-1)*d*2*pi*cos(source_doa(1))/l);exp(-j*(0:sensor_num-1)*d*2*pi*cos(source_doa(2))/l);exp(-j*(0:sensor_num-1)*d*2*pi*cos(source_doa(3))/l)].';%8*3

s=sqrt(10.^(snr/10))*exp(j*w*[0:N-1]);%  signal   3*1024

x=A*s+(1/sqrt(2))*(randn(sensor_num,N)+j*randn(sensor_num,N));%signal+noise   8*1024
%the majority problem lies here   T why not x?
T=1/sqrt(sensor_num)*[exp(-j*(0:sensor_num-1)'*pi*sin(1/sensor_num)) exp(-j*(0:sensor_num-1)'*pi*sin(2/sensor_num)) exp(-j*(0:sensor_num-1)'*pi*sin(4/sensor_num)) exp(-j*(0:sensor_num-1)'*pi*sin(6/sensor_num))];%8*4              
%T=1/sqrt(sensor_num)*x;%8*1024
T1=T'*T;%4*4 or 6*6
disp(T1);
y=T'*x;%4*8-8*1024  6*8-8*1024
R=y*y'/N;%4*4 6*6
[V,D]=eig(R);
D=diag(D);
disp(D);
Un=V(:,1:B-source_num);
disp(Un);
Gn=Un*Un';
disp(Gn);

searching_doa=-89:0.1:90;
 for i=1:length(searching_doa)
   a_theta=exp(-j*(0:sensor_num-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/l);
   Pmusic(i)=1./abs((a_theta)'*T*Gn*T'*a_theta);
 end
plot(searching_doa,10*log(Pmusic));
grid;

%t=0:0.1:8;
%theta(n)=[30 45 60 90]
%S(n,t)=[2*sin(pi*t) 2*sin(pi*t) 2*sin(pi*t) 2*sin(pi*t)];%related with theta?
%for i=1:16
%   N(i,t)=randn(size(t));
%end
%j=sqrt(-1);
%X(i,t)=zeros(i,t);
%for i=1:16        %l=lamda/2  array unit---i
%   for n=1:4     %signal unit --- n
%       X(i,t)=X(i,t)+S(n,t)*exp(-j*pi*(i-1)*cos(theta(n)))
%   end
%    X(i,t)=X(i,t)+N(i,t)
%end
