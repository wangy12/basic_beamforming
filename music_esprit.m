%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%design by Hangzhou Dianzi University
%hxj
%AR model parameter estimate 
clear;clc;close all;
m=sqrt(-1);
delta=0.101043;
a1=-0.850848;
sample=32;      %number of sample spot
p=10;           %number of sample spot in coef method;
f1=0.05; f2=0.40; f3=0.42;
fstep=0.01;
fstart=-0.5;
fend=0.5;  
f=fstart:fstep:fend;
nfft=(fend-fstart)/fstep+1;

%un=urn+juin
urn= normrnd(0,delta/2,1,sample);
uin= normrnd(0,delta/2,1,sample);
un=urn+m*uin;

%计算 zn
for n=1:sample-1
    zn(1)=un(1);
    zn(n+1)=-a1*zn(n)+un(n+1);
end
    
%计算 xn
    for n=1:sample
        xn(n)=2*cos(2*pi*f1*(n-1))+2*cos(2*pi*f2*(n-1))+2*cos(2*pi*f3*(n-1))+sqrt(2)*real(zn(n));
    end
    
%计算 rxx
for k=0:p
    s=0;
    for n=1:sample-k
        s=s+1/sample*(conj(xn(n))*xn(n+k));
    end
    rxx(k+1)=s;
end
rx=toeplitz(rxx(1:p));
b=reshape(rxx,p+1,1);
a=-inv(rx)*b(2:p+1);

%计算 psd of xn
for i=1:length(f)
   sum=0;
   for k=1:p
       sum=sum+a(k)*exp(-m*2*pi*f(i)*k);
   end
  Pacoef(i)=delta/(abs(1+sum))^2;
end
figure
semilogy(f,Pacoef);
title('人工计算自相关');

