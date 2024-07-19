clc;
clear;
close all;
format rat;
k=2;
h=1/k;
%Computing  Newton-Cotes rules coefficients for k
syms y r
R=zeros(k,k+1);
coefficients=zeros(k,k+1);
g=@(y,r) y.^r;
I=0:1/k:1;
C=I;
for i=0:k
    R(i+1,:)=g(I,i);
end
for j=1:k
    
    for i=0:k
        C(i+1)=int(g(y,i),y,0,j*1/k);
    end
    coefficients(j,:)=(R\C')';
end
coefficients=coefficientscoefficients*k
