clc;
clear;
close all;
format rat;
k=10;
h=1/k;

syms t f  factor(t)
ww=zeros(k,k+1);
factor(t)=1;
   
for i=0:k
    factor(t)=factor(t)*(t-i);
end

for q=1:k
   
for i=0:k
    ww(q,i+1)=(-1)^(k-i)/(factorial(i)*factorial(k-i))*int(factor(t)/(t-i),0,q);
end
end
save('k10','ww');
