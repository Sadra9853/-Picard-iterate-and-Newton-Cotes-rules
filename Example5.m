clc;
clear;
close all;
format short;
format compact;
global I K n  k;
syms t s u
%Example5
syms t s u
exact=@(t) (4-t^2)^2/16;
f=@(t) 1+t^2;
K=@(t,s,u) -s*(sqrt(u)+2);
%Example5
n =40;
k=10;
%Computing  Newton-Cotes rules coefficients for k
syms y r
R=zeros(k,k+1);
ww=zeros(k,k+1);
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
    ww(j,:)=(R\C')';
end
a=0;
b=1;
h=(b - a)/n;
coefficients=k*h*ww;
%Computing  Newton-Cotes rules coefficients for k
SS=0;
xs=a:h:b;%Numerical slotion
ys=xs.*k;
I=a:h:b;
it1=0;
L=a:h:b;
xs(1)=f(I(1));
%Main loop
tic
while max(abs(ys-xs))>10^(-15)
    ys=xs;
    it1=it1+1
    for s=2:n+1
        a=mod(s,k);
        if a==0 || a==1
            a=a+k;
        end
        xs(s)=f(I(s))+NC(xs(1:s-a+1),s,s-a+1,k,coefficients(k,:))+NCEND(xs(s-a+1:s-a+k+1),s,I(s-a+1:s-a+k+1),coefficients(a-1,:));
    end
end
toc
%Main loop
for s=1:n+1
    L(s)=exact(I(s));
end
e=abs(L-xs);
M=max(e)
%save(['Voltteraexample.mat','e','M','it1','xs','L','n','k');
 plot(I,e,'r*');
%Computing  Re-NCR n,r,q (K(t r ,·,x(·)),[t 0 = a,t r = a + rh])Iterative method based on inside Newton-Cotes rules (equation 16) for q=0
function s3 = NC(Array,location_s,n,k,w)
global  I K
SS=0;
for j1=1:k:n-k
    for i1=0:k
        
        SS=SS+w(i1+1)*K(I(location_s),I(j1+i1),Array(j1+i1));
        
    end
end
s3= SS;
end
%Computing  Re-NCR n,r,q (K(t r ,·,x(·)),[t 0 = a,t r = a + rh])Iterative method based on inside Newton-Cotes rules (equation 16) for q=0

%Computing  Re-NCR n,r,q (K(t r ,·,x(·)),[t 0 = a,t r = a + rh])Iterative method based on outside Newton-Cotes rules (equation 16) for q=1,...,k
function s3 = NCEND(Array,location_s,Inew,w2)
global  I K k
SS=0;

for i1=0:k
    
    SS=SS+w2(i1+1)*K(I(location_s),Inew(i1+1),Array(i1+1));
    
end

s3= SS;
end
%Computing  Re-NCR n,r,q (K(t r ,·,x(·)),[t 0 = a,t r = a + rh])Iterative method based on outside Newton-Cotes rules (equation 16) for q=1,...,k
