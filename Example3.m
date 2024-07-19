clc;
clear;
close all;
format short;
format compact;
global I K n  k;
syms t s u
%Example3
exact=@(t) -(6*t + 1)/(3*t + 2);
f=@(t) (exp(-4*t)/2-1) ;
K=@(t,s,u) -exp(-4*(t-s))*u^2;
%Example3
n =20; %n must be a multip multiple of k
k=2;
%Computing  Newton-Cotes rules coefficients for k
syms y r
coefficients=zeros(k,k+1);
w=zeros(k,k+1);
g=@(y,r) y.^r;
I=0:1/k:1;
C=I;
a=0;
b=1;
for i=0:k
    coefficients(i+1,:)=g(I,i);
end
for j=1:k
    
    for i=0:k
        C(i+1)=int(g(y,i),y,0,j*1/k);
    end
    w(j,:)=(coefficients\C')';
end
h=(b - a)/n;
coefficients=k*h*w;
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
% save('C:\09127944942\amanuscript\Fractional\Volltera  one\Matlab\Example3\n=20\k=2\Voltteraexample.mat','e','M','it1','xs','L','n','k');
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
%Computing  Re-NCR n,r,q (K(t r ,·,x(·)),[t 0 = a,t r = a + rh])Iterative method

