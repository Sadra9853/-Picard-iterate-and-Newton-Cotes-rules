clc; clear; close all;
format short; format compact;

%% =====================  Global Variables  ============================
global I K n k;

syms t s u

%% =====================  Example 6 Definition  ========================
exact = @(t) 5*t;
f     = @(t) log(26/25)/10 - log(t.^2 + 1/25)/10 + 5*t;
K     = @(t,s,u) 1./(1 + 25*s.^2) * u;

n = 100;
k = 10;


%% ===================== Interval and Step Size ========================
a = -1;
b = 1;
h = (b - a) / n;
%% ===================== Computing Newton–Cotes Coefficients ===========
syms y
R  = zeros(k, k+1);
ww = zeros(k, k+1);
g = @(y,r) y.^r;
Igrid = 0 : 1/k : 1;
C = Igrid;

for i = 0:k
    R(i+1, :) = g(Igrid, i);
end

for j = 1:k
    for i = 0:k
        C(i+1) = int(g(y, i), y, 0, j/k);
    end
    ww(j,:) = (R \ C')';
end

coefficients = k * h * ww;

%% ===================== Numerical Solution Initialization =============

I  = a:h:b;      % global I used in NC(), NCEND()
xs = f(I);      %intial numerical solution
ys = xs .* k;
xs(1) = f(I(1));
L = I;

%% ===================== Main Iterative Loop ============================
tic
it1 = 0;
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

%% ===================== Compute Exact and Error ========================
L = exact(I);
e = abs(L - xs);
M = max(e)
figure;
plot(I, e, 'r*');
title('Error Plot');
xlabel('t');
ylabel('Error');
%% ===================== Save ========================
% save(['D:\09127944942\amanuscript\Volltera  one\Matlab\Example6\k=100\n=10\Voltteraexample.mat'], ...
    % 'e','M','it1','xs','L','n','k');

%% ===================== Inside Newton–Cotes Function ==================
function s3 = NC(Array, location_s, n, k, w)
global I K
SS = 0;
for j1 = 1:k:n-k
    for i1 = 0:k
        SS = SS + w(i1+1) * K( I(location_s), I(j1+i1), Array(j1+i1) );
    end
end
s3 = SS;
end
%% ===================== Outside Newton–Cotes Function =================
function s3 = NCEND(Array, location_s, Inew, w2)
global I K k
SS = 0;

for i1 = 0:k
    SS = SS + w2(i1+1) * K( I(location_s), Inew(i1+1), Array(i1+1) );
end

s3 = SS;
end
