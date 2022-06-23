clear;
clc;
format shortE %important announcement, because the default accuracy in MATLAB is not enough

d = [0.2388,0.2319,0.2252,0.2188,0.2126,0.2066,0.2008,0.1952,0.1898,0.1846,0.1795,0.1746,0.1699,0.1654,0.1610,0.1567,0.1526,0.1486,0.1447,0.1410]';
y = (0.025:0.05:0.9750)';
% discretize the x
x = (0.025:0.05:0.9750)';
dx = 0.05;

% caculate G
G = zeros(20,length(x));
for a = 1:length(y)
    for j = 1:length(x) %he column of G matrix is related to the length of x
        G(a,j) = dx*x(j)*exp(-x(j)*y(a));
    end
end

%calculate the unknown m(x)
m = G\d;
m1 = pinv(G)*d;
d1 = G*m;
d2 = G*m1;
z = (1:20);%number of node
figure(1)
plot(z,d,'g+',z,d1,'ro',z,d2,'b.'),legend('d(original)','d_1(m=G^-1d)','d_2(m=G_g^-1d)')
xlabel('i')
ylabel('d')
figure(2)
plot(z,m,'o')
xlabel('i')
ylabel('m(i)')
figure(3)
plot(z,m1,'o')
xlabel('i')
ylabel('m(i)')

%to calculate the singular value and the condition number of G
[U,S,V] = svd(G);
Singular_matrix = S;
condition_number = cond(G,2);

%make the plot of picard condition 
for a=1:20
    Ud(a) = abs(U(:,a)'*d); 
end
Si = diag(S); 
figure(8)
plot(z,Si,'o'),xlabel('i'),ylabel('S_i'),set(gca,'Yscale','log')
for c = 1:20
   Ud1(c) = Ud(c)/Si(c); 
end

figure(4)
plot(z,Si,'r+',z,Ud,'b.',z,Ud1,'g*'),legend('S_i','|U^T_.,i*d|','|U^T_.,i*d|/S_i')
set(gca,'Yscale','log')
title('Discrete Picard Condition')
xlabel('i')
set(gca,'yminortick','on') 

% using TSVD
for ii = 1:4
    S1(ii,ii) = S(ii,ii);
end

for jj = 1:4
    V1(:,jj) = V(:,jj);
end

for kk = 1:4
    U1(:,kk) = U(:,kk);
end

mest = V1*inv(S1)*U1'*d;
figure(5)
subplot(211)
plot(z,mest,'o'),legend('m(i)')
xlabel('i'),ylabel('m(i)')
d3 = G*mest;
subplot(212)
plot(z,d,'.',z,d3,'o'),legend('original d','d by TSVD')
xlabel('i'),ylabel('d')

% choose an arbitrary e to calculate 
e = 10^(-2);
 mest2 = inv(G'*G+e^2*I)*G'*d;
 d4 = G*mest2;
 figure(9)
 subplot(211)
 plot(z,d,'.',z,d3,'o'),legend('d(original data)','d(TSVD-form)')
 xlabel('i'),ylabel('d')
 subplot(212)
plot(z,d,'.',z,d4,'o'),legend('d(original data)','d(SVD-form e = 10^-2)')
xlabel('i'),ylabel('d')

