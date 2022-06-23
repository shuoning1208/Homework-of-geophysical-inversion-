clear;
clc;
%initial parameter
dx=0.05;
p=1:20;
N=20;

x=(0.025:0.05:0.975);
y=(0.025:0.05:0.975);
d = [0.2388 0.2319 0.2252 0.2188 0.2126 0.2066 0.2008 0.1952 0.1898 0.1846 0.1795 0.1746 0.1699 0.1654 0.1610 0.1567 0.1526 0.1486 0.1447 0.1410]';
G=zeros(20,20); %initialize the G
% construct the G matrix
for j=1:length(y)
    for i=1:length(x)
        G(i,j)=dx*x(j)*exp(-x(j)*y(i));
    end
end

e = logspace(-6,-3,100);
SN = zeros(length(e),1);%to initialize the solutional norm
RN = zeros(length(e),1);%to initialize the residual norm
%calculate the 2-norm of the two parameters
for i=1:length(e)
    ml = (G'*G+e(i)^2*eye(N))\(G'*d);
    dl = G*ml;
    SN(i) = norm(ml);
    RN(i) = norm(dl-d);
end
corner = 65;
SNc = SN(corner);
RNc = RN(corner);
e1 = e(corner);
figure(2)
plot(RN,SN,'r.',RNc,SNc,'bo','MarkerSize',10),legend('L-Curve','corner point')
% plot(RN,SN,'r.')
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
xlabel('Residual norm');
ylabel('Solutional norm');
title('L-Curve');

%plot the GCV curve
e2 = linspace(1e-6,1e-4,1000);%has the better plot
V = zeros(length(e2),1);
for i=1:length(e2)
    m2 = (G'*G+e2(i)^2*eye(N))\(G'*d);
    d2 = G*m2;
    V(i) = N*((d-d2)'*(d-d2)/(trace(eye(N)-G*((G'*G+e2(i)^2*eye(N))\G')))^2);
end
[minimumV,point] = min(V);
egcv = e2(point);
figure(3);
loglog(e2,V,'r.',e2(point),minimumV,'ko','MarkerSize',10);
legend('GCV curve','Minimum Point');
xlabel('\epsilon');
ylabel('V(\epsilon^2)');
title('GCV-Curve');

%exammine the GCV and L-curve method
mL = (G'*G + e1^2 * eye(N))\(G'*d);
dL = G*mL;
mgcv = (G'*G + egcv^2 * eye(N))\(G'*d);
dgcv = G * mgcv;

figure(4);
subplot(2,1,1)
plot(p,mL,'ro',p,mgcv,'g*'),xlabel('i'),ylabel('m(i)')
legend('L-curve method','GCV method')
subplot(2,1,2)
plot(p,d,'k.',p,dL,'r*',p,dgcv,'go'),xlabel('i'),ylabel('d')
legend('original data','L-curve','GCV')
