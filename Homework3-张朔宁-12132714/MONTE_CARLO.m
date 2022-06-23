%% MONTE CARLO METHOD
clear;
clc;
m=[0:0.01:1;0:0.01:1]';
x=(0:10)';
dobs=[0,0.456,0.7614,0.8596,0.7445,0.4661,0.1045,-0.2472,-0.5073,-0.6233,-0.5816]';
E=zeros(101,101);
for i=1:101
    for j=1:101
        yg=exp(-m(i,1)*x).*sin(m(j,2)*x);
        Eg1=abs((yg-dobs)'*(yg-dobs));
        E(i,j)=Eg1;
    end
end
figure(1)
imagesc(0:0.01:1,0:0.01:1,E');
colorbar;
mg = [0.25,0.25]';
dg = exp(-mg(1)*x).*sin(mg(2)*x);
Eg = (dobs - dg)'*(dobs - dg);
hold on
plot(mg(1),mg(2),'ko',"LineWidth",3);

%randomly generate pairs of model parameters and check if they further
%minimize the error
Niter = 10000;
ma = zeros(2,1);
choosek = [100,200,300,1000,2000,3000,4000,5000,6000,7000,8000,9000];
for k = 1:Niter
    ma(1) = random('Unif',0,1);
    ma(2) = random('Unif',0,1);
    %compute the error
    da = exp(-ma(1)*x).*sin(ma(2)*x);
    Ea = (dobs - da)'*(dobs-da);

    %adopt it if it is better
    if (Ea < Eg)
        mg = ma;
        Eg = Ea;
    end

    %save history
    Ehis = zeros(length(k),1);
    m1his = zeros(length(k),1);
    m2his = zeros(length(k),1);
    Ehis(1+k) = Eg;
    m1his(1+k) = mg(1);
    m2his(1+k) = mg(2);
    hold on
    plot([m1his(1+k-1),m1his(1+k)],[m2his(1+k-1),m2his(1+k)],'r','linewidth',2)


    fprintf('已经执行%d次\n',k);

end
hold on
plot(mg(1),mg(2),'go','LineWidth',3);
xlabel('\alpha','FontSize',15),ylabel('\beta','FontSize',15)
figure(2)
dg1 = exp(-mg(1)*x).*sin(mg(2)*x);
plot(x,dobs,'r*',x,dg1,'bo','LineWidth',2,'MarkerSize',10),xlabel('x','FontSize',15),ylabel('d','FontSize',15),legend('data-obs','data-prediction')