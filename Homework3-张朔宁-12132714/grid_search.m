% grid search method
clear;
clc;
m=[0:0.02:2;0:0.02:2]';
x=(0:10)';
dobs=[0,0.456,0.7614,0.8596,0.7445,0.4661,0.1045,-0.2472,-0.5073,-0.6233,-0.5816]';
E=zeros(101,101);
for i=1:101
    for j=1:101
        yg=exp(-m(i,1)*x).*sin(m(j,2)*x);
        Eg=(yg-dobs)'*(yg-dobs);
        E(i,j)=Eg;
    end
end
imagesc(0:0.02:2,0:0.02:2,E');
colorbar
xlabel('\alpha');
ylabel('\beta');

mgo=zeros(2,1);
N=length(x);
ma=[0:0.5:2;0:0.5:2]';% set different initial guess
n=length(ma);

for i=1:n
    for j=1:n  
        mgo(1)=ma(i,1);
        mgo(2)=ma(j,1);
        hold on;
        plot(mgo(1),mgo(2),'ko','MarkerSize',10);
        ygo=exp(-mgo(1)*x).*sin(mgo(2)*x);
        Ego=(ygo-dobs)'*(ygo-dobs);
        dydmo=zeros(N,2);
        dydmo(:,1)=-x.*exp(-mgo(1)*x).*sin(mgo(2)*x);
        dydmo(:,2)=exp(-mgo(1)*x).*cos(mgo(2)*x).*x;
        dEdmo=2*dydmo'*(ygo-dobs);

        alpha = 0.05;c1=0.0001;tau=0.5;Niter=500;
        for k=1:Niter
            v=-dEdmo/sqrt(dEdmo'*dEdmo);
            for kk=1:10
            mg=mgo+alpha*v;
            hold on
            plot(mg(1),mg(2),'r*','LineWidth',1);
            yg=exp(-mg(1)*x).*sin(mg(2)*x);
            Eg=(yg-dobs)'*(yg-dobs);
            dydm=zeros(N,2);
            dydm(:,1)=-x.*exp(-mg(1)*x).*sin(mg(2)*x);
            dydm(:,2)=exp(-mg(1)*x).*cos(mg(2)*x).*x;
            dEdm=2*dydm'*(yg-dobs);
            if((Eg<=(Ego+c1*alpha*v'*dEdmo)))
                break;
            end
            alpha=tau*alpha;
            end
        Dmg=sqrt((mg-mgo)'*(mg-mgo));
        mgo=mg;
        ygo=yg;
        Ego=Eg;
        dydmo=dydm;
        dEdmo=dEdm;
        if(Dmg<1.0e-6)
            break;  
        end
        end
    end
end
hold on
plot(mg(1),mg(2),'go','MarkerSize',10,'LineWidth',2)




