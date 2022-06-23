%% GRADIENT DESCENT METHOD
clear;
clc;
x = (0:1:10)';
N = length(x);
mgo = [0.25,0.25]';%give the initial guess of alpha and beta
ygo = exp(-mgo(1)*x).*sin(mgo(2)*x);
dobs = [0,0.4560,0.7614,0.8586,0.7445,0.4661,0.1045,-0.2472,-0.5073,-0.6233,-0.5816]';
Ego = (ygo-dobs)'*(ygo-dobs);

% calculate the Gp
dydmo = zeros(N,2);
dydmo(:,1) = -mgo(1)*exp(-mgo(1)*x).*sin(mgo(2)*x);
dydmo(:,2) = mgo(2)*exp(-mgo(1)*x).*cos(mgo(2)*x);

dEdmo = 2*dydmo'*(ygo-dobs);% calculate the vector b

alpha = 0.05;c1 = 0.0001;tau = 0.5;Niter = 500;

for k = 1:Niter
    v = -dEdmo/sqrt(dEdmo'*dEdmo);
    for kk = 1:10
        mg = mgo + alpha*v; % add the value of the step number to update the mg
        yg = exp(-mg(1)*x).*sin(mg(2)*x);
        Eg = (yg - dobs)'*(yg - dobs); %???? 
        dydm = zeros(N,2);
        dydm(:,1) = -mg(1)*exp(-mg(1)*x).*sin(mg(2)*x);
        dydm(:,2) = mg(2)*exp(-mg(1)*x).*cos(mg(2)*x);
        dEdm = 2*dydm'*(yg - dobs);
        if(Eg <= (Ego + c1*alpha*v'*dEdmo))
            break;
        end
        alpha = tau*alpha;
    end
    Dmg = sqrt((mg - mgo)'*(mg - mgo));
    mgo = mg;ygo = yg; Ego = Eg;
    dydmo = dydm;dEdmo = dEdm;
    if(Dmg < 1.0e-6)
        break;
    end
end
fprintf("The optimal parameters: alpha = %.2f, beta = %.2f\n",mg(1),mg(2));
figure(1)
plot(x,dobs,'ro',x,yg,'b.','LineWidth',2,'MarkerSize',10),xlabel('x','FontSize',15),ylabel('d','FontSize',15),legend('data-obs','data-prediction')


