clear;
clc;
%construct a vector z with 11 elements equally spaced 0.1
z = [1:0.1:2];
%randomly assigns elements of mtrue
m = [0.1,0.3,-0.2,0.5]';

%build the appropriate kernel G
G = zeros(11,4);
G(:,1) = ones(1,11)';
G(:,2) = z';
G(:,3) = z'.^2;
G(:,4) = z'.^3;

%create synthetic data with Gaussian random numbers
mean = 0;
sigmad = 0.05;
dobs = G*m + normrnd(mean,sigmad,11,1);

%solve the inverse problem by simple least squares
mest = (G'*G)\(G'*dobs);

dpre = G*mest;
figure(3)
clf;
x = [1:11]';
plot(x,dobs,x,dpre,'r-.','linewidth',2),xlabel("parameter number"),ylabel("value of d"),legend('dobs','dpre')

