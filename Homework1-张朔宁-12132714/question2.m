clear;
clc;

 m = rand(1,100);%randomly assigns masses mtrue in the range of 0-1kg
 % build the appropriate kernel G
 % G is a 100*100 matrix
 G1 = tril(ones(100)); 
 G = sparse(G1);
 % create synthetic observed data by adding Gaussian random vector
 sigmad = 0.01;
 dobs = G1*m' + normrnd(0,sigmad,100,1);
%solve the inverse problem by simple least squares
 mest = (G'*G)\(G'*dobs);
 % to calculate the variance of each of the estimated model parameters
 var = std2(mest);

sigmam = sqrt(var);
count = 0;
for j = (1:length(mest))
   a = mest(j);
   b = m(j);
   if abs(a-b) <= 2*sigmam
       count = count + 1;
   end
end
fprintf("the number of estimated model parameters that are within 2σ of their true value is %d\n",count)

%draw the picture
figure(2)
x = [1:100];
clf;
plot(x,m,x,mest,'r-','linewidth',2),title('Fitting result(the 2nd measurement)'),xlabel('parameter number'),ylabel('mass/kg'),legend('initial value','inversed value')