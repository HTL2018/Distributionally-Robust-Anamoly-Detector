%% Upper Bound SDP
clear all; clc; close all;

n     = 2; % dimension
% x_bar = [0.2 0.3]';
x_bar  = [0 0]';

b     = [0 1]';
c     = -8.7184;
% S     = A;
% S     = [0.2 0.06; 0.06 0.11];
% S = [0.045  -0.011               % System Noise Covariance Matrix
%      -0.011 0.02];
S = [2.086 0.134
     0.134 2.230];              
A = eye(size(S,1));

attack_tol       = 0.00001;           
attack_high      = 1000;
attack_low       = 0;    
attack_threshold = 100;
iter_counter     = 1;       
alarm_rate       = 0.05;
    
% Get Attack_Threshold by solving below SDP using bisection algorithm
while (attack_high - attack_low > attack_tol) 
    
    attack_threshold = (attack_high + attack_low)/2   

    cvx_begin sdp quiet
        variable Z(n,n) symmetric
        variable z(n,1)
        variable lambda

        minimize (1 - lambda)
        subject to
            [Z z; z' lambda] <= [S x_bar; x_bar' 1];
            [Z z; z' lambda] >= 0;
            trace(A*Z) - attack_threshold*lambda >= 0;
    cvx_end

    if(cvx_optval > 1 - alarm_rate)
        attack_high = attack_threshold;            
    elseif(cvx_optval < 1 - alarm_rate)
        attack_low  = attack_threshold;             
    end

    threshold(iter_counter) = attack_threshold;
    iter_counter            = iter_counter + 1;

end        

x0 = 1/(1-lambda)*(x_bar - lambda*z/lambda);
S0 = 1/(1-lambda)*(S - lambda*Z/lambda);
 
%% 
Z = Z/lambda;
z = z/lambda;
 
[V,E] = eig(Z - z*z');
lam = z'*z + c;
W = V*sqrt(E);
w1 = W(:,1);
w2 = W(:,2);
 
mu1 = w1'*w1;
mu2 = w2'*w2;
 
beta1 = roots([mu1 2*w1'*z lam]);
beta2 = roots([mu2 2*w2'*z lam]);
 
v1  = z + beta1(1)*w1;
v1r = z + beta1(2)*w1;
v2  = z + beta2(1)*w2;
v2r = z + beta2(2)*w2;
 
alpha1  = mu1/(1-beta1(1)/beta1(2))/(mu1+mu2);
alpha2  = mu2/(1-beta2(1)/beta2(2))/(mu1+mu2);
alpha1r = -alpha1*beta1(1)/beta1(2);
alpha2r = -alpha2*beta2(1)/beta2(2);
 
alpha1*v1 + alpha1r*v1r + alpha2*v2 + alpha2r*v2r;
alpha1*v1*v1' + alpha2*v2*v2' + alpha1r*v1r*v1r' + alpha2r*v2r*v2r';

%% Plot The Figures
 
figure;
% Plot x_bar
f(1) = plot(x_bar(1),x_bar(2),'O', 'Color', 'm', 'MarkerSize',15, 'MarkerFaceColor', 'm');
hold on;
% S^{-1} Ellipse
syms xx yy 
ellipse_function  = [xx yy]*inv(S)*[xx;yy] -1;
ellipse_function_handle = matlabFunction(ellipse_function);
f(2) = fimplicit(ellipse_function_handle,[-3 3 -3 3],'b','lineWidth',4);
hold on;
% Unit Circle
syms xx yy 
ellipse_function  = [xx yy]*eye(n)*[xx;yy] + c;
ellipse_function_handle = matlabFunction(ellipse_function);
f(3) = fimplicit(ellipse_function_handle,[-3 3 -3 3],'r','lineWidth',4);
hold on;
f(4) = plot(x0(1), x0(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
text(x0(1), x0(2)+0.25, num2str(1-lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
plot(v1(1), v1(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
text(v1(1)-0.25, v1(2)+0.25, num2str(alpha1*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
plot(v2(1), v2(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
text(v2(1), v2(2)+0.25, num2str(alpha2*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
plot(v1r(1), v1r(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
text(v1r(1)-0.25, v1r(2)+0.25, num2str(alpha1r*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
plot(v2r(1), v2r(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
text(v2r(1), v2r(2)+0.25, num2str(alpha2r*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
legend(f(1:4),'$\bar{x}$','$S^{-1}$ Ellipse','Circle with radius $\sqrt{\alpha^{*}}$','Discrete Supports', 'Interpreter', 'latex')
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 24);
set(gca,'TickLabelInterpreter','latex')
hold off


%% discrete distribution: x0 with probability 1-lambda, v2 with probability alpha2*lambda, v2r with probability alpha2r*lambda
p = [1-lambda alpha1*lambda alpha1r*lambda alpha2*lambda alpha2r*lambda];
X = [x0 v1 v1r v2 v2r];
 
sum(diag(p)*X')  % should be equal to x_bar
X*diag(p)*X'  % should be equal to S
 
Ns = 100000;
x = zeros(2, Ns);
 
% sample just to check
for i=1:Ns
    x(:,i) = X(:,sample(p));
end
 
% sample mean (should be close to x_bar with many samples)
mean(x')
 
% sample second moment (should be close to S with many samples)
1/Ns*x*x'

 

