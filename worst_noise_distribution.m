
function sample = worst_noise_distribution(input_param)
% Computes the worst case noise distribution
    
    n          = size(input_param.Sigma_noise,2); % dimension
    x_bar      = input_param.mu_noise;
    S          = input_param.Sigma_noise;
    noise_type = input_param.noise_type;
    A          = inv(S);
    
    out_param = compute_noise_threshold(input_param);
    z      = out_param.z;
    Z      = out_param.Z;
    lambda = out_param.lambda;
    c      = -out_param.noise_threshold; 
    x0     = 1/(1-lambda)*(x_bar - lambda*z/lambda);
    S0     = 1/(1-lambda)*(S - lambda*Z/lambda);

    %% Construction of Random variable
    Z = Z/lambda;
    z = z/lambda;

    [V,E] = eig(Z - z*z');
    lam   = z'*z + c;
    W     = V*sqrt(E);
    w1    = W(:,1);
    w2    = W(:,2);

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
    
    
    figure;
    % Plot x_bar
    f(1) = plot(x_bar(1),x_bar(2),'O', 'Color', 'm', 'MarkerSize',15, 'MarkerFaceColor', 'm');
    hold on;    
    % Unit Circle
    syms xx yy 
    ellipse_function  = [xx yy]*eye(n)*[xx;yy] + c;
    ellipse_function_handle = matlabFunction(ellipse_function);
    f(2) = fimplicit(ellipse_function_handle,[-20 20 -20 20],'r','lineWidth',4);
    hold on;
    f(3) = plot(x0(1), x0(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
    text(x0(1), x0(2)+0.5, num2str(1-lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
    plot(v1(1), v1(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
    text(v1(1)-0.25, v1(2)+0.5, num2str(alpha1*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
    plot(v2(1), v2(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
    text(v2(1), v2(2)+0.5, num2str(alpha2*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
    plot(v1r(1), v1r(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
    text(v1r(1)-0.25, v1r(2)+0.5, num2str(alpha1r*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
    plot(v2r(1), v2r(2), 'p', 'Color', 'k', 'MarkerSize',15, 'MarkerFaceColor', 'k');
    text(v2r(1), v2r(2)+0.5, num2str(alpha2r*lambda),'Color','blue','FontSize',30, 'Interpreter', 'latex');
    if noise_type == 1
        legend(f(1:3),'$\bar{x}$','Circle with radius $\sqrt{\bar{w}^{*}}$','Discrete Supports', 'Interpreter', 'latex');
        xlabel('$w_1$', 'Interpreter', 'latex');
        ylabel('$w_2$', 'Interpreter', 'latex');
    else
        legend(f(1:3),'$\bar{x}$','Circle with radius $\sqrt{\bar{v}^{*}}$','Discrete Supports', 'Interpreter', 'latex');
        xlabel('$v_1$', 'Interpreter', 'latex');
        ylabel('$v_2$', 'Interpreter', 'latex');
    end
    a = findobj(gcf, 'type', 'axes');
    h = findobj(gcf, 'type', 'line');
    
    set(h, 'linewidth', 4);
    set(a, 'linewidth', 4);
    set(a, 'FontSize', 30);
    set(gca,'TickLabelInterpreter','latex');
    axis square;
    hold off


    %% discrete distribution: x0 with probability 1-lambda, v2 with probability alpha2*lambda, v2r with probability alpha2r*lambda
    sample.probabilities = [1-lambda alpha1*lambda alpha1r*lambda alpha2*lambda alpha2r*lambda];
    sample.support       = [x0 v1 v1r v2 v2r];
    sample.threshold     = out_param.noise_threshold;
end