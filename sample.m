function z = sample(p)
% draws a sample from a discrete pdf

% form cdf
c = zeros(length(p),1);

for i=1:length(p)
    c(i) = ones(1,i)*p(1:i)';
end

z = find(rand <= c, 1);