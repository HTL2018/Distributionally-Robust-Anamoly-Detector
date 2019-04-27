function s = uniform_random_z (sigma,alpha,n,false_alarm)
beta=(6*sigma-alpha)/(alpha*false_alarm);

%%% sigma should be bigger than alpha/6 otherwise you can not cover all of
%%% 0<x<alpha with uniform distribution and have something extra for beyond alpha.

b=((sqrt(1+4*beta)-1)/2)^2;
s=zeros(1,n);
for i=1:n
    xrand=rand;
    if xrand<1-false_alarm
        s(i)=alpha*xrand/(1-false_alarm);
    else
        s(i)=alpha*(b-1)*(  (xrand/false_alarm) - ( ((1-false_alarm)/false_alarm)  -...
            (1/(b-1)))   );
    end
end
% histogram(s,10,'Normalization','pdf','DisplayStyle','stairs');
end