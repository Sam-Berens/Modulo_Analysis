function [x] = getLearnOffset(a,b)
% a is the learning onset for a given question type
% b is the associated learning rate
% x is the learning offset for that same question type

%function that takes a value of kappa (conc param of vonMis)
%and gives the probability of the target bin being chosen
 f = @(k)exp(k)/sum(exp(k*cos((0:5)*pi/3)));

 %estimate the value of kappa required to get a probability of 0.99
 % for the target bin
cK = fmincon(@(k)(f(k)-0.99).^2,1);

%convert your cK (estimated kappa) value to the resultant vector,R, for the subject's response
%as mapped on a unit circle?

% %% For any value of cK computation would be:
% if cK < 7e-3
%     cR = cK;
% elseif (cK < 3.648) && (cK >= 1.2516)
%     cR = (-sqrt((10000*cK^2 - 19800*cK + 33709)) +100*cK + 179) / 278;
% elseif cK < 1.2516
%     cR = roots([5,0,6,0,12,(-6*cK)]);
%     cR = abs(cR(end,1));
% else
%     cR = roots([1,-4,3,(-1/cK)]);
%     cR = abs(cR(2));
% end
%we only need this bit because we know that kappa is 10.5736 for a P(corr) of 0.99
cR = roots([1,-4,3,(-1/cK)]); 
cR = abs(cR(2));


%to get from our resultant vector to the subject's learning offset, x, 
%we use the following formula, which is a rearraned
%version of our softplus learning model...

x = log(exp(atanh(cR)/b)-1) + a;
end 


