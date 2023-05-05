close all
clear all
x = linspace(0,1,100);
% - arcsine pdf
% y = 1./(pi.*sqrt(x.*(1-x)));
% plot(x,y)

% - ZEROTH MOMENT
% - arcsine cdf on wikipedia
y = (2/pi).*asin(sqrt(x));
plot(x,y,'r--')
% - arcsine cdf computed with wolfram alpha
hold on
y = 1-(2*asin(sqrt(1 - x)))/pi;
plot(x,y,'b*')
hold off

% - FIRST MOMENT
% - wolfram alpha: integrate x/(pi sqrt(x(1-x))). Need to add 1/2 at the
% end to get the first moment that at x = 1 gives the E[x], for x inn [0,1]
% which is indeed 0.5 as the distribution is symmetric

y = 1/2-(sqrt(x-x.^2)+asin(sqrt(1-x)))./pi;
plot(x,y)
