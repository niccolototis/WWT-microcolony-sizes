function xres = analyticalMoments(tspan,m10,GrowthRate,BreakageRate)
% Analytical defivation of the fist moment for (cell size cell growth and division model)
global scenario

if isequal(GrowthRate,'constant') && isequal(BreakageRate,'constant') 
 scenario = 1; % 1) f(x) = g, b(x) = k
 [~,xres] = ode15s(@odeFun, tspan, m10);
else
if isequal(GrowthRate,'linear') && isequal(BreakageRate,'constant')
 scenario = 2;% 2) f(x) = g*x, b(x) = k    
 [~,xres] = ode15s(@odeFun, tspan, m10);
else 
if isequal(GrowthRate,'constant') && isequal(BreakageRate,'linear')
 scenario = 3;% 3) f(x) = g, b(x) = k*x
 [~,xres] = ode15s(@odeFun, tspan, m10);
else
if isequal(GrowthRate,'linear') && isequal(BreakageRate,'linear')
 scenario = 4;% 4) f(x) = g*x, b(x) = k*x
 [~,xres] = ode15s(@odeFun, tspan, m10);
else
if isequal(GrowthRate,'linear') && isequal(BreakageRate,'zero')
 scenario = 5;% 4) f(x) = g*x, b(x) = 0
 [~,xres] = ode15s(@odeFun, tspan, m10);
end
end
end
end
end
return

function dm1 = odeFun(t,m1)
global scenario g k
switch scenario
 case 1
    dm1 = -k*m1+g;
 case 2
    dm1 = -k*m1;
 case 3
    dm1 = -k*m1^2+g;
 case 4
    dm1 = -k*m1^2+g*m1;
 case 5
    dm1 = +g*m1;
end
return
