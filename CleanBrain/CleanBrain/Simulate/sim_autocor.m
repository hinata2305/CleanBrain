function [y] = sim_autocor(a,b,sigma_par,time)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sigma=sigma_par;

y0 = a / (1 - b); %eg initialize to unconditional mean of stationary time series

y = zeros(time,1);
y(1) = a + b * y0 + randn() * sigma;

for t = 2:time
    y(t) = a + b * y(t-1) + randn() * sigma;
end


end

