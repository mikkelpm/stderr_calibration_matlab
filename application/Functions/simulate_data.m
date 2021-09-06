function [change, price] = simulate_data(moments, quant_grid, quantiles, n)

    % Simulate data on price changes

    change = (rand(n,1)<moments(1)); % Simulate price change indicators
    price = pchip(quant_grid, quantiles, rand(n,1)); % Simulate price changes (interpolate quantile function outside grid)

end