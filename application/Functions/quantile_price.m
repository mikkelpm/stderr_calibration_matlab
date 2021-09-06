function quantiles = quantile_price(num_prod, y_bar, quant_grid)

    % Compute quantiles of price change distribution in Alvarez & Lippi (ECMA 2014) model

    w = @(y) (1-y.^2/y_bar).^((num_prod-3)/2)/(beta((num_prod-1)/2,0.5)*sqrt(y_bar)); % Density of price changes (Equation 13)
    quantiles = nan(size(quant_grid));
    
    for j=1:length(quant_grid) % For each quantile in grid...
        if j==1
            x0 = -sqrt(y_bar);
        else
            x0 = quantiles(j-1);
        end
        quantiles(j) = fzero(@(y) integral(w, -sqrt(y_bar), y)-quant_grid(j), [x0 sqrt(y_bar)+sqrt(eps)]); % Invert CDF
    end

end