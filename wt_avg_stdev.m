function [xavg_wt,dx_wt] = wt_avg_stdev(xavg,dx)

% 2020-08-16 Roy Ready. Calculates a weighted mean
% and standard deviation for a given vector of measurements and errors.
% Planned: calculate standard error and chi-square too.

    x = [xavg ; dx];

    sample_size = length(xavg);

    %calculate weights & product of weight and measurement
    wt_x = zeros(sample_size,1);

    xwt = zeros(sample_size,1);

    for i = 1:sample_size
        
        
        wt_x(i) = x(2,i)^(-2);
        
        xwt(i) = x(1,i)*wt_x(i);
        
    end

    %calculate weight sum

    wt_x_sum = sum(wt_x(:));

    xavg_wt = sum(xwt(:))/wt_x_sum;

    dx_wt = (sum(wt_x(:)))^(-1/2);



    
    

