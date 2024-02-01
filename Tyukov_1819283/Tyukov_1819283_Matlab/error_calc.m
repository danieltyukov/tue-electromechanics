function error = error_calc(simulated_x, simulated_y, measured_x, measured_y, tolerance)
    % This function calculates the median absolute error between
    % simulations and measured data using a specified method for calculating errors.
    % Input: 
    %   simulated_x - x values of the simulated data
    %   simulated_y - y values of the simulated data
    %   measured_x - x values of the measured data
    %   measured_y - y values of the measured data
    %   tolerance - The tolerance within which x values are considered a match
    %   errorMethod - Function handle to the method used for calculating errors
    % Output: 
    %   median_absolute_error - The median absolute error between the datasets

    % Preallocate the errors array for maximum possible size
    error = zeros(1, size(simulated_x, 2));
    error_index = 1;

    measured_x = measured_x((measured_y ~= 0) & (~isnan(measured_y)));
    measured_y = measured_y((measured_y ~= 0) & (~isnan(measured_y)));
    % Loop through the measured array
    for i = 1:length(measured_x)
        % Find simulated_x values within the tolerance of the current measured_x value
        x_index = find(abs(simulated_x - measured_x(i)) <= tolerance);
        
        % Loop through x_index to calculate errors for all matching simulated_x values
        for j = 1:length(x_index)
            % Calculate error using the specified method
            err = 100*abs((simulated_y(x_index(j)) - measured_y(i))/measured_y(i));
            if (err > 100 | isnan(err))
                fprintf("%f at %f, %f at %f\n", simulated_y(x_index(j)), simulated_x(x_index(j)), measured_y(i), measured_x(i));
            end
            error(error_index) = err; % Store error
            error_index = error_index + 1;
        end
    end

    % Remove unused preallocated space
    error = error(1:error_index-1);
end
