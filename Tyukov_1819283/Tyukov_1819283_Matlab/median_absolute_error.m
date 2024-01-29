function median_absolute_error = median_absolute_error(simulated_x, simulated_y, measured_x, measured_y, tolerance)
    errors = [];

    % Remove possible zero values in measured data
    measured_x = measured_x(measured_y ~= 0);
    measured_y = measured_y(measured_y ~= 0);

    for i = 1:length(measured_x)
        ind = find(abs(simulated_x - measured_x(i)) <= tolerance);

        % Loop through all matching simulated_x values
        for j = 1:length(ind)
            error = 100*abs((simulated_y(ind(j)) - measured_y(i))/measured_y(i));
            errors = [errors, error]; % Append to errors array
        end
    end
    % Calculate median absolute error
    median_absolute_error = median(errors);
end