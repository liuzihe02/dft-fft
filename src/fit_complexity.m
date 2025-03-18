function results = fit_complexity(N_values, runtime_data, alg_names)
    % COMPLEXITY Analyzes time complexity of algorithms
    %   Fits various complexity models to runtime data and determines best fit
    %
    %   Inputs:
    %     N_values - Array of input sizes
    %     runtime_data - Cell array of runtime measurements for each algorithm
    %     alg_names - Cell array of algorithm names
    %
    %   Output:
    %     results - Structure with analysis results
    
    % Define the models to test as anonymous functions
    models = {
        @(b,x) b(1)*x + b(2),                  % O(n)
        @(b,x) b(1)*x.^2 + b(2),               % O(n²)
        @(b,x) b(1)*x.^3 + b(2),               % O(n³)
        @(b,x) b(1)*log2(x) + b(2),            % O(log n)
        @(b,x) b(1)*x.*log2(x) + b(2)          % O(n log n)
    };
    model_names = {'O(n)', 'O(n²)', 'O(n³)', 'O(log n)', 'O(n log n)'};
    
    % Initialize results
    num_algs = length(alg_names);
    num_models = length(models);
    r_squared = zeros(num_algs, num_models);
    best_fits = zeros(num_algs, 2); % [model_index, R²]
    fitted_params = cell(num_algs, num_models);
    
    % Set optimization options for lsqcurvefit (suppress display)
    options = optimoptions('lsqcurvefit','Display','off');
    
    % Fit each model to each algorithm's data
    fprintf('\n--- Complexity Analysis (R² values) ---\n');
    fprintf('%-12s', 'Algorithm');
    for m = 1:num_models
        fprintf('%-12s', model_names{m});
    end
    fprintf('\n');
    
    % Loop over algorithms
    for a = 1:num_algs
        fprintf('%-12s', alg_names{a});
        best_R2 = -Inf;
        best_model = 0;
        
        % Ensure N_values and runtime_data are column vectors
        xdata = N_values(:);
        ydata = runtime_data{a}(:);
        
        for m = 1:num_models
            % Initial guess for parameters [slope, intercept]
            b0 = [1e-6, 0];
            % Fit the model using lsqcurvefit
            try
                b_fit = lsqcurvefit(models{m}, b0, xdata, ydata, [], [], options);
            catch ME
                warning('Fit failed for algorithm %s with model %s: %s', ...
                    alg_names{a}, model_names{m}, ME.message);
                b_fit = [NaN, NaN];
            end
            fitted_params{a, m} = b_fit;
            
            % Compute the fitted values and R² manually
            yfit = models{m}(b_fit, xdata);
            SS_res = sum((ydata - yfit).^2);
            SS_tot = sum((ydata - mean(ydata)).^2);
            if SS_tot == 0
                rsq = 1; % Perfect fit if ydata is constant
            else
                rsq = 1 - SS_res/SS_tot;
            end
            r_squared(a, m) = rsq;
            fprintf('%-12.4f', rsq);
            
            % Track the best model
            if rsq > best_R2
                best_R2 = rsq;
                best_model = m;
            end
        end
        
        best_fits(a,:) = [best_model, best_R2];
        fprintf('\n');
    end
    
    % Report best models
    fprintf('\n--- Best Complexity Model for Each Algorithm ---\n');
    for a = 1:num_algs
        fprintf('%s: Best model is %s (R² = %.4f)\n', alg_names{a}, ...
            model_names{best_fits(a,1)}, best_fits(a,2));
    end
    
    % Prepare output results
    results = struct();
    results.r_squared = r_squared;
    results.best_fits = best_fits;
    results.model_names = model_names;
    results.fitted_params = fitted_params;
    results.models = models;
end
