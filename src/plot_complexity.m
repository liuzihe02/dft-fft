% First, modify your complexity plotting function:
function plot_complexity(N_values, mean_runtime_data, std_runtime_data, alg_names, results)
    % Define consistent colors to use across all plots
    plot_colors = get(gca, 'ColorOrder');
    
    % Simple marker options
    markers = {'o', 's', '^'};
    
    % Create figure
    fig = figure;
    hold on;
    
    % Create a fine grid for smooth curve plotting
    N_fine = logspace(log10(min(N_values)), log10(max(N_values)), 100)';
    
    % Create legend entries
    legend_entries = cell(1, 2 * length(alg_names));
    
    % Loop over algorithms
    for a = 1:length(alg_names)
        % Get color for this algorithm (consistent across plots)
        alg_color = plot_colors(mod(a-1, size(plot_colors, 1))+1, :);
        
        % Plot data points with error bars with specified color
        h = errorbar(N_values, mean_runtime_data{a}, std_runtime_data{a}, ...
            'Marker', markers{a}, ...
            'MarkerSize', 8, ...
            'LineStyle', 'none', ...
            'Color', alg_color, ...
            'MarkerFaceColor', alg_color);
        
        % Create legend entry for data points
        legend_entries{2*a-1} = [alg_names{a}, ' data'];
        
        % Retrieve the best-fit model info
        best_model_idx = results.best_fits(a, 1);
        best_params = results.fitted_params{a, best_model_idx};
        best_model_name = results.model_names{best_model_idx};
        model_func = results.models{best_model_idx};
        
        % Generate the fitted curve
        fitted_curve = model_func(best_params, N_fine);
        
        % Plot fitted curve with the same color
        loglog(N_fine, fitted_curve, ...
               'LineWidth', 1.5, ...
               'Color', alg_color);
               
        % Create legend entry for fit line
        legend_entries{2*a} = [alg_names{a}, ' fit (', best_model_name, ')'];
    end
    
    % Set plot properties
    xlabel('Sequence Length (N)');
    ylabel('Runtime (seconds)');
    title('Algorithm Complexity with Best-Fit Models');
    
    % Apply the legend
    legend(legend_entries, 'Location', 'NorthWest');
    
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    % Save the figure
    print(fig, 'complexity_fitted.png', '-dpng', '-r300', '-painters');
    fprintf('Complexity fit plot saved\n');
end