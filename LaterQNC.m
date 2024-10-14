%%  1. Get the data
data = later_getData([], [], 0.2);
RTs = data{1};
clear data

%%  2. Define the objective function
laterErrFcn = @(fits) -sum(log(normpdf(RTs, fits(1), fits(2))));

%%  3. Define initial conditions
lowerBounds = [0.001 0.001];
upperBounds = [1000 1000]; 
% Initial values: use mean and standard deviation of the reciprocal of RTs
initialValues = [mean(1 ./ RTs) * 1.5, std(1 ./ RTs)];  % [muR, deltaS]

% Step 4: Run the Fits 
opts = optimoptions(@fmincon, 'Algorithm', 'active-set', 'MaxIter', 3000, 'MaxFunEvals', 3000); 
problem = createOptimProblem('fmincon', 'objective', laterErrFcn, 'x0', initialValues, ... 'lb', lowerBounds, 'ub', upperBounds, 'options', opts); 
gs = GlobalSearch;
[fits, nllk] = run(gs, problem); % Optimize and get fitted parameters

%%  5. Evaluate the fits
histogram(RTs, 'Normalization', 'pdf');
hold on;

muR = fits(1);
deltaS = fits(2);
x = linspace(min(RTs), max(RTs), 100);
predictedPDF = normpdf(x, muR, deltaS);
plot(x, predictedPDF, 'r-', 'LineWidth', 2);

legend('Observed RT', 'Model Prediction');
xlabel('Reaction Time (RT)');
ylabel('Probability Density');
title('Model Fit to Reaction Time Data');
hold off;

