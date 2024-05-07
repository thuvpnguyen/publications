function inside = funcInsideOutside_v2(params, inputs)
%(TN 2022/09/13) Calculating the expected inside MOI for the entry dynamics of co-infecting phages
%   params = [eta, k, tau]
%   inputs = an N-by-2 matrix, with each row being [n, t], where n = outside MOI, and t = timepoint of interest
%   This is a vectorized function, returning a column vector.
%       eta, dimensionless
%       k, per unit of time
%       tau, unit of time
%       n, dimensionless
%       t, unit of time

%Argument validation
if numel(params) ~= 3 | size(inputs, 2) ~= 2
    error('Dimensions of params or inputs are invalid.');
    inside = NaN;
    return
end

%Parameter specifications
eta = params(1);
k = params(2);
tau = params(3);
nCells = size(inputs, 1);

inside = zeros(nCells, 1);
for index_cell = 1:nCells
    n = inputs(index_cell, 1);
    t = inputs(index_cell, 2);

    %Initialization
    iVector = 0:1:n; %For all possible values of inside MOI

    iteratingTerms = zeros(1, numel(iVector));
    for index_i = 1:numel(iVector)
        i = iVector(index_i);

        iteratingTerms(index_i) = i * funcProbInsideOutside_v2([eta, k, tau], [n, i, t]);
    end

    %Final calculation
    inside(index_cell) = sum(iteratingTerms);
end
end