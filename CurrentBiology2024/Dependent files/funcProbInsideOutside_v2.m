function prob = funcProbInsideOutside_v2(params, inputs)
%(TN 2022/09/13) Calculating the probability for a given inside MOI for the entry dynamics of co-infecting phages
%   params = [eta, k, tau]
%   inputs = an N-by-3 matrix, with each row being [n, i, t], where n = outside MOI, i = entry order, and t = timepoint of interest
%   This is a vectorized function, returning a column vector.
%       eta, dimensionless
%       k, per unit of time
%       tau, unit of time
%       n, dimensionless
%       i, dimensionless
%       t, unit of time

%Argument validation
if numel(params) ~= 3 | size(inputs, 2) ~= 3
    error('Dimensions of params or inputs are invalid.');
    prob = NaN;
    return
end

%Parameter specifications
eta = params(1);
k = params(2);
tau = params(3);
nCells = size(inputs, 1);

prob = zeros(nCells, 1);
for index_cell = 1:nCells
    n = inputs(index_cell, 1);
    i = inputs(index_cell, 2);
    t = inputs(index_cell, 3);

    %Initialization
    mVector = i:1:n; %Meaningless if #injecting phages < Inside MOI

    iteratingTerms = zeros(1, numel(mVector));
    for index_m = 1:numel(mVector)
        m = mVector(index_m);

        iteratingTerms(index_m) = nchoosek(n, m) * eta^m * (1-eta)^(n-m)...
            * (funcCDFEntryTimeForM_v2([eta, k, tau], [m, i, t])...
            - funcCDFEntryTimeForM_v2([eta, k, tau], [m, i+1, t]));
    end

    %Final calculation
    prob(index_cell) = sum(iteratingTerms);
end
end