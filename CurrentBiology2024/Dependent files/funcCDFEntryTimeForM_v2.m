function cdf = funcCDFEntryTimeForM_v2(params, inputs)
%(TN 2022/09/13) Calculating the CDF by a given time t for the entry dynamics of co-infecting phages
%   params = [eta, k, tau]
%   inputs = an n-by-3 matrix, with each row being [m, i, t], where m = #injecting phages, i = entry order, and t = timepoint of interest
%   This is a vectorized function, returning a column vector.
%       eta, dimensionless
%       k, per unit of time
%       tau, unit of time
%       m, dimensionless
%       i, dimensionless
%       t, unit of time

%Argument validation
if numel(params) ~= 3 | size(inputs, 2) ~= 3
    error('Dimensions of params or inputs are invalid.');
    cdf = NaN;
    return
end

%Parameter specifications
eta = params(1);
k = params(2);
tau = params(3);
nCells = size(inputs, 1);

cdf = zeros(nCells, 1);
for index_cell = 1:nCells
    m = inputs(index_cell, 1);
    i = inputs(index_cell, 2);
    t = inputs(index_cell, 3);

    %Boundary case
    if i == 0 %The zero-th entry has always happened
        cdf(index_cell) = 1;
        return
    end

    %Stop the function if the inputs are invalid
    if t < tau | m < i
        cdf(index_cell) = 0;
        return
    end

    %Initialization
    jVector = 1:1:i;

    %Waiting time terms
    waitingTerms = zeros(1, numel(jVector));
    for index_j = 1:numel(jVector)
        j = jVector(index_j);

        waitingTerms(index_j) = (-1)^j * nchoosek(i-1, j-1)...
            * (1 - exp(-(m-(j-1)) * k * (t-tau)))/(m - (j-1));
    end
    waitingSum = sum(waitingTerms);

    %Final calculation
    cdf(index_cell) = nchoosek(m, i) * i * (-1)^i * waitingSum;
end
end