function pdf = funcPDFEntryTime_v2(params, inputs)
%(TN 2022/09/13) Calculating the PDF of a given time t for the entry dynamics of co-infecting phages
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
    pdf = NaN;
    return
end

%Parameter specifications
eta = params(1);
k = params(2);
tau = params(3);
nCells = size(inputs, 1);

pdf = zeros(nCells, 1);
for index_cell = 1:nCells
    n = inputs(index_cell, 1);
    i = inputs(index_cell, 2);
    t = inputs(index_cell, 3);

    %Boundary case
    if i == 0 %The zero-th entry has always happened
        if t == 0
            pdf(index_cell) = 1;
        else
            pdf(index_cell) = 0;
        end
        continue
    end

    %Stop the function if the inputs are invalid
    if t < tau | n < i
        pdf(index_cell) = 0;
        continue
    end

    %Initialization
    mVector = i:1:n;

    %Normalization term
    normalizationTerms = zeros(1, numel(mVector));
    for index_m = 1:numel(mVector)
        m = mVector(index_m);
        normalizationTerms(index_m) = nchoosek(n, m) * eta^m * (1-eta)^(n-m);
    end
    normalizationSum = sum(normalizationTerms);

    %Waiting time terms
    waitingTerms = zeros(1, numel(mVector));
    for index_m = 1:numel(mVector)
        m = mVector(index_m);

        waitingTerms(index_m) = nchoosek(n, m) * eta^m * (1-eta)^(n-m)...
            * nchoosek(m, i) * i * k * exp(-(m - (i-1)) * k * (t-tau))...
            * (1 - exp(-k * (t-tau)))^(i-1);
    end
    waitingSum = sum(waitingTerms);

    %Final calculation
    pdf(index_cell) = normalizationSum^(-1) * waitingSum;
end
end
