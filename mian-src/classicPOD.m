% Proper Orthogonal Decomposition (POD) without using SVD
% X is the input data matrix n*N, n is the dimension, N is the snapshots
function [PODMs, PODEvals, SortPOD, PODReTMs] = classicPOD(X)
    % Input:
    % X - Data matrix (n x N) where n is the number of spatial points,
    %     N is the number of time snapshots.
    %
    % Output:
    % PODMs - POD spatial modes (n x N)
    % PODEvals - Eigenvalues from the covariance matrix
    % SortPOD - Sorted norms based on eigenvalues
    % PODReTMs - Temporal modes based on the projection of data onto the POD modes
    
    % Step 1: Calculate the covariance matrix
    C = X * X';  % Covariance matrix (n x n)

    % Step 2: Eigenvalue decomposition of the covariance matrix
    [PODMs, D] = eig(C);  % POD spatial modes and eigenvalues
    PODEvals = diag(D);  % Extract eigenvalues as a vector

    % Step 3: Sort the eigenvalues and corresponding modes in descending order
    [PODEvals, sortedIdx] = sort(PODEvals, 'descend');  % Sort eigenvalues in descending order
    PODMs = PODMs(:, sortedIdx);  % Sort POD modes accordingly

    % Step 4: Compute the norms of the POD modes (based on eigenvalues)
    SortPOD = sqrt(PODEvals);  % Norms based on sorted eigenvalues
    
    % Step 5: Compute the POD reconstructed temporal modes
    PODReTMs = PODMs' * X;  % (N x N) temporal modes
    
    % Normalize the temporal modes (optional step, depending on application)
    for i = 1:size(PODReTMs, 1)
        PODReTMs(i, :) = PODReTMs(i, :) / norm(PODReTMs(i, :));
    end
end
