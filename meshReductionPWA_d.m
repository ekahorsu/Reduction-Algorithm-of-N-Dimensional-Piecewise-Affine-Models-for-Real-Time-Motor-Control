function [union_2d, Rempoints] = meshReductionPWA_d(input,output,N)
    % Initial setup
    pwa2 = getMPT3Union(input', output');
    points = input;
    numPoints = size(points, 1);
    remainingPoints = points; % Initialize remaining points with all points
    n = size(input,2)
    
    % Treat boundary points as fixed points
    if n == 1
        fixedPoints = [min(input); max(input)];
    else
        boundary_indices = convhull(input, 'Simplify', true);
        boundary_indices(end) = [];
        fixedPoints = input(boundary_indices, :)
    end

    % Generate list of points to remove
    l2Norm = zeros(size(remainingPoints, 1), 1);
    for i = 1:size(remainingPoints, 1)
        % Check if the point is a boundary point
        isBoundaryPoint = any(ismember(fixedPoints(:, :), remainingPoints(i, :), 'rows'));
        if isBoundaryPoint
            l2Norm(i) = inf; % Set L2 norm to infinity for fixed points
        else    
            % Get the indices of polyhedra common to that the point       
            polyhedraIndices = [];
            for j = 1:numel(pwa2)
                if any(ismember(pwa2(j).V(:, 1:n), points(i, :), 'rows'))
                    polyhedraIndices = [polyhedraIndices, j];
                end
            end

            % Compute the sum of pairwise gradient vector Euclidean
            % distances for the point
            numPairs = numel(polyhedraIndices);
            pairwiseDistancesSum = 0;
            for j = 1:numPairs
                for k = (j+1):numPairs
                    pairwiseDistancesSum = pairwiseDistancesSum + norm(pwa2(polyhedraIndices(j)).getFunction('f').F - pwa2(polyhedraIndices(k)).getFunction('f').F);
                end
            end

            % Calculate the L2 norm (similarity value) as the average of the sum of pairwise distances
            numCombinations = (numPairs * (numPairs - 1)) / 2; % Calculate the number of combinations
            if numCombinations > 0
                l2Norm(i) = pairwiseDistancesSum / numCombinations;
            else
                % Handle the case with only one vector (Euclidean distance between vector and zero vector)
                l2Norm(i) = norm(pwa2(polyhedraIndices(1)).getFunction('f').F);
            end
        end
    end
    
    % Sort the points based on the lowest representative L2 norm values
    [sortedValues, sortedIndices] = sort(l2Norm, 'ascend');
    bestPoints = remainingPoints(sortedIndices, :);
    
    % Select the points to delete
    numPointsToDelete = size(remainingPoints, 1) - N;
    remainingPoints = [bestPoints(numPointsToDelete+1:end,:)];

    % Generate the remaining input and output matrices
    remainingInput = input(ismember(input(:, 1:n), remainingPoints, 'rows'), :);
    remainingOutput = output(ismember(input(:, 1:n), remainingPoints, 'rows'), :);

    % Create the simplified PWA objects with the remaining points
    simpPWA = getMPT3Union(remainingInput', remainingOutput');

    % Combining PWA objects into PWA map
    union_2d = PolyUnion(simpPWA);
    Rempoints = size(remainingPoints, 1);
end
