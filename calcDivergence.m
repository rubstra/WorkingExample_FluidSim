function [ divergence ] = calcDivergence( velocity, obstacle )
%CALCDIVERGENCE Calculate the divergence from a velocity field.
%   return divergence: matrix(yNodes, xNodes)
%   velocity: matrix(yNodes, xNodes, dim)
%   obstacle: boolean matrix(yNodes, xNodes)

    ySize = size(velocity, 1);% number of rows
    xSize = size(velocity, 2);% number of columns
    
    % make look-up at borders easier
    % Padarray 'replicate', replicates the border of the matrix, increasing
    % the matrix dimension by two in both dimensions.
    paddedVelocity = padarray(velocity, [1 1 0], 'replicate');
    paddedObstacle = padarray(obstacle, [1 1], 'replicate');
    
    % Find neighboring velocity
%Think of the matrice as constituting the entire grid, so that 4 adjacent
%elements define a cell. A method used to control the environment of cells
%or to calculate differences across cells, involve the generation of a
%matrice with shifted element positions. If the below matrix A defines 4 
%cells, we could calculate the differences at the top and bottom of the
%cell by defining a matrix B, extracted from A, and just compute A-B. As
%can be seen, B is just A excluding the top row.
%
% A  1 2 3 4  B 5 4 3 2
%    5 4 3 2    6 7 8 9
%    6 7 8 9    1 9 8 7
%    1 9 8 7 
%
%Obviously A and B have different dimensions, so to be able to combine them
%we have already padded the matrices, to get larger matrices. The only 
%thing to do then is as below, extract the different parts so that we can 
%combine them to calculate differences between top and bottom and left- and
%right side. The left,right,up and down refer to the parts of the cell 
%being excluded from the matrix. Below we do this operation for the
%velocity fields, as well as the obstacle matrix. But we only define the
%matrices, it is not until the bottom differencing that we use the matrices
%to compute the diffenreces as explained above.

% Defning the matrices to be used for calculating differences across cells
    velocityLeft = paddedVelocity(2 : ySize + 1, 3 : xSize + 2, :);
    velocityRight = paddedVelocity(2 : ySize + 1, 1 : xSize, :);
    velocityUp = paddedVelocity(3 : ySize + 2, 2 : xSize + 1, :);
    velocityDown = paddedVelocity(1 : ySize, 2 : xSize + 1, :);
    
% Find neighboring obstacles
    obstacleLeft = logical(paddedObstacle(2 : ySize + 1, 3 : xSize + 2));
    obstacleRight = logical(paddedObstacle(2 : ySize + 1, 1 : xSize));
    obstacleUp = logical(paddedObstacle(3 : ySize + 2, 2 : xSize + 1));
    obstacleDown = logical(paddedObstacle(1 : ySize, 2 : xSize + 1));
    
% Set velocities to 0 for solid cells
    velocityUp(obstacleUp) = 0;
    velocityDown(obstacleDown) = 0;
    velocityRight(obstacleRight) = 0;
    velocityLeft(obstacleLeft) = 0;
    
% Calcualte the differences of x + difference of y around the cell
% Horizontal velocities are stored within velocity(:,:,1) while vertical
% y-velocites within velocity(:,:,2)
    divergence = 0.5 * (velocityRight(:,:,1) - velocityLeft(:,:,1) + ...
        velocityUp(:,:,2) - velocityDown(:,:,2));
    
end

