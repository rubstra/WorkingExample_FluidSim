function [ gradient ] = calcGradient( pressure, obstacle )
%CALCGRADIENT Calculate the pressure gradient

% based on: https://github.com/candycat1992/2DFluidSim

% return newVelocityField: matrix(yNodes, xNodes, dim)
% pressure: matrix(yNodes, xNodes)
% obstacle: boolean matrix(yNodes, xNodes)

    mSize = size(pressure, 1);
    nSize = size(pressure, 2);
    
% make look-up at borders easier
% Padarray 'replicate', replicates the border of the matrix, increasing
% the matrix dimension by two in both dimensions.
    paddedPressure = padarray(pressure, [1 1]);
    paddedObstacle = padarray(obstacle, [1 1]);
    
% Find neighboring pressure

% Think of the matrice as constituting the entire grid, so that 4 adjacent
% elements define a cell. A method used to control the environment of cells
% or to calculate differences across cells, involve the generation of a
% matrice with shifted element positions. If the below matrix A defines 4 
% cells, we could calculate the differences at the top and bottom of the
% cell by defining a matrix B, extracted from A, and just compute A-B. As
% can be seen, B is just A excluding the top row.
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
    pressureLeft = paddedPressure(2 : mSize + 1, 3 : nSize + 2);
    pressureRight = paddedPressure(2 : mSize + 1, 1 : nSize);
    pressureUp = paddedPressure(3 : mSize + 2, 2 : nSize + 1);
    pressureDown = paddedPressure(1 : mSize, 2 : nSize + 1);
    
% Find neighboring obstacles
    obstacleLeft = logical(paddedObstacle(2 : mSize + 1, 3 : nSize + 2));
    obstacleRight = logical(paddedObstacle(2 : mSize + 1, 1 : nSize));
    obstacleUp = logical(paddedObstacle(3 : mSize + 2, 2 : nSize + 1));
    obstacleDown = logical(paddedObstacle(1 : mSize, 2 : nSize + 1));
    
% Use center pressure for solid cells
    pressureLeft(obstacleLeft) = pressure(obstacleLeft);
    pressureRight(obstacleRight) = pressure(obstacleRight);
    pressureUp(obstacleUp) = pressure(obstacleUp);
    pressureDown(obstacleDown) = pressure(obstacleDown);
    
% Note where forward / backward difference is used
% divide by 2 only for central differences
    differenceFactorHorizontal = 2 .* ones(size(pressure));
    differenceFactorHorizontal(obstacleLeft | obstacleRight) = 1;
    differenceFactorVertical = 2 .* ones(size(pressure));
    differenceFactorVertical(obstacleUp | obstacleDown) = 1;

% calc gradient with central differences
    gradient = zeros(mSize, nSize, 2);
    gradient(:, :, 1) = (pressureRight(:, :) - pressureLeft(:, :)) ./ ...
        differenceFactorHorizontal; 
    gradient(:, :, 2) = (pressureDown(:, :) - pressureUp(:, :)) ./ ...
        differenceFactorVertical;
    
% set gradient in obstacles to 0
% velocity in obstacles will be set to 0 anyway, helps with testing
    gradient(repmat(logical(obstacle), [1 1 2])) = 0; 
    
end
