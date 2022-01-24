function [velocityField, pressureField, obstacle, forceField ] = loadScenariotest
%LOADSCENARIO Create blobs and obstacles as defined in a scenario
%   return distanceField: matrix(yNodes, xNodes), distances to interface
%       velocityField: matrix(yNodes, xNodes, dim), velocity of
%           fluid in m/s
%       pressureField: matrix(yNodes, xNodes), pressure at each node
%       obstacle: bool matrix(yNodes, xNodes), node is in obstacle 
%           (no fluid or air)
%       forceField: matrix(yNodes, xNodes, dimension), force applied to
%           fluid in m/s^2
%   simulationSettings: struct with these fields:
%       xCells: grid cells in x direction
%       yCells: grid cells in y direction
%       cellSize: grid cell size in m
%       scenario: string, one of the construction functions below
%       deltaTime: passed time per step in seconds
    
    % never initial pressure
    xNodes = 100;
    yNodes = 100;
    pressureField = zeros(yNodes, xNodes);
    
    %which scenario
    [velocityField, obstacle, forceField] = singleBlocker;
  
    
    % cells to nodes
    %obstacle = logical(cellsToNodes(obstacle));
    %fluidNodes = cellsToNodes(fluidNodes);
    
    % calculate distanceField from waterbody
    %distanceField = fluidNodesToDistanceField(fluidNodes);

end




function [velocityField, obstacle, forceField] = ...
    singleBlocker

    xCells = 99%simulationSettings.xCells;
    yCells = 99%simulationSettings.yCells;
    xNodes = 100%simulationSettings.xCells + 1;
    yNodes = 100%simulationSettings.yCells + 1;
    
    % move everything down
%     velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 2) = 0.5;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
%     % apply gravity
%     forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % blob above box
    %waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.3, 0.25 );
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    % border around the fluid
%     obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    % box in middle
    obstacle(uint32(yCells * 0.65 : yCells * 0.85), ...
        uint32(xCells * 0.4 : xCells * 0.6)) = true;
    
end

