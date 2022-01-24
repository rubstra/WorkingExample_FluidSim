function [ fluidCells ] = getFluidCells( distanceField )
%
%
%GETFLUIDCELLS Convert a distance field, based on nodes, to a boolean 
%       matrix of cells. The conversion from node to cells involve matrix
%       reduction. It does so by adding 4 boolean matrices, each 
%       representing the internal matrix of dimension yCells*xCells, taken 
%       from the four individual corners of the larger node-based
%       distancefield. This result in all nodes being included, so that all
%       nodes affect their neighbouring cells. 
%       For a cell to contain fluid it must contain at least one fluid node.
%       No difference is made wih regards to the number of nodes with water
%       that result in water within a cell. That is, mixed cells, air and
%       water, are treated as completely water filled.
%
%OUTPUT:
%       fluidCells: boolean matrix (yCells, xCells) indicating water by
%       use of 1's (true)
%INPUT:
%       distanceField: matrix(yNodes, xNodes)

    xCells = size(distanceField, 1) - 1;
    yCells = size(distanceField, 2) - 1;
    
    % surrounded by 4 nodes
    topLeftFluidNodes = ...
        (distanceField(1 : xCells, 1 : yCells)>=0);
    topRightFluidNodes = ...
        (distanceField(1 : xCells, 2 : yCells+1)>=0);
    bottomLeftFluidNodes = ...
        (distanceField(2 : xCells+1, 1 : yCells)>=0);
    bottomRightFluidNodes = ...
        (distanceField(2 : xCells+1, 2 : yCells+1)>=0);
    
    % mark water filled cells
    % treat mixed cells like completely water filled
    fluidCells = topLeftFluidNodes | bottomLeftFluidNodes | ...
                   topRightFluidNodes | bottomRightFluidNodes;
            
end

