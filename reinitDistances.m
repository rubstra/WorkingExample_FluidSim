function [ newDistanceField ] = reinitDistances( distanceField, obstacle )
%REINITDISTANCES Re-initialize the distance field. Nodes at the interface 
%       keep their distance. Without this function, the interpolation in 
%       the advection would cause the distances to blend.
%   return newDistanceField: matrix(yNodes, xNodes)
%   distanceField: matrix(yNodes, xNodes)
    
    % get fluid cells and interface nodes
    interfaceNodes = getInterfaceNodes(distanceField, obstacle, false);
    
    % get water filled areas
    fluidNodes = (distanceField>=0);
    
    % assemble new distance field
    % these are the very same lines as describen in loadscenario.m, see 
    % mentioned file for further description.      
        newDistanceField = zeros(size(fluidNodes));%create matrix
        distancesToWater = (bwdist(fluidNodes) - 0.5) .* -1;
        distancesToAir = (bwdist(~fluidNodes) - 0.5);
        newDistanceField(~fluidNodes) = distancesToWater(~fluidNodes);
        newDistanceField(fluidNodes) = distancesToAir(fluidNodes);
    
    
    % remove infinite values (created by bwdist if no water exists)
    newDistanceField(isinf(newDistanceField)) = -1;
    
    % keep distance at the interface
    newDistanceField(interfaceNodes) = distanceField(interfaceNodes);

end

