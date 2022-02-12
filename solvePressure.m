function [ pressure ] = solvePressure( xCells,yCells, velocity, distanceField, obstacle )
%SOLVEPRESSURE Apply pressure solver with finite elements method on 
%   first order quadrilateral elements, squares or rectangles depending on 
%   cell sizes.
%   
%   
%   For more information on the solver and finite elements methods a
%   recommended book isThe Finite Element Method: Theory, Implementation,
%   and Applications by Larson and Bengzon. The book contain several minor
%   errors in the examples, but explain the theory in a very understandable
%   manner.
%
%   OUTPUT: pressure matrix (yNodes, xNodes) containing the pressure at 
%           each node, solved by the finite element method
%
%    INPUT: DistanceField: A matrix with indices representing nodes(yNodes, 
%           xNodes). The nodal values denote the distances to interfaces, see
%           further descriptions in loadscenario or other.

%           Velocity: velocity field matrix(yNodes, xNodes,[dim]), velocity 
%           of fluid in m/s at the nodes.
%
%           Obstacle: Boolean matrix, obstacles are represented with
%           1's, these are domain borders and obstructions


%%%%%%%%%%%%%%%%%%%%   
%% initialization %%
%%%%%%%%%%%%%%%%%%%%


% Note that we use y,x only when accessing matrizes from caller, internal 
% is x,y
    numCellsX = xCells;
    numCellsY = yCells;
    
    velocityFieldX = velocity(:,:,1);
    velocityFieldY = velocity(:,:,2);
    
    % gauss quadrature sample locations
    s = zeros(4, 2);
    s(1,:) = [-sqrt(1/3); -sqrt(1/3)];
    s(2,:) = [-sqrt(1/3); +sqrt(1/3)];
    s(3,:) = [+sqrt(1/3); -sqrt(1/3)];
    s(4,:) = [+sqrt(1/3); +sqrt(1/3)];
    zeta = s(:, 1);
    eta = s(:, 2);

    % Allocating matrices
    % sparse matrix constructor parameters
    Kx = [];
    Ky = [];
    Kval = [];
    
    % pressure vector as sparse matrix
    fx = [];
    fval = [];
    
    % Get fluid cells and interface nodes
    % Convert a distance field, based on nodes, to a boolean matrix of 
    % cells. See getFluidCells for further description
    fluidCells = getFluidCells(distanceField);
    % See getInterfaceNodes for further description
    interfaceNodes = getInterfaceNodes(distanceField, obstacle, true);
    
    %Defining and allocating 
    nodeIndex = 0;
    nodeCoordToIndex = zeros(size(distanceField, 1), size(distanceField, 2));
    nodeIndexToCoord = zeros(size(distanceField, 1) * size(distanceField, 2), 2);
    
    % if we have a cell corresponding to 2,2 (row, column) then the node
    % offset would give the nodes surrounding that cell. Take the first row
    % vector 0 0 and add by elements to 2 2 it equals 2 2. Then node 2 2 is
    % the upper left node in the cell. Likewise 1,1 result in node 3,3
    % which is the bottom right node of the cell. The offsets refer to a
    % real coordinate system, not the matrix elements
    nodeOffsets = [0 0; 0 1; 1 1; 1 0];
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %% pressure solving %%
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    % for each cell in y-dir
    for cellY = 1 : numCellsY 
        % for each cell in x-dir
        for cellX = 1 : numCellsX
            % surrounded by 4 nodes
            % skip cells without water
            if ~fluidCells(cellY, cellX)
                continue;
            end
            
            %coordinates for the current iteration
            currentCoord = [cellX cellY];
            
            %Allocating/defining vectors
            wTildeEX = [0; 0; 0; 0];
            wTildeEY = [0; 0; 0; 0];
            
            % gather all information for cell
            for corner = 1 : 4 % run onces per corner/node of the cell
                % gives corner coordinates, see nodeOffsets described above
                cornerCoord = currentCoord + nodeOffsets(corner, :);
                
                % Temporarily store nodal x and y components of velocity in
                % respective vectors
                wTildeEX(corner) = velocityFieldX(cornerCoord(2), cornerCoord(1));
                wTildeEY(corner) = velocityFieldY(cornerCoord(2), cornerCoord(1));
                
                % add this node to indexing
                if nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) == 0
                    nodeIndex = nodeIndex + 1;
                    
                    nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) = nodeIndex;
                    nodeIndexToCoord(nodeIndex, :) = cornerCoord;
                end
            end
            
            % integrate with gauss quadrature
            Ke = 0; % 4x4 matrix
            fe = 0; % 4 dim vector
            wi = 1; % weights for gauss quadrature, here always 1
            for corner = 1 : 4
                currG = G(zeta(corner), eta(corner));
                Bsi = B(cellX, cellY, currG);
                JsiDet = det(J(cellX, cellY, currG));
                shapeFunctions = N(zeta(corner), eta(corner));
                BsiX = Bsi(:,1)';
                BsiY = Bsi(:,2)';
                
                Ke = Ke + wi * JsiDet * (Bsi * Bsi');
                fe = fe + wi * JsiDet * ...
                     (BsiX * wTildeEX + BsiY * wTildeEY) * shapeFunctions;
            end
            
            % add this element to the global matrix and vector
            newKx = zeros(16, 1);
            newKy = zeros(16, 1);
            newKval = zeros(16, 1);
            newFx = zeros(4, 1);
            
            vectorIndex = 1;
            for fromCorner = 1 : 4
                fromCornerCoord = currentCoord + nodeOffsets(fromCorner, :);
                fromCornerNodeIndex = nodeCoordToIndex(fromCornerCoord(2), fromCornerCoord(1));
                
                for toCorner = 1 : 4
                    toCornerCoord = currentCoord + nodeOffsets(toCorner, :);
                    toCornerNodeIndex = nodeCoordToIndex(toCornerCoord(2), toCornerCoord(1));
                
                    newKx(vectorIndex) = fromCornerNodeIndex;
                    newKy(vectorIndex) = toCornerNodeIndex;
                    newKval(vectorIndex) = Ke(fromCorner, toCorner);
                    vectorIndex = vectorIndex + 1;
                end
                
                newFx(fromCorner) = fromCornerNodeIndex;
            end
            Kx = [Kx; newKx];%adding newKx to the groing columnvector Kx
            Ky = [Ky; newKy];
            Kval = [Kval; newKval];
            
            fx = [fx; newFx];
            fval = [fval; fe];
            
        end % y cells
    end % x cells
    
    [Kx, Ky, Kval, fx, fval] = removeBoundaryNodes(...
          Kx, Ky, Kval, fx, fval, nodeIndexToCoord, interfaceNodes);
      
    % map node indieces to matrix collumn indieces
    activeNodeIndieces = unique(fx);
    numActiveNodes = size(activeNodeIndieces, 1);
    nodeIndexToCollumn(activeNodeIndieces) = 1:numActiveNodes;
    Kx = nodeIndexToCollumn(Kx);
    Ky = nodeIndexToCollumn(Ky);
    fx = nodeIndexToCollumn(fx);
    
    % get pressurefield
    K = sparse(Kx, Ky, Kval);
    f = sparse(fx, ones(size(fx)), fval);
    
    pressureAtNodes = K \ f;    % solve
    
    % should never be necessary in a working system
    %pressureAtNodes(isnan(pressureAtNodes)) = 0; 
    %pressureAtNodes(isinf(pressureAtNodes)) = 0;
    
    % save resulting pressure values from node vector to matrix
    pressure = zeros(size(distanceField));
    currNodeIndex = 1 : numActiveNodes;
    currNodeY = nodeIndexToCoord(activeNodeIndieces(currNodeIndex), 2);
    currNodeX = nodeIndexToCoord(activeNodeIndieces(currNodeIndex), 1);
    matrixIndices = sub2ind(size(pressure), currNodeY, currNodeX);
    pressure(matrixIndices) = pressureAtNodes(currNodeIndex);
    
    %maxabspressure = max(max(abs(pressure)))
    %minpressure = min(min(pressure))
end

% shape functions of a 4-node quad element, parameters are natural
% (local) coordinates of the element
% corners are at (-1,-1), (-1,1), (1,1), (1,-1)
function res = N1(zeta, eta)
    res = 0.25 * (1-zeta) * (1-eta);
end
function res = N2(zeta, eta)
    res = 0.25 * (1-zeta) * (1+eta);
end
function res = N3(zeta, eta)
    res = 0.25 * (1+zeta) * (1+eta);
end
function res = N4(zeta, eta)
    res = 0.25 * (1+zeta) * (1-eta);
end

% vector of shape functions
function res = N(zeta, eta)
    res = [N1(zeta, eta); N2(zeta, eta); N3(zeta, eta); N4(zeta, eta)];
end

% % linear interpolation in quad element, bottom
% function res = qBottom(eta, q1, q4)
%     res = 0.5 * ((1 - eta) * q1 + (1 + eta) * q4);
% end
% 
% % linear interpolation in quad element, top
% function res = qTop(eta, q2, q3)
%     res = 0.5 * ((1 - eta) * q2 + (1 + eta) * q3);
% end
% 
% % bi-linear interpolation in quad element
% function res = q(eta, zeta, q1, q2, q3, q4)
%     res = 0.5 * ((1 - zeta) * qBottom(eta, q1, q4) + ...
%                  (1 + zeta) * qTop(eta, q2, q3));
% end

% matrix of partially derived shape functions
function res = G(zeta, eta)
    res = 0.25 .* [(+eta-1), (+zeta-1); 
                   (-eta-1), (-zeta+1); 
                   (+eta+1), (+zeta+1); 
                   (-eta+1), (-zeta-1)];
end

% global coordinates of nodes
function res = Xtilde(cellX, cellY)
    res = [cellX, cellX,   cellX+1, cellX+1;
           cellY, cellY+1, cellY+1, cellY   ];
end

% jakobi matrix of cell
function res = J(cellX, cellY, currG)
    res = Xtilde(cellX, cellY) * currG;
end

% inverse jakobi matrix of cell
function res = invJ(cellX, cellY, currG)
    res = inv(J(cellX, cellY, currG));
end

% partially derived shape functions with respect to global coordinates
function res = B(cellX, cellY, currG)
    res = currG * invJ(cellX, cellY, currG);
end

function [newKx, newKy, newKval, newfx, newfval] = removeBoundaryNodes(...
          Kx, Ky, Kval, fx, fval, nodeIndexToCoord, interfaceNodes)

    % remove boundary nodes from matrix and vector, pressure must be 0
    % keep them if they are in obstacles
    fromNodeCoords = nodeIndexToCoord(Kx, :);
    toNodeCoords = nodeIndexToCoord(Ky, :);

    fromNodeInBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        fromNodeCoords(:,2), fromNodeCoords(:,1)));
    toNodeInBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        toNodeCoords(:,2), toNodeCoords(:,1)));
    
    sparseIndiecesToRemove = fromNodeInBoundary | toNodeInBoundary;
        
    Kval(sparseIndiecesToRemove) = [];
    Kx(sparseIndiecesToRemove) = [];
    Ky(sparseIndiecesToRemove) = [];
    
    nodeCoords = nodeIndexToCoord(fx, :);
    inBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        nodeCoords(:,2), nodeCoords(:,1)));
    sparseIndiecesToRemove = inBoundary;
    
    fx(sparseIndiecesToRemove) = [];
    fval(sparseIndiecesToRemove) = [];
    
    newKx = Kx;
    newKy = Ky;
    newKval = Kval;
    newfx = fx;
    newfval = fval;
    
end
