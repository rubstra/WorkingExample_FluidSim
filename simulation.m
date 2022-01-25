function [divergence, fluidPixels] = simulation( ...
    frames,...
    xCells,...
    yCells,...
    timestep,...
    cellsize,...
    dissipation,...
    visualSettings )
%SIMULATION Run the simulation and solve the fluid problem. The overall
%solution algorithm can be withdrawn from the primary steps of this file

%the simulation saves images every timestep, combined in the end to make a
%video presentation of the flow problem
%
%OUTPUT :
%   divergence: vector(frames - 1) contains sum of divergence computed for
%   within the domain every timestep, that is, every frame.
%
%   fluidPixels: vector(frames - 1) contains number of fluid pixels * 1000
%   it is related to the divergence in the sense of loosing matter (fluid)
%   from the simulation. If we have divergence the number og water fluid
%   pixels will change (likely decrease).
%
%INPUT:
%       frames: total number of frames, or, total number of timesteps to
%       reach overall simulation time.
%
%       cellSize: grid cell size in m
%       timestep: passed time per step in seconds, same as deltatime
%
%       advectionInterpolation: string, can be chosen to be 
%       'spline', 'cubic', 'linear' or 'nearest'
%
%       dissipation: friction within fluid, set to 1.
%       
%       visualSettings: struct with multiple fields, see visualsettings.m
%
%
%


%loadscenario contain the description of the physical domain, as well as
%initial conditions
    display('Loading scenario...');
    
    [ distanceField,... 
      velocityField,...
      pressureField,...
      obstacle,...
      forceField] = loadScenario(xCells,yCells);
    
    display('Simulating...');
    
    
    
    % save start condition
    frame = calcFrame(visualSettings,...
                      distanceField, ...
                      obstacle,...
                      velocityField,...
                      pressureField);
                  
    saveImage(frame, 0, visualSettings.outputFolder);
        
     % Create matrices
     divergence = zeros(1, frames - 1);
     fluidPixels = zeros(1, frames - 1);
        
    % delete velocities in obstacles (set obstacle velocity = 0). 
    % The built-in Matlab function repmat is utilized here. obstacle must 
    % be a logical matrix (boolean). Below line compare positions in both
    % matrices and for any position where obstacle correspond to value 1,
    % the value 0 is placed in velocityField.
    velocityField(repmat(obstacle, [1 1 2])) = 0;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %% SIMULATION ALGORITHM %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %The following loop is repeated once per timestep, one loop will result
    %in a new frame.
    
    for i = 1 : frames - 1
        %% simulation
        
        % STEP 1 - apply forces
        velocityField = velocityField + forceField * timestep;
        
        % STEP 2 - again, delete velocities in obstacles, as above
        velocityField(repmat(obstacle, [1 1 2])) = 0;
                
        % STEP 3 - calculate pressurefield, pressure correction, by use of 
        % Finite element method
        pressureField = solvePressureISO(xCells, yCells, velocityField, distanceField, obstacle);
        
        % STEP 4 - apply pressure correction (projection) to keep the flow 
        % divergence free, that is, keep as much water volume as possible
        gradient = calcGradient(pressureField, obstacle );
        velocityField = velocityField - gradient; 
                
        % STEP 5 - correct air-velocities, necessary for advection
        velocityField = setAirVelocity(velocityField, distanceField, obstacle);
        
        % STEP 6 - Advect density and velocities
        %     Advect distancefield
        distanceField = advectSemiLagrange(cellsize, timestep, dissipation, velocityField, distanceField);
        %     Advect velocity
        velocityField = advectSemiLagrange(cellsize, timestep, dissipation, velocityField, velocityField);
        
        %% visualization
        % re-init distance field after n frames
        if visualSettings.distanceRefreshInterval > 0 && ...
                mod(i, visualSettings.distanceRefreshInterval) == 0
            distanceField = reinitDistances(distanceField, obstacle);
        end

        
        [frame, fluidPixels(i)] = calcFrame(visualSettings, distanceField, ...
            obstacle, velocityField, pressureField);
        saveImage(frame, i, visualSettings.outputFolder);
        
        % store overall divergence (should be 0 after perfect pressure solving)
        waterBody = (distanceField>=0);
        frameDivergence=calcDivergence( velocityField, obstacle );
        frameDivergence(obstacle)=0; % care only for fluid
        frameDivergence(~waterBody)=0;
        divergence(i) = sum(sum(frameDivergence));
    
        disp(strcat('Frame:', num2str(i),'/',num2str(frames-1)));
    end
    
    display('Simulation finished!');
end

