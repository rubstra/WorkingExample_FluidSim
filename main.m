 
% Fluid Simulation with Finite Element Method


% Based on the work of Philip Erler and Christian Hafner, a further
% development of the implementation made by Erler.


%This is a working example developed for students studying fluid
%simulations and/or related topics. Many textbooks on fluid simulations are
%comprehensive for beginners, depending on your background. Understanding
%the theory could be challenging itself, and it often become more
%disturbing when this theory shall be put to practice. This working example
%aim at that gap, by presenting a code that solves a relatively simple
%fluid problem in both space and time. The code focus on describing what
%the different steps are performing, related to the theory, so that 
%in-depth knowledge of both the solution algorithm and mathematical
%programming can be gained.

%The example is related to the book by Robert Bridson, which is a very nice
%introduction to fluid simulations. Its called:
%               "Fluid Simulation For Computer Graphics"
%By reading the first small four chapters, and combining with this code,
%you should be able to get a quite good knowledge, step by step, both in
%theory and in practice.

%This example is developed for further developing, for this reason it is
%kept quite simple, so that it may introduce the basic and serve as a
%platform for new implementations and practice.


%% CLEAR AND CLOSE %%
clear global; % no "clear all" or breakpoints are gone
close all;
clc;


%% DISPLAY RUNTIME %%
% Can be set to true or false. If true you get the runtime of each routine
% within the simulation.
enableProfiling = true;

if enableProfiling
    profile on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% SIMULATION PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMAIN AND CELLS %

domainsizeXdir = 1; %m
domainsizeYdir = 1; %m
xCells = 99;
yCells = 99;
%only uniform cells
cellsize = domainsizeXdir/xCells;


% GENERAL %
simulationtime = 10; %s
fps = 120; %frames per second
frames = simulationtime*fps; %total number of frames
timestep = 1/fps; %(1 / frames per second) = seconds per frame = timestep
dissipation = 1;%can be ignored the first time, refer the book of Robert 
%Bridson when you want to learn more

% VISUALS &
%The visual settings are related to the postprocessor (what creates the
%visual from the simulation) it is not commented or elaborated in this
%working example.
[visualSettings,outputString] = visualsettings(fps,xCells,yCells);


%%%%%%%%%%%%%%%%%%%%%%
  %% SIMULATION %%
%%%%%%%%%%%%%%%%%%%%%%

disp('Based on the implementation from Philip Erler:');
disp(['Start simulation at: ', datestr(now)]);

% The below simulation function is the primary function for computing the 
% fluid problem. It is this separate .m file that you need to investigate 
% to understand the overall algorithm and solution procedure.
[divergence, fluidPixels] = simulation(...
    frames,...
    xCells,...
    yCells,...
    timestep,...
    cellsize,...
    dissipation,...
    visualSettings );
disp(['Finish simulation at: ', datestr(now)]);



%% Save and display runtime data %%
if enableProfiling
    profsave
end


%% PLOTTING %%
% This function plot the divergence and number of fluid pixels as a
% function of the frame number
plotresult(visualSettings, outputString, frames, divergence, fluidPixels);

