function plotresultARTICLE(visualSettings, outputString,frames, divergenceORG, fluidpixelsORG, fluidpixelsISO, divergenceISO)

FramesToSec = linspace(0,10,1199);

figureFileName = [visualSettings.outputFolder '\' outputString];
figure('Name','Divergence Plot','NumberTitle','off')
plot(FramesToSec, divergenceORG)
title('Sum of divergence')
xlabel('Time [s]')
ylabel('Divergence')
hold on
plot(FramesToSec, divergenceISO)
legend('Quadrilateral','Triangle')
%savefig([figureFileName ' divergence.fig'])
print([figureFileName 'divergenceORGvsISO.png'],'-dpng')

figure('Name','Fluid Pixels Plot','NumberTitle','off')
plot(FramesToSec, fluidpixelsORG ./ 1000)
title('Mass loss')
xlabel('Time [s]')
ylabel('Fluid Pixels [*1000]')
hold on
plot(FramesToSec, fluidpixelsISO ./ 1000)
legend('Quadrilateral','Triangle')   
%savefig([visualSettings.outputFolder '\' outputString ' fluid pixels.fig'])
print([figureFileName ' fluid pixelsORGvsISO.png'],'-dpng')


%% save resulting video
%saveVideo(visualSettings, false);