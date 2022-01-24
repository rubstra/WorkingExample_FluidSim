function plotresult(visualSettings, outputString, frames, divergence, fluidPixels)


figureFileName = [visualSettings.outputFolder '\' outputString];
figure('Name','Divergence Plot','NumberTitle','off')
plot(1:frames - 1, divergence)
title('Sum of divergence')
xlabel('Frame')
ylabel('Divergence')
savefig([figureFileName ' divergence.fig'])
print([figureFileName ' divergence.png'],'-dpng')

figure('Name','Fluid Pixels Plot','NumberTitle','off')
plot(1:frames - 1, fluidPixels .* 0.001)
title('Fluid Pixels *1000')
xlabel('Frame')
ylabel('Fluid Pixels')
savefig([visualSettings.outputFolder '\' outputString ' fluid pixels.fig'])
print([figureFileName ' fluid pixels.png'],'-dpng')
    
%% save resulting video
saveVideo(visualSettings, false);