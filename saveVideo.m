function saveVideo( visualSettings, deleteImages )
%SAVEVIDEO Create a video out of all .png images in the output folder.
%   visualisationSettings: struct with these fields:
%       outputFolder: folder where the resulting video should be saved
%       fps: frames per second of the video
%       outputFile: the name of the resulting video
%   deleteImages: bool delete or keep source images

    display('Saving video...');
    
    [~,~,~] = mkdir(visualSettings.outputFolder);
    
    if (exist([visualSettings.outputFolder '\' visualSettings.outputFile], 'file'))
        delete([visualSettings.outputFolder '\' visualSettings.outputFile]);
    end
    
    videoFile = [visualSettings.outputFolder '\' visualSettings.outputFile];
    
    % initialize video writer
    writerObj = VideoWriter(videoFile);
    writerObj.FrameRate = visualSettings.fps;
    open(writerObj);
    
    % open every single output frame and add it to the video
    file_list_frames = dir([visualSettings.outputFolder '/image*.png']);
    for j = 1:numel(file_list_frames)
        frame = imread([visualSettings.outputFolder '/' file_list_frames(j).name]);
        writeVideo(writerObj, frame);
    end
    
    if deleteImages
        delete([visualSettings.outputFolder '\image*.png']);
    end
    
    display('Finished!');
    
end

