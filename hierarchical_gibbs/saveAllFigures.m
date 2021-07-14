function [] = saveAllFigures(path)
%SAVEALLFIGURES Summary of this function goes here
%   Detailed explanation goes here
    if not(isfolder(path))
        mkdir(path)
    end

    figList = findobj(allchild(0), 'flat', 'Type', 'figure');

    for i = 1:length(figList)
        handle = figList(i);
        
        figName = get(handle, 'Name');
        
        if isempty(figName)
            figName = ['figure' num2str(get(handle, 'Number'))];
        end
        
        set(handle, 'Units', 'Inches');
        pos = get(handle, 'Position');
        set(handle, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches',...
            'PaperSize', [pos(3), pos(4)])
        
        saveFile = fullfile(path, figName);
        disp(saveFile)
        print(handle, saveFile, '-dpdf', '-r0');
        savefig(handle, saveFile);
    end
end

