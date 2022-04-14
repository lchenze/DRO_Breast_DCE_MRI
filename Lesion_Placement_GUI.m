function location = Lesion_Placement_GUI(images,ss)

% display the images with scrolling and allow input
% prep
if nargin<2
    ss = round(size(images,3)/2);
end
maxS = size(images,3);
minS = 1;

if max(images(:)) > 1
    normimages = images;
    for count = 1:maxS
        im = images(:,:,count);
        imMax = max(im(:));
        normimages(:,:,count) = im./imMax; 
    end    
else
    normimages = images;
end
% LPfig = findobj('Type','Figure','Name','Lesion Placement');

% % create gui
% if ~isempty(LPfig)
%     figure(LPfig)
%     slider1 = findobj('Style','slider');
% else
LPfig = figure('Menubar','none',...
                'Name','Lesion Placement','NumberTitle','off',...
                'IntegerHandle','off','units','normalized',...
                'outerposition',[0.15 0 .7 1],'Visible','off');

% LPfig = figure('Menubar','none',...
%                 'Name','Lesion Placement','NumberTitle','off',...
%                 'IntegerHandle','off','units','normalized',...
%                 'outerposition',[0.15 0 .7 1]);

% image selection slider
slider1 = uicontrol(LPfig,'Style','slider','Callback',@slider1_callback);
set(slider1,'value',ss);
set(slider1,'max',maxS);
set(slider1,'min',minS);
set(slider1,'units','normalized','Position',[.4 .05 .2 .03]);
set(slider1, 'SliderStep', [1/maxS, 5/maxS]);
% end

% initialize display
slice = get(slider1,'Value');
imshow(normimages(:,:,slice));
titl = 'Click to place a lesion';
title([titl ' \newline ' 'Slice ' num2str(slice)],'FontSize',16)
drawnow
set(LPfig,'Visible','on');

% get corrections
clicked = false;
while ~clicked
    pH = impoint;
    try
        pos = pH.getPosition;
        if max(pos) < max(size(images(:,:,1)))
            clicked = true;
            location = round([pos get(slider1,'Value')]);
        end
    catch
    end
end

close(LPfig)

%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
    function slider1_callback(slider1,~,~)
        disp(get(slider1,'Value'))
        slice = round(get(slider1,'Value'));
        imshow(normimages(:,:,slice))
        title(['Click to place a lesion'...
            ' \newline ' 'Slice ' num2str(slice)],'FontSize',16)
        drawnow
    end
end