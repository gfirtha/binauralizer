function [handles,hObject] = make_layout(handles, hObject, Nch)

mono = 0;
stereo = [30,-30];
triplet = [30,-30, 0];
dolby_51 = [30,-30,0,110,-110];
R = 1.5;

axes(handles.axes1);
switch Nch
    case 1
        [x,y] = pol2cart(mono*pi/180,R);
        handles.rois = drawpoint(handles.axes1,'Position',[x,y],'Color','blue');
        handles.sourcePosition = addlistener(handles.rois ,'MovingROI',@allevents);
        handles.output = hObject;
        guidata(hObject, handles);        
    case 2
        for n = 1 : Nch
            [x,y] = pol2cart(stereo(n)*pi/180,R);
            handles.rois{n} = drawpoint(handles.axes1,'Position',[x,y],'Color','blue');
            handles.sourcePosition(n) = addlistener(handles.rois{n} ,'MovingROI',@allevents);
            handles.output = hObject;
            guidata(hObject, handles);
        end
    case 3
        for n = 1 : Nch
            [x,y] = pol2cart(triplet(n)*pi/180,R);
            handles.rois{n} = drawpoint(handles.axes1,'Position',[x,y],'Color','blue');
            handles.sourcePosition(n) = addlistener(handles.rois{n} ,'MovingROI',@allevents);
            handles.output = hObject;
            guidata(hObject, handles);
        end
    case 5
        for n = 1 : Nch
            [x,y] = pol2cart(dolby_51(n)*pi/180,R);
            handles.rois{n} = drawpoint(handles.axes1,'Position',[x,y],'Color','blue');
            handles.sourcePosition(n) = addlistener(handles.rois{n} ,'MovingROI',@allevents);
            handles.output = hObject;
            guidata(hObject, handles);
        end
        
end

