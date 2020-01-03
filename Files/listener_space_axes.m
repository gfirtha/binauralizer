classdef listener_space_axes < handle
    %LISTENER_SPACE_AXES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        axes
        virtual_source_points
        binaural_source_points
    end
    
    methods
        function obj = listener_space_axes()
            %LISTENER_SPACE_AXES Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function draw_gui(obj,gui_handle)
            obj.axes = gui_handle;
            axes(obj.axes);
            plot(0,0,'ok','MarkerSize',25,'MarkerFaceColor',[255,206,180]/255)
            hold on
            quiver(0,0,0.2,0, 0, 'AutoScale','off','MaxHeadSize',2,'LineWidth',2,'LineStyle','-','Color',[255,206,180]/255)
            grid on
            xlabel('x -> [m]')
            ylabel('y -> [m]')
            xlim([-3,3]);
            ylim([-3,3]);
        end
        
        function [rois,sourcePosition] = update_sources(obj,setup)
            for n = 1 : size(setup.Default_source_position,1)
                rois{n} = drawpoint(obj.gui_handle,...
                    'Position',setup.Default_source_position(n,:),'Color','blue');
                sourcePosition(n) = addlistener(rois{n} ,'MovingROI',@allevents);
            end
        end
        
        function [rois,sourcePosition] = add_sources(obj,setup,gui_handle)
            for n = 1 : size(setup.Default_source_position,1)
                rois{n} = drawpoint(gui_handle,...
                    'Position',setup.Default_source_position(n,:),'Color','blue');
                sourcePosition(n) = addlistener(rois{n} ,'MovingROI',@allevents);
            end
        end
    end
end

