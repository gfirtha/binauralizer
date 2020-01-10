classdef listener_space_axes < handle
    %LISTENER_SPACE_AXES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        axes
        receiver
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
            grid on
            xlabel('x -> [m]')
            ylabel('y -> [m]')
            xlim([-3,3]);
            ylim([-3,3]);
        end
        function receiver = draw_head(obj,pos,R)
            N = 20;
            fi = (0:2*pi/N:2*pi*(1-1/N));
            x_head = cos(fi)*R;
            y_head = sin(fi)*R;
            fi_nose =  [10 5   0    -5  -10 ];
            A_nose = R*[1  1.15 1.2 1.15   1 ];
            x_nose = cosd(linspace(fi_nose(1),fi_nose(end),N)).*interp1(fi_nose,A_nose,linspace(fi_nose(1),fi_nose(end),N),'spline');
            y_nose = sind(linspace(fi_nose(1),fi_nose(end),N)).*interp1(fi_nose,A_nose,linspace(fi_nose(1),fi_nose(end),N),'spline');
            fi_lear =  fliplr([18.5 18   16   0    -5   -15 ]+90);
            A_lear  = fliplr(R*[1  1.08 1.1 1.04  1.06   1 ]);
            x_lear = cosd(linspace(fi_lear(1),fi_lear(end),N)).*interp1(fi_lear,A_lear,linspace(fi_lear(1),fi_lear(end),N));
            y_lear = sind(linspace(fi_lear(1),fi_lear(end),N)).*interp1(fi_lear,A_lear,linspace(fi_lear(1),fi_lear(end),N));
            fi_rear =  -[18.5 18   16   0    -5   -15 ]-90;
            A_rear =  R*[1  1.08 1.1 1.04  1.06   1 ];
            x_rear = cosd(linspace(fi_rear(1),fi_rear(end),N)).*interp1(fi_rear,A_rear,linspace(fi_rear(1),fi_rear(end),N));
            y_rear = sind(linspace(fi_rear(1),fi_rear(end),N)).*interp1(fi_rear,A_rear,linspace(fi_rear(1),fi_rear(end),N));
            x_torso = cos(fi)*R*0.7-R/7;
            y_torso = sin(fi)*R*1.7;
            x_rec = [x_torso;x_head;x_lear;x_rear;x_nose]';
            y_rec = [y_torso;y_head;y_lear;y_rear;y_nose]';
            x_rec = x_rec - mean(mean(x_rec));
            y_rec = y_rec - mean(mean(y_rec));
            c = [37, 160, 217;
                77, 41, 14;
                255 206 180;
                255 206 180;
                255 206 180]/255;
            receiver = patch(gca, x_rec+pos(1) ,y_rec + pos(2),[0;1;1;1;1]);
            set(receiver,'FaceVertexCData',c);
            receiver.UserData = struct( 'Label', 1,...
                'Origin', [ mean(receiver.Vertices(:,1)), mean(receiver.Vertices(:,2))  ],...
                'Orientation', 0 );
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

