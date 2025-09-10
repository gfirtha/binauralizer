classdef brir_extrapolator < handle
    %HRTF_EXTRAPOLATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        measurement_grid
        HRTF_database
    end

    methods
        function obj = brir_extrapolator(hrtf)
            %HRTF_EXTRAPOLATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.HRTF_database = hrtf;

            name = strcat(obj.HRTF_database.GLOBAL_DatabaseName,'_',obj.HRTF_database.GLOBAL_ListenerShortName);
            if all(abs(hrtf.SourcePosition(:,2))<1e-9)
                obj.measurement_grid = circular_design(obj.HRTF_database.EmitterPosition,name);
            else
                obj.measurement_grid = spherical_design('Direct',obj.HRTF_database.EmitterPosition,name);
            end

        end

        function hrtf_out = extrapolate_hrtf(obj,source, receiver)
            xr = [receiver.position,0];
            xs = [source.position, source.height];
            o = [receiver.orientation,0];
            [~,ix_o] = min(sum(abs(obj.HRTF_database.ListenerView-o),2));
            xe = (get_rotation_mx(o',[0,1,0]')*(xs-xr)')';
            BRIR0 =  squeeze(obj.HRTF_database.Data.IR(ix_o,:,:,:));
            hrtf_out = obj.measurement_grid.interpolate(xe',permute( BRIR0,[2,1,3]),'nearest')*...
                mean(obj.measurement_grid.Phi0(:,end))/norm(xe);
        end
    end
end

