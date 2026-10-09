classdef wall < handle
    %WALL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        vertices
        normal
        impedance
        reflection_transform
    end
    
    methods
        function obj = wall(ix, vertices )
            %WALL Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = ix;
            obj.vertices = vertices;
            
            v_w =  (vertices(:,1)-vertices(:,2))/norm(vertices(:,1)-vertices(:,2));
            obj.normal = [-v_w(2,:); v_w(1,: )]; % normal vector
            
            D = norm( vertices(:,1) - (vertices(:,1)'*v_w)*v_w  );
            
            Tr_mx  = eye(2) - 2 * obj.normal*obj.normal';
            v_trans = (Tr_mx-eye(2))*(D*obj.normal);
            
            obj.reflection_transform = struct('Transform_matrix',Tr_mx,'Translation_vector',v_trans);
            obj.impedance = 1;
        end
        
        function xsr = reflect(obj, xs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xsr = obj.reflection_transform.Transform_matrix*xs +...
                  obj.reflection_transform.Translation_vector;
        end
        function xsr = update(obj, dxs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xsr = xs + obj.reflection_transform.Transform_matrix*dxs;
        end
    end
end

