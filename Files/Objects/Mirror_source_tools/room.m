classdef room < handle
    %ROOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        walls         % Normal vector of walls
        adjacency_matrix  % e.g to check closed volumes or corner effects
        reflection_transform
    end
    
    methods
        function obj = room(room_vertices, wall_vixs)
            %ROOM Construct an instance of this class
            %  room_vertices: vetices present in the room in matrix of column vectors
            %  wall_indices:  index of vertices for each wall
            %                    (2 row and N column for N 2D walls )
            for n = 1 : size( wall_vixs, 2 )
                obj.walls{n} = wall(n, room_vertices(:,wall_vixs(:,n)) );
            end
        end
        
        function mother_source = generate_mirror_source_tree(obj, mother_source, N)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cnt = 0;
            for n = 1 : N
                %%
                if n == 1
                    ixs = (1:length(obj.walls));
                    mirros_source_positions = cellfun(@(x) x.reflect(mother_source.position), obj.walls(ixs),'UniformOutput',false);
                    for j = 1 : length(ixs)
                        cnt = cnt + 1;
                        mother_source.image_sources{cnt} = image_source(j,mirros_source_positions{j},[0,0], n, 0, obj.walls{j}.index);
                    end
                else
                    ix_prev_layer = find(cell2mat(cellfun( @(x) find(x.depth==n-1),mother_source.image_sources,'UniformOutput',false)));
                    for k = 1 : length(ix_prev_layer)
                        curr_mother_source = mother_source.image_sources{ix_prev_layer(k)};
                        
                        % Reflect an image source only in planes that the image source
                        % is in front of
                        mw_ixs = curr_mother_source.mother_wall;
                        
                        % Image source is reflected to every wall except for its
                        % mother wall
                        o_w = cell2mat(cellfun(@(y) sum(y,2),(cellfun( @(x) x.vertices, obj.walls,'UniformOutput',false)),'UniformOutput',false))/2;
                        front_ixs = find(sum(bsxfun(@minus,o_w,curr_mother_source.position).*cell2mat(cellfun( @(x) x.normal, obj.walls,'UniformOutput',false)),1)>0);
                        
                        % Walls to reflect
                        ixs = setdiff((1:length(obj.walls)), union(front_ixs, mw_ixs));
                        
                        % Check reflected positions
                        mirros_source_positions = cellfun(@(x) x.reflect(curr_mother_source.position), obj.walls(ixs),'UniformOutput',false);
                        
                        for j = 1 : length(mirros_source_positions)
                            cnt = cnt + 1;
                            % TODO: save all reflecting wall history is
                            % better?
                            mother_source.image_sources{cnt} = image_source(length(mother_source.image_sources)+1,mirros_source_positions{j},[0,0], n, curr_mother_source.index, obj.walls{ixs(j)}.index);
                        end
                        
                    end
                end
            end
            
        end
    
        function reverberant_source = update_mirror_source_positions(obj, reverberant_source)
            for n = 1 : length(reverberant_source.image_sources)
                if reverberant_source.image_sources{n}.mother_source == 0
                     reverberant_source.image_sources{n}.position = ...
                         obj.walls{reverberant_source.image_sources{n}.mother_wall}.reflect(reverberant_source.position);              
                else
                     reverberant_source.image_sources{n}.position = ...
                         obj.walls{reverberant_source.image_sources{n}.mother_wall}.reflect(...
                         reverberant_source.image_sources{reverberant_source.image_sources{n}.mother_source}.position);              
                end
            end
        end
    end
end

