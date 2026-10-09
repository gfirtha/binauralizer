classdef OmnidirectionalDirectivity < Directivity
    % OMNIDIRECTIONALDIRECTIVITY  Uniform radiation in all directions.
    %   D(theta, f) = 1 for all angles and frequencies.

    methods
        function D = evaluate(~, theta, f)
            % theta : [M×1], f : [1×Nf]  →  D : [M × Nf]
            theta = theta(:);
            f = f(:)';
            D = ones(numel(theta), numel(f));
        end
    end
end
