classdef avendano_ring_separator < avendano_separator
    %AVENDANO_RING_SEPARATOR  Experimental 5.1 pair-bank Avendano/Jot mode.
    %
    %   Uses four adjacent non-center speaker pairs:
    %     L/R, R/Rs, Rs/Ls, Ls/L
    %   The center channel is not used in this experimental output. LFE is
    %   ignored. Each non-center channel is split by sqrt power between its
    %   two adjacent pairs, producing four direct pair outputs and four
    %   pair ambience outputs.

    methods
        function obj = avendano_ring_separator(varargin)
            obj@avendano_separator(varargin{:});
            obj.EnablePairBank = true;
        end
    end
end
