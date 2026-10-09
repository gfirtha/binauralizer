classdef avendano_5ch_separator < avendano_separator
    %AVENDANO_5CH_SEPARATOR  Experimental 5-channel pair-bank analysis.
    %
    %   Uses adjacent 5.0/5.1 speaker pairs:
    %     L-C, C-R, R-Rs, Rs-Ls, Ls-L
    %   LFE is ignored. Each channel is split by sqrt power between its
    %   adjacent pairs. The output is five direct pair signals and five
    %   pair ambience/residual signals, intended for listening and
    %   diagnostics before more selective bin-level assignment.

    methods
        function obj = avendano_5ch_separator(varargin)
            obj@avendano_separator(varargin{:});
            obj.EnablePairBank = true;
            obj.PairBankMode = 'five_channel';
        end
    end
end
