classdef Physics < handle
    % Global physics constants (singleton, observable)
    properties (SetObservable)
        c (1,1) double {mustBeFinite,mustBePositive} = 343.1;  % m/s
    end

    methods (Access = private)
        function obj = Physics, end
    end

    methods (Static)
        function h = get()
            % Access the unique instance anywhere: audioapp.util.Physics.get()
            persistent s
            if isempty(s) || ~isvalid(s), s = audioapp.util.Physics; end
            h = s;
        end

        function val = speedOfSound()
            val = audioapp.util.Physics.get().c;
        end

        function setSpeedOfSound(val)
            audioapp.util.Physics.get().c = val;
        end

        function setFromTemperature(T_celsius)
            % handy approximate (dry air)
            val = 331.3 + 0.606*T_celsius;       % m/s
            audioapp.util.Physics.setSpeedOfSound(val);
        end
    end
end
