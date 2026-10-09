classdef brir_convolver < handle
    % BRIR_CONVOLVER  Zero-latency uniformly partitioned convolution of a mono
    % input with a two-ear filter, with a one-block crossfade on filter change.
    %
    % Both ears share one FFT/IFFT: the two real filters are packed as
    % H_L + 1i*H_R, so ifft(X.*(H_L + 1i*H_R)) = y_L + 1i*y_R.

    properties (SetAccess = protected)
        block_size
        N_fft            % 2*block_size
        H                % [N_fft x P] packed partition spectra (current)
        H_next = []      % pending filter, faded in on the next block
        FDL              % [N_fft x P] frequency-domain delay line
        prev             % previous input block
    end

    methods
        function obj = brir_convolver(block_size)
            obj.block_size = block_size;
            obj.N_fft = 2*block_size;
            obj.H = zeros(obj.N_fft, 1);
            obj.FDL = zeros(obj.N_fft, 1);
            obj.prev = zeros(block_size, 1);
        end

        function set_filter(obj, h)
            % h: [L x 2] two-ear impulse response
            Hp = obj.partition(h);
            if ~any(obj.H(:))
                obj.H = Hp;          % first filter: no fade needed
                obj.H_next = [];
            else
                obj.H_next = Hp;
            end
            P = max(size(obj.H,2), size(Hp,2));
            if size(obj.FDL,2) < P
                obj.FDL(:, end+1:P) = 0;
            end
        end

        function y = process(obj, x)
            % x: [block_size x 1] -> y: [block_size x 2]
            X = fft([obj.prev; x]);
            obj.prev = x;
            obj.FDL = [X, obj.FDL(:, 1:end-1)];
            y = obj.apply(obj.H);
            if ~isempty(obj.H_next)
                w = (0.5 - 0.5*cos(pi*(1:obj.block_size)'/obj.block_size));
                y = (1-w).*y + w.*obj.apply(obj.H_next);
                obj.H = obj.H_next;
                obj.H_next = [];
            end
        end
    end

    methods (Access = private)
        function y = apply(obj, Hp)
            P = size(Hp, 2);
            Y = sum(obj.FDL(:, 1:P) .* Hp, 2);
            yc = ifft(Y);
            yc = yc(obj.block_size+1:end);
            y = [real(yc), imag(yc)];
        end

        function Hp = partition(obj, h)
            B = obj.block_size;
            P = max(1, ceil(size(h,1)/B));
            h(end+1:P*B, :) = 0;
            hp = h(:,1) + 1i*h(:,2);
            Hp = fft(reshape(hp, B, P), obj.N_fft, 1);
        end
    end
end
