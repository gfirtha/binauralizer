classdef time_delay_renderer < base_renderer
    %TIME_DELAY_RENDERER
    % Simple "time-delay stereo" renderer:
    % - Chooses the two most "frontal" SSD loudspeakers as a stereo base.
    % - Applies opposite phase ramps → time delays that steer towards the
    %   virtual source as seen from a reference point (sweet spot).
    %
    % Renderer_setup fields (with defaults):
    %   TDNumTaps      : FIR length (default 512)
    %   TDUseLinear    : true→use linearized phase ramp, false→atan() exact (default true)
    %   TDRefPoint     : [x y] sweet-spot in meters (default [0 0])
    %   TDRefRadiusR0  : 'auto' (mean distance of SSD from ref point) or scalar (default 'auto')

    properties
        fs                 % sampling rate (from virtual source)
        c                  % speed of sound
        N_taps = 2048       % FIR length
        useLinear = true   % linearized phase ramp (H_sp*_lin)
        refPoint = [0 0]   % sweet-spot position
        R0_override = []   % [] means auto
    end

    properties (SetAccess=protected)
        renderer_filter    % cell{Nls} of OLS_convolver
    end

    methods
        function obj = time_delay_renderer(virtual_source, loudspeaker_array, setup)
            obj = obj@base_renderer(virtual_source, loudspeaker_array);

            obj.fs = obj.virtual_source.source_signal.fs;
            obj.c  = audioapp.util.Physics.speedOfSound();

            % ----- read Renderer_setup (with safe defaults)
            RS = getfielddef(setup,'Renderer_setup', struct());
            obj.N_taps     = getfielddef(RS,'TDNumTaps',   512);
            obj.useLinear  = logical(getfielddef(RS,'TDUseLinear', true));
            obj.refPoint   = getfielddef(RS,'TDRefPoint', [0 0]);
            obj.R0_override = getfielddef(RS,'TDRefRadiusR0', 'auto');

            % ----- build filters once
            obj.build_or_update_filters();
        end

        function info = get_renderer_info(obj)
            lines = [
                "Time-delay stereo renderer — overview"
                ""
                "PARAMETERS"
                "• FIR taps (N):       " + string(obj.N_taps)
                "• Linear phase ramp:  " + string(obj.useLinear)
                "• Ref point [x y]:    [" + num2str(obj.refPoint(1),'%0.3f') + "  " + num2str(obj.refPoint(2),'%0.3f') + "]"
                "• R0 mode:            " + string(ifempty(obj.R0_override,'auto'))
                ""
                "NOTES"
                "• Picks two most frontal loudspeakers at the sweet spot and applies"
                "  opposite delays to steer towards the virtual source."
                "• If a single loudspeaker is perfectly aligned, it outputs an impulse there."
            ];
            info = strjoin(lines, newline);
        end

        % ------------------------------------------------------------------
        % Integration API
        % ------------------------------------------------------------------
        function update_renderer(obj, evt)
            %#ok<*INUSD>
            % Recompute the two-speaker filters (VS or SSD moved/rotated)
            obj.build_or_update_filters();
        end

        function update_renderer_settings(obj, RS)
            % Accept live changes from Setup.Renderer_setup
            needRebuild = false;

            if isfield(RS,'TDNumTaps') && isfinite(RS.TDNumTaps) && RS.TDNumTaps>0 ...
                    && RS.TDNumTaps ~= obj.N_taps
                obj.N_taps = RS.TDNumTaps;
                needRebuild = true;
            end
            if isfield(RS,'TDUseLinear')
                v = logical(RS.TDUseLinear);
                if v ~= obj.useLinear
                    obj.useLinear = v; needRebuild = true;
                end
            end
            if isfield(RS,'TDRefPoint') && numel(RS.TDRefPoint)>=2
                if any(RS.TDRefPoint(1:2) ~= obj.refPoint(1:2))
                    obj.refPoint = RS.TDRefPoint(1:2);
                    needRebuild = true;
                end
            end
            if isfield(RS,'TDRefRadiusR0')
                % allow [] / 'auto' / numeric
                v = RS.TDRefRadiusR0;
                if (ischar(v) || isstring(v)) && strcmpi(string(v),'auto')
                    if ~isempty(obj.R0_override), obj.R0_override = []; needRebuild = true; end
                elseif isempty(v)
                    if ~isempty(obj.R0_override), obj.R0_override = []; needRebuild = true; end
                elseif isnumeric(v) && isfinite(v) && v>0
                    if isempty(obj.R0_override) || obj.R0_override ~= v
                        obj.R0_override = v; needRebuild = true;
                    end
                end
            end

            if needRebuild
                obj.build_or_update_filters();
            end
        end

        function render(obj)
            x = obj.virtual_source.source_signal.time_series;
            for n = 1:numel(obj.renderer_filter)
                obj.output_signal{n}.set_signal( obj.renderer_filter{n}.convolve(x) );
            end
            obj.add_output_to_ssd_signal;
        end
    end

    % ======================================================================
    % Implementation
    % ======================================================================
    methods (Access=private)
        function build_or_update_filters(obj)
            % Collect geometry
            x0 = cell2mat(cellfun(@(s) s.position,    obj.ssd.loudspeakers,'uni',0)');
            n0 = cell2mat(cellfun(@(s) s.orientation, obj.ssd.loudspeakers,'uni',0)'); %#ok<NASGU> (kept for possible future use)
            xs = obj.virtual_source.position;

            % Compute all FIRs (Nx x Nls)
            H = obj.compute_td_filters(xs, x0);

            % Allocate/update per-speaker OLS convolver
            Nls = size(x0,1);
            if isempty(obj.renderer_filter) || numel(obj.renderer_filter) ~= Nls
                obj.renderer_filter = cell(1,Nls);
                for n = 1:Nls
                    obj.output_signal{n} = signal;  % ensure exists
                    obj.renderer_filter{n} = OLS_convolver(H(:,n), ...
                        length(obj.virtual_source.source_signal.time_series));
                end
            else
                for n = 1:Nls
                    obj.renderer_filter{n}.update_coefficients(H(:,n));
                end
            end
        end

        function td_filters = compute_td_filters(obj, xs, x0)
            % Core of your original get_td_filters(), generalized:
            % - ref point is obj.refPoint
            % - R0 can be overridden or auto = mean distance of SSD from ref

            xr = obj.refPoint(:)';                 % 1x2
            N  = obj.N_taps;
            w  = (0:N-1)'/N * 2*pi*obj.fs;         % Nx1 rad/s

            % Direction from ref→SSD (unit), projection on ref→source dir
            v  = bsxfun(@minus, xr, x0);                      % Nls x 2
            v  = bsxfun(@times, v, 1./sqrt(sum(v.^2,2)));     % unit
            dir_xs = (xr - xs);                               % 1x2
            dn = norm(dir_xs); if dn>0, dir_xs = dir_xs/dn; end
            vn = sum(v .* dir_xs, 2);                         % Nls x 1

            % Mean radius (ref→SSD), unless overridden
            if isempty(obj.R0_override) || (ischar(obj.R0_override) && strcmpi(string(obj.R0_override),'auto'))
                R0 = mean( sqrt(sum( (x0 - xr).^2, 2 )) );
            else
                R0 = obj.R0_override;
            end

            td_filters = zeros(N, size(x0,1));

            % Case A: a single loudspeaker perfectly aligned → pure impulse
            if max(abs(vn)) >= 1 - 1e-12
                [~,ix] = max(abs(vn));
                % IFFT{1} = delta at n=0 -> centered by fftshift
                td_filters(:,ix) = fftshift(ifft(ones(N,1), N, 'symmetric'));
                return;
            end

            % Case B: pick two most frontal loudspeakers
            [~,ixAll] = sort(vn,'descend');
            ixs = ixAll(1:2);

            % Base angle (between the two SSD rays from ref point), halved
            u1 = unitrow(x0(ixs(1),:) - xr);
            u2 = unitrow(x0(ixs(2),:) - xr);
            cosBase = max(-1,min(1, dot(u1,u2)));
            fi_base = acos(cosBase) / 2;

            % Base direction (bisector) and angle to the source direction
            base_vec = unitrow( (x0(ixs(1),:) + x0(ixs(2),:) - 2*xr) );
            vec_xs   = unitrow(xs - xr);
            cosF0    = max(-1,min(1, dot(base_vec, vec_xs)));
            fi0      = acos(cosF0);

            % Phase ramps
            if obj.useLinear
                % Linearized: φ = ± (1/4) * (ω/c) * R0 * tan(fi0)/tan(fi_base)
                phi = (1/4) * (w/obj.c) * R0 * (tan(fi0)/tan(fi_base));
            else
                % Exact: φ = ± atan( (ω/c) * R0 * tan(fi0)/tan(fi_base) )
                phi = atan( (w/obj.c) * R0 * (tan(fi0)/tan(fi_base)) );
            end
            H1 = exp( 1i*phi );
            H2 = exp(-1i*phi );

            td_filters(:,ixs(1)) = real( fftshift( ifft(H1, N, 'symmetric') ) );
            td_filters(:,ixs(2)) = real( fftshift( ifft(H2, N, 'symmetric') ) );
        end
    end
end

% ----- tiny helpers -------------------------------------------------------
function v = getfielddef(S, f, d)
if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function u = unitrow(v)
n = norm(v);
if n>0, u = v./n; else, u = [1 0]; end
end

function s = ifempty(val, repl)
if isempty(val), s = repl; else, s = val; end
end
