function [amp, phase_in] = get_amp_and_phase(sphere, input, mode)
amp = abs(input);

switch mode
    case 'normal'
        phase_in = unwrap(angle(input),[],3);
    case 'unwrap'
        phase_in = unwrap(angle(input),[],3);   
        phase_left = squeeze(phase_in(:,1,:));
        phase_right = squeeze(phase_in(:,2,:));
        [~,ix_min_l] = min(abs(phase_left*(0:size(phase_in,3)-1)'));
        [~,ix_min_r] = min(abs(phase_right*(0:size(phase_in,3)-1)'));

        %% Spiral tracking
        figure;sphere.plot_data(phase_left(:,200));
        shading interp
        [phase_out_l] = unwrap_sphere(sphere,phase_left,ix_min_l);
        [phase_out_r] = unwrap_sphere(sphere,phase_right,ix_min_r);
        phase = zeros(size(phase_in));
        phase(:,1,:) = phase_out_l;
        phase(:,2,:) = phase_out_r;
end

