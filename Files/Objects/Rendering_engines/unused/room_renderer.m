classdef room_renderer < base_renderer

    properties
        receiver
        modalFEM_simulator
        rayTracer
        output_filters
        t0
        fs_modal
        freq_rt
        crossover_win
    end

    methods
        function obj = room_renderer(virtual_source, receiver, SSD, simulators)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.receiver = receiver;
            obj.modalFEM_simulator = simulators{1};
            obj.rayTracer = simulators{2};
            obj.t0 = 1;
            obj.fs_modal = ceil( max(obj.modalFEM_simulator.modes.modalOmega)/2/pi*4 );
            Nt = obj.t0*obj.fs_modal;
            freqModal = (0:Nt/2)'/Nt*obj.fs_modal;
            tria_ixs = find(obj.modalFEM_simulator.room.roomMesh.Elements(:,2) == 23);
            tetra_ixs = find(obj.modalFEM_simulator.room.roomMesh.Elements(:,2) == 34);

            obj.modalFEM_simulator.getDampingMatrix(freqModal, 'full');
            Nout = ceil(obj.virtual_source.source_signal.fs/obj.fs_modal*(Nt/2+1));
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_filters{n} = OLS_convolver(zeros(Nout,1), length(obj.virtual_source.source_signal.time_series));
            end

            fs_out = obj.virtual_source.source_signal.fs;
            Nt_rt = obj.t0*fs_out;
            obj.freq_rt = (0:Nt_rt/2)'/Nt_rt*fs_out;

            crossover_freq =  max(obj.modalFEM_simulator.modes.modalOmega)/2/pi;
            obj.crossover_win = zeros(size(obj.freq_rt,1),1);
            [~,fc_ix] = min(abs(crossover_freq-obj.freq_rt));
            Ncross = round(2*fc_ix/10)/2;
            Ncross = Ncross-mod(Ncross,2);
            win0 = hann(2*Ncross);
            fadeout = win0(end/2+1:end);

            obj.crossover_win(fc_ix-Ncross/2+1:fc_ix+Ncross/2) = fadeout;
            obj.crossover_win(1:fc_ix-Ncross/2) = 1;
            obj.update_renderer('boot');
        end

        function obj = update_renderer(obj, event)

            fs_out = obj.virtual_source.source_signal.fs;
            Nt_rt = obj.t0*fs_out;
            Nt = obj.t0*obj.fs_modal;
            freqModal = (0:Nt/2)'/Nt*obj.fs_modal;
            switch event
                case 'receiver_moved'
                    xr =  [obj.receiver.position, 1.7];
                    receiverOut = room_receiver(xr, obj.modalFEM_simulator.room.roomMesh );
                    obj.modalFEM_simulator.setReceiver( receiverOut );
                    obj.rayTracer.setReceiver(receiverOut);
                case 'virtual_source_moved'
                    xs = [obj.virtual_source.position,1.7];
                    sourceIn = room_source( xs, [1,0,0], 0, 'monopole',[],[] );
                    obj.modalFEM_simulator.setSource( sourceIn );
                    obj.rayTracer.setSource(sourceIn);
                    obj.rayTracer.traceRays(3);
                case 'boot'
                    xr =  [obj.receiver.position, 1.7];
                    receiverOut = room_receiver(xr, obj.modalFEM_simulator.room.roomMesh );
                    obj.modalFEM_simulator.setReceiver( receiverOut );
                    obj.rayTracer.setReceiver(receiverOut);
                    xs = [obj.virtual_source.position,1.7];
                    sourceIn = room_source( xs, [1,0,0], 0, 'monopole',[],[] );
                    obj.modalFEM_simulator.setSource( sourceIn );
                    obj.rayTracer.setSource(sourceIn);
                    obj.rayTracer.traceRays(3);
            end

            P_modal = conj( obj.modalFEM_simulator.evaluateHarmonicField( freqModal, 'simple' ) ).';
            P_fem = interp1(freqModal, abs(P_modal), obj.freq_rt,'linear',0).*...
                exp(1i*interp1(freqModal, unwrap(angle(P_modal)), obj.freq_rt,'linear',0))*fs_out/obj.fs_modal;
            P_rt = obj.rayTracer.evaluateRays(obj.freq_rt);

            P_out = obj.crossover_win.*P_fem + (1-obj.crossover_win).*P_rt;
            p = ifft( P_out,Nt_rt, 'symmetric' );
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_filters{n}.update_coefficients(p);
            end

        end

        function render(obj)
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n}.set_signal( obj.output_filters{n}.convolve(obj.virtual_source.source_signal.time_series) );
            end
            obj.add_output_to_ssd_signal;
        end


    end
end

