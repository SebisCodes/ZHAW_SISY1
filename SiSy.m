
classdef SiSy
    %SISY Class for the ZHAW module SiSy1
    %   A class to analyze signals
    %   Created by Sebastian Grünwald, sebiscodes@gmail.com
    %   14.12.2023, Winterthur
    %   Github: https://github.com/SebisCodes/
    
    properties
        % Variables
        o_s (1,:) double % Array of signal values
        o_fs double % Frequency of your signal
        o_Tp double  % Time for one period
        o_Ts double % Time between two signal values
        o_x (1,:) uint64 % The x-array
        o_t (1,:) double % The time calculated by the x array and Ts
        o_N uint64 % The amount of values in the measurment
        o_Np uint64 % The amount of values in a period
        o_offset uint64 % The offset of the period part

        % Characteristics
        o_maxS double % Maximal value of the signal
        o_minS double % Minimum value of the signal
        o_rms double % RMS Value
        o_P double % Power of the function
        o_integral_t (1,:) double % The x-axis of the integral of the signal
        o_integral_s (1,:) double % The integral of the signal
        o_area double
        o_absarea double

        % Fourier
        fft_N uint64 % Amount of samples x(start:start+N-1);
        fft_start uint64 % Start value of array to create fft x(start:start+N-1)
        fft_s (1,:) double % The part of the function to create the fft
        fft_t (1,:) double % The part of the time value to create the fft
        fft_y (1,:) double % The fft function value
        fft_f (1,:) double % The fft frequency value

        % Convolution
        conv_t (1,:) double % The time vector of the convolution
        conv_s (1,:) double % The signal vector of the convolution

        % Fourier-Coefficients
        coeff_s (1,:) double % The part of the function which we get the fourier-coefficients from
        coeff_t (1,:) double % The part of the time value which we get the fourier-coefficients from
        coeff_Ak (1,:) double % Ak coefficients
        coeff_Bk (1,:) double % Bk coefficients
        coeff_Mk (1,:) double % Mk coefficients
        coeff_N uint64 % Amount of coefficients
        coeff_y (1,:) double % Approximated function
        coeff_yA (1,:) double % Approximated function using only A-coefficients
        coeff_yB (1,:) double % Approximated function using only B-coefficients
    end
    
    methods
        function obj = SiSy()
            %SISY Construct an instance of this class
        end
        
        function obj = addWav(obj, path, samplesPerPeriod, offset)
            arguments
                obj
                path string
                samplesPerPeriod uint64 = 0
                offset uint64 = 0
            end
            [obj.o_s,obj.o_fs] = audioread(path);
            obj = obj.setSignal(obj.o_s, obj.o_fs,samplesPerPeriod, offset);
        end

        %% Store a signal to this class
        function obj = setSignal(obj, signal_values,sampleFrequency, samplesPerPeriod, offset)
            arguments
                obj
                signal_values (1,:) string
                sampleFrequency double
                samplesPerPeriod uint64 = length(signal_values)
                offset uint64 = 0
            end
            % Signal calculations
            obj.o_offset = offset;
            obj.o_s = signal_values;
            obj.o_N = length(obj.o_s);
            obj.o_x = (0:(obj.o_N-1));
            obj.o_fs = sampleFrequency;
            obj.o_Ts = 1/obj.o_fs;
            obj.o_Np = samplesPerPeriod;
            obj.o_Tp = samplesPerPeriod*obj.o_Ts;
            obj.o_t = double(obj.o_x)*obj.o_Ts;
            if obj.o_Tp == 0
                obj.o_Tp = obj.o_t(1,end);
            end 

            % Get the characteristics of the function
            obj.o_minS = min(obj.o_s);
            obj.o_maxS = max(obj.o_s);
            [obj, obj.o_rms, obj.o_P] = obj.getRMSandP();
            [obj, obj.o_integral_t, obj.o_integral_s, obj.o_area, obj.o_absarea] = obj.getIntegral();
            [obj, obj.fft_t, obj.fft_s, obj.fft_f, obj.fft_y] = obj.getFFT();
        end

        %% Get rms and power of the signal by index offset
        % (use getIndexOffsetByTime to convert time to an offset)
        function [obj, rms, P] = getRMSandP(obj, offset, N)
            arguments
                obj
                offset uint64 = obj.o_offset % The offset to your signal
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.o_P = 0;
            for n = offset:offset+N-1
                obj.o_P = obj.o_P + obj.o_s(n)*obj.o_s(n);
            end
            obj.o_P = obj.o_P / obj.o_Tp;
            obj.o_rms = sqrt(obj.o_P);
            rms = obj.o_rms;
            P = obj.o_P;
        end

        %% Get the x and y of the signal
        function [x,y, fs, N] = getSignal(obj)
            x = obj.o_t;
            y = obj.o_s;
            fs = obj.o_fs;
            N = obj.o_N;
        end

        %% Get an integal of one period
        function [obj, t_integral, s_integral, area, absarea] = getIntegral(obj, offset, N)
            arguments
                obj
                offset uint64 = obj.o_offset % The offset to your signal
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            obj.o_area = 0;
            obj.o_absarea = 0;
            obj.o_integral_s = zeros(1,N);
            obj.o_integral_t = zeros(1,N);
            obj.o_integral_t(1) = obj.o_t(1);
            obj.o_integral_s(1) = 0;
            %obj.o_integral_t(1) = obj.o_t(1+offset);
            %obj.o_integral_s(1) = obj.o_s(1+offset);
            for n = 2:N
                obj.o_integral_t(n) = obj.o_t(n);
                n_offset = n+offset;
                addVal =  obj.o_s(n_offset-1)*(obj.o_t(n_offset)-obj.o_t(n_offset-1))+(obj.o_s(n_offset)-obj.o_s(n_offset-1))*(obj.o_t(n_offset)-obj.o_t(n_offset-1))/2;
                obj.o_integral_s(n) = obj.o_integral_s(n-1) + addVal;
                obj.o_area = obj.o_area + addVal;
                obj.o_absarea = obj.o_absarea + abs(addVal);
            end
            t_integral = obj.o_integral_t;
            s_integral = obj.o_integral_s;
            area = obj.o_area;
            absarea = obj.o_absarea;
        end

        %% Get the x and y of the signal by index offset 
        % (use getIndexOffsetByTime to convert time to an offset)
        function [xp,yp] = getPSignal(obj, offset, N)
            arguments
                obj
                offset uint64 = obj.o_offset % The offset to your signal
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            xp = obj.o_t(1+offset:1+offset+N);
            yp = obj.o_s(1+offset:1+offset+N);
        end
        
        %% Convert time offset to index offset
        function [offset] = getIndexOffsetByTime(obj, time)
            arguments
                obj
                time double % The time value of the beginning of your signal
            end
            offset = round(time/obj.o_Ts);
        end

        %% Make the fft out of the signal
        function [obj, fft_t, fft_s, fft_f, fft_y] = getFFT(obj, offset, N)
            arguments
                obj
                offset uint64 = obj.o_offset % The offset to your signal x(offset:offset+N-1)
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.fft_start = offset;
            obj.fft_N = N;
            obj.fft_t = obj.o_t(obj.fft_start:obj.fft_start+obj.fft_N-1);
            obj.fft_f = [0:obj.fft_N-1]*obj.o_fs/obj.fft_N;
            obj.fft_s = obj.o_s(obj.fft_start:obj.fft_start+obj.fft_N-1);
            obj.fft_y = abs(fft(obj.fft_s));
            fft_t = obj.fft_t;
            fft_s = obj.fft_s;
            fft_f = obj.fft_f;
            fft_y = obj.fft_y;
        end
    
        function [obj, conv_t, conv_s] = getConvolution(obj, h, offset, N)
            arguments
                obj
                h (1,:) double % The impulse response
                offset uint64 = obj.o_offset % The offset to your signal x(offset:offset+N-1)
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.conv_t = obj.o_t(offset:offset+N);
            obj.conv_s = conv(obj.o_s,h);
            obj.conv_s = obj.conv_s(offset:offset+N);
            conv_t = obj.conv_t;
            conv_s = obj.conv_s;
        end

        function [obj, coeff_t, coeff_s, approx_s, approx_sA, approx_sB, Ak, Bk, Mk] = getFourierCoefficients(obj, N_coefficients, offset, N)
            arguments
                obj
                N_coefficients uint64 = 11 % Amount of fourier-coefficients to calculate
                offset uint64 = obj.o_offset % The offset to your signal x(offset:offset+N-1)
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.coeff_N = N_coefficients;
            obj.coeff_t = obj.o_t(offset:offset+N-1);
            obj.coeff_s = obj.o_s(offset:offset+N-1);

            obj.coeff_Ak = zeros(1,obj.coeff_N);
            obj.coeff_Bk = zeros(1,obj.coeff_N);
            obj.coeff_Mk = zeros(1,obj.coeff_N);
            %obj.coeff_Ak(1) = mean(obj.coeff_s)*2;
            for k = 1:obj.coeff_N
                for n = 1:N
                    obj.coeff_Ak(k) = obj.coeff_Ak(k) + obj.coeff_s(n)*cos(2.0*pi*double((k-1)*n)/double(N));
                    obj.coeff_Bk(k) = obj.coeff_Bk(k) + obj.coeff_s(n)*sin(2.0*pi*double((k-1)*n)/double(N));
                end
                obj.coeff_Ak(k)=(2.0/double(N))*obj.coeff_Ak(k);
                obj.coeff_Bk(k)=(2.0/double(N))*obj.coeff_Bk(k);
                obj.coeff_Mk(k)=sqrt(obj.coeff_Ak(k)^2+obj.coeff_Bk(k)^2);
            end
            obj.coeff_y = zeros(1,N);
            obj.coeff_yA = zeros(1,N);
            obj.coeff_yB = zeros(1,N);

            obj.coeff_y = obj.coeff_Ak(1)*(1.0/2.0);
            obj.coeff_yA = obj.coeff_Ak(1)*(1.0/2.0);
            for n = 1:obj.coeff_N-1
                %obj.coeff_y = obj.coeff_y - obj.coeff_yA - obj.coeff_yB;
                %obj.coeff_yA = obj.coeff_yA+obj.coeff_Ak(n+1)*cos(2.0*pi*double(n)*obj.o_fs*obj.coeff_t);
                %obj.coeff_yB = obj.coeff_yB+obj.coeff_Bk(n+1)*sin(2.0*pi*double(n)*obj.o_fs*obj.coeff_t);
                Ak_temp = obj.coeff_Ak(n+1)*cos(pi/2.0*(1.0/obj.o_fs)*double(n)*(double(N))*obj.coeff_t);
                Bk_temp = obj.coeff_Bk(n+1)*sin(pi/2.0*(1.0/obj.o_fs)*double(n)*(double(N))*obj.coeff_t);
                obj.coeff_yA = obj.coeff_yA+Ak_temp;
                obj.coeff_yB = obj.coeff_yB+Bk_temp;
                obj.coeff_y = obj.coeff_y + Ak_temp+Bk_temp;
            end
            
            coeff_t = obj.coeff_t;
            coeff_s = obj.coeff_s;
            Ak = obj.coeff_Ak;
            Bk = obj.coeff_Bk;
            Mk = obj.coeff_Mk;
            approx_s = obj.coeff_y;
            approx_sA = obj.coeff_yA;
            approx_sB = obj.coeff_yB;
        end
    end
end

