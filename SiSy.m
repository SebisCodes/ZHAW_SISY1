
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
    end
    
    methods
        function obj = SiSy()
            %SISY Construct an instance of this class
        end
        
        function obj = addWav(obj, path, time_of_one_period, offset)
            arguments
                obj
                path string
                time_of_one_period double = 0
                offset uint64 = 0
            end
            [obj.o_s,obj.o_fs] = audioread(path);
            obj = obj.setSignal(obj.o_s, obj.o_fs,time_of_one_period, offset);
        end

        %% Store a signal to this class
        function obj = setSignal(obj, signal_values,frequency, time_of_one_period, offset)
            arguments
                obj
                signal_values (1,:) string
                frequency double
                time_of_one_period double = 0
                offset uint64 = 0
            end
            % Signal calculations
            obj.o_s = signal_values;
            obj.o_N = length(obj.o_s);
            obj.o_x = (0:(obj.o_N-1));
            obj.o_fs = frequency;
            obj.o_Ts = 1/obj.o_fs;
            obj.o_Tp = time_of_one_period;
            obj.o_t = double(obj.o_x)*obj.o_Ts;
            if obj.o_Tp == 0
                obj.o_Tp = obj.o_t(1,end);
            end 

            obj.o_Np = obj.o_Tp/obj.o_Ts;

            % Get the characteristics of the function
            obj.o_minS = min(obj.o_s);
            obj.o_maxS = max(obj.o_s);
            [obj, obj.o_rms, obj.o_P] = obj.getRMSandP(offset);
            [obj, obj.o_integral_t, obj.o_integral_s, obj.o_area, obj.o_absarea] = obj.getIntegral(offset);
            [obj, obj.fft_t, obj.fft_s, obj.fft_f, obj.fft_y] = obj.getFFT(offset);
        end

        %% Get rms and power of the signal by index offset
        % (use getIndexOffsetByTime to convert time to an offset)
        function [obj, rms, P] = getRMSandP(obj, offset, N)
            arguments
                obj
                offset uint64 % The offset to your signal
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.o_P = 0;
            for n = offset+1:offset+N
                obj.o_P = obj.o_P + obj.o_s(n)*obj.o_s(n);
            end
            obj.o_P = obj.o_P / obj.o_Tp;
            obj.o_rms = sqrt(obj.o_P);
            rms = obj.o_rms;
            P = obj.o_P;
        end

        %% Get the x and y of the signal
        function [x,y, f, N] = getSignal(obj)
            x = obj.o_t;
            y = obj.o_s;
            f = obj.o_fs;
            N = obj.o_N;
        end

        %% Get an integal of one period
        function [obj, t_integral, s_integral, area, absarea] = getIntegral(obj, offset, N)
            arguments
                obj
                offset uint64 = 0 % The offset to your signal
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            obj.o_area = 0;
            obj.o_absarea = 0;
            obj.o_integral_s = obj.o_s(1:N)*0;
            obj.o_integral_t = obj.o_s(1:N)*0;
            %obj.o_integral_t(1) = obj.o_t(1+offset);
            %obj.o_integral_s(1) = obj.o_s(1+offset);
            for n = 2+offset:obj.o_N+offset
                obj.o_integral_t(n-offset) = obj.o_t(n);
                addVal =  obj.o_s(n-1)*(obj.o_t(n)-obj.o_t(n-1))+(obj.o_s(n)-obj.o_s(n-1))*(obj.o_t(n)-obj.o_t(n-1))/2;
                obj.o_integral_s(n-offset) = obj.o_integral_s(n-1) + addVal;
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
                offset uint64 = 0 % The offset to your signal
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
                offset uint64 = 0 % The offset to your signal x(offset:offset+N-1)
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.fft_start = offset;
            obj.fft_N = N;
            obj.fft_t = obj.o_t(obj.fft_start:obj.fft_start+obj.fft_N);
            obj.fft_f = [0:obj.fft_N]*obj.o_fs/obj.fft_N;
            obj.fft_s = obj.o_s(obj.fft_start:obj.fft_start+obj.fft_N);
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
                offset uint64 = 0 % The offset to your signal x(offset:offset+N-1)
                N uint64 = obj.o_Np % The amount of samples from the offset x(offset:offset+N-1)
            end
            offset = offset + 1;
            obj.conv_t = obj.o_t(offset:offset+N);
            obj.conv_s = conv(obj.o_s,h);
            obj.conv_s = obj.conv_s(offset:offset+N);
            conv_t = obj.conv_t;
            conv_s = obj.conv_s;
        end
    end
end

