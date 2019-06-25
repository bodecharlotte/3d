classdef Spectral 
    properties
       % define the properties of the class here, (like fields of a struct)
           initial
           wave
           wold
           density
           k 
           w
           k_squared
           w_squared
           R
           K
           Rd
           Kd
       
    end
    methods
        function obj = Spectral(initial, P, Po, Co, m)
               obj.initial = initial;
               obj.wave = obj.initial;
               obj.wold = obj.wave;
               obj.density = real(conj(obj.wave).*obj.wave);
               obj.k = cat(2, 0:P.res/2 - 1, -P.res/2 : -1) * 2*pi/(P.interval_x(2)-P.interval_x(1));
               obj.w = cat(2, 0:P.res/2 - 1, -P.res/2 : -1) * 2*pi/(P.interval_y(2)-P.interval_y(1));  
               [obj.k,obj.w] = meshgrid(obj.k, obj.w); 
               obj.k_squared = obj.k.^2;
               obj.w_squared = obj.w.^2;
               obj.R = exp(-1i*(Po.pot + P.nonlinear*abs(obj.wave).^2)*P.timestep_imag/2);
               obj.K = exp(-1i * P.timestep_imag * (obj.k_squared+obj.w_squared)/2); 
               obj.Rd = exp(-1i*(Po.pot + P.nonlinear*abs(obj.wave).^2)*P.timestep/2);
               obj.Kd = exp(-1i * P.timestep * (obj.k_squared+obj.w_squared)/2); 
               
        end
    end
end

        