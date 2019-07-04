classdef Energy
   % write a description of the class here.
       properties
       % define the properties of the class here, (like fields of a struct)
           Kinetic = 0
           PotentialE = 0
           Nonlinear = 0
           initial
           wave 
           wold
           density 
           matrix
           enmatrix
           enmu
           k
           w
           v
           k_squared
           w_squared
           v_squared
       end
       methods
       % methods, including the constructor are defined in this block
           function obj = Energy(initial, P)
               obj.initial = initial;
               obj.wave = obj.initial;
               obj.wold = obj.wave;
               obj.density = real(conj(obj.wave).*obj.wave);
               obj.k = cat(2, 0:P.res/2 - 1, -P.res/2 : -1) * pi/(P.interval_x(2));
               obj.w = cat(2, 0:P.res/2 - 1 , -P.res/2 : -1) * pi/(P.interval_y(2));
               obj.v = cat(2, 0:P.res/2 - 1 , -P.res/2 : -1) * pi/(P.interval_z(2));
               %fftshift(-P.res/2:P.res/2-1*2*pi/P.interval(2));
               obj.k_squared = obj.k.^2;
               obj.w_squared = obj.w.^2;
               obj.v_squared = obj.v.^2;
           end
          
       end
   end