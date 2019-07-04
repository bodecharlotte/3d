classdef Parameters 
% write a description of the class here.
   properties
   % define the properties of the class here, (like fields of a struct)
       timestep;
       timestep_imag;
       maxiterations;
       interval_x;
       interval_y;
       interval_z;
       nonlinear;
       hinitial_x;
       hinitial_y;
       hinitial_z;
       h;
       a;
       x;
       y;
       z;
       X;
       Y;
       Z;
       res;
       hbar;
       mass;


   end
   methods
   % methods, including the constructor are defined in this block
       function obj = Parameters(timestep, maxIterations, interval_x, interval_y, interval_z, res, a, Co)
       % class constructor
           if(nargin > 0)
             obj.timestep = timestep;
             obj.timestep_imag = -1i*obj.timestep;
             obj.maxiterations   = maxIterations;
             obj.interval_x    = interval_x;
             obj.interval_y    = interval_y;
             obj.interval_z = interval_z;
             obj.nonlinear  = Co.g;
             obj.res   = res;
             obj.hinitial_x = (obj.interval_x(2)-obj.interval_x(1))/(obj.res);
             obj.hinitial_y = (obj.interval_y(2)-obj.interval_y(1))/(obj.res);
             obj.hinitial_z = (obj.interval_z(2)-obj.interval_z(1))/(obj.res);
             obj.h = [];
             obj.a = a;
             obj.x = (obj.interval_x(1)+obj.hinitial_x*(0:obj.res-1)); 
             obj.y = (obj.interval_y(1)+obj.hinitial_y*(0:obj.res-1)); 
             obj.z = (obj.interval_z(1)+obj.hinitial_z*(0:obj.res-1)); 
             obj.hbar = 1;
             obj.mass = 1;

           end
       end
   end
end