classdef Conservation
   % write a description of the class here.
       properties
       % define the properties of the class here, (like fields of a struct)
           Mass = 0
           En = 0
           Enold
           mu = 0
           xrms = 0
           yrms = 0
       end
       methods
       % methods, including the constructor are defined in this block
           function obj = Conservation(E, P, Po, method, der)
               obj = con(obj, E, P, Po, method, der);
           end
          
       end
   end