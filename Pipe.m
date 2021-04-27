% *** Class Pipe ***

classdef Pipe < handle
    
    %/////////////////////////////////////////////////////////////////////////////////////
    
    properties
        %---------------------------------------------------------------------------------
        length               double        % pipe length
        lengthAdjusted       double        % pipe length adjusted due to Courant Condition
        diameter             double        % pipe diameter
        roughness            double        % pipe roughness
        wavespeed            double        % pipe wavespeed
        %---------------------------------------------------------------------------------
        area                 double        % pipe cross section area
        f                    double        % Darcy-Weisbach friciton coeff
        fAdjusted            double        % adjusted f
        %---------------------------------------------------------------------------------
        initFlow             double        % pipe initial flow
        %---------------------------------------------------------------------------------
        transientFlow        cell          % pipe transient flow
        transientHead        cell          % pipe transient head
        FRMFlowUp            double
        FRMFlowDown          double
        %---------------------------------------------------------------------------------
    end
    
    %/////////////////////////////////////////////////////////////////////////////////////
    
    methods
        %---------------------------------------------------------------------------------
        function setArea(obj)
            obj.area = 0.25 * pi * obj.diameter.^2;
        end
        %---------------------------------------------------------------------------------
        
        %---------------------------------------------------------------------------------
    end
    
    %/////////////////////////////////////////////////////////////////////////////////////
    
end

