classdef TissueLayer < handle
% Classdef for the tissueLayer object.  
%    This class is for use by the Tissue class to store details of the 
%    different tissue layers making up 
%    a tissue type.  For example, a single water based tissue may contain 
%    layers defining fibroglandular tissue, pectoralis muscle and lesions.  

    properties (SetAccess = {?TissueLayer, ?Tissue})

        %layerName: user defined tissue name.  Optional.
        %
        layerName = '';
        
        %mask: Mask defining tissue extent in phantom object.  Can be
        %logical or values can range from 0 to 1.
        %mask = [];
        mask = qt_image.empty(1,0);
        
        % t10: Nominal T1 relaxation time for the tissue in ms
        t10 = []; 
        
        %Optional parameters:
        
        %t2: T2 relaxation time for the tissue in ms
        t2 = [];
        
        %modelParam: Enhancement model parameters.  The choice of these
        %parameters will depend upon the model chosen by the user for
        %enhancement simulation.  
        modelParams = struct;
        

    end
    
    properties (SetAccess = 'protected', Dependent)
        
        %Indicaites if tissue layer has all components filled in
        %tissueLayerValid;

    end
    
    properties (SetAccess = 'protected')
        
        defaultTissueTypes={'blood', 'fat', 'fibroglandular'};
        
    end
    
    % ----------- Class Constructor -------------%
    methods
        function obj = TissueLayer(varargin)
            
            % Do nothing with zero inputs...this is required by MATLAB for
            % smooth operation
            if ~nargin
                return
            end
            
            %First input must be the qt_image mask object:
            if isa(varargin{1}, 'qt_image')
                if islogical(varargin{1}.value)
                    obj.mask = qt_image(varargin{1}.value);
                else
                    obj.mask = varargin{1};
                end
            else
                ERRMSG = ['Error: Expecting  first input to TissueLayer constructor to be of class qt_image'];
                error(ERRMSG);
            end
            
            the_properties = properties(TissueLayer);
            counter = 2;
            while counter <= size(varargin, 2)
                thevar = varargin{counter};
                foundProp = false;   
                %If the input is a string, figure out what property it
                %corresponds to and call that set function:
                if ischar(thevar) 
                    for prop_num = 1:length(the_properties)
                        if strcmp('modelParams', thevar)
                            obj.setModelParams(varargin{counter + 1}, varargin{counter + 2});
                            counter = counter +1;
                            break;
                        elseif strcmp(the_properties{prop_num}, thevar)
                            obj.(thevar) = varargin{counter + 1};
                            foundProp = true;
                            break;
                        end
                    end
                end
                
                
                if foundProp
                    counter = counter + 2;
                else
                    counter = counter + 1;
                end
            end
            
            if isempty(obj.layerName)
                obj.layerName = 'No Tissue Name Assigned';
            end

        end % End Class Constructor
    end  
    
    %---------------Set Methods ----------------------%
    methods 
        function set.layerName(obj, str)
        %layerName: user defined tissue layer name.  Optional.
        %
            validateattributes(str, {'char'}, {'2d'});
            obj.layerName = str;
        end %Tissue.set.layerName
               
        function set.t10(obj, value)
        % t10: Nominal T1 relaxation time for the tissue in ms
        
            if isempty(value)
                obj.t10 = [];
                return;
            end
        
            if ischar(value)
                obj.setDefaultTissue('t10', value);
                return;
            end
            
            %Validate the input:
            validateattributes(value,{'double'}, {'nonempty', 'size', [1, 1], 'nonnegative'}); 
            obj.t10 = value;

        end %Tissue.set.t10
        
        function set.t2(obj, value)
        % t2: Nominal T2 relaxation time for the tissue in ms
        
            if isempty(value)
                obj.t2 = [];
                return;
            end
        
            if ischar(value)
                obj.setDefaultTissue('t10', value);
                return;
            end
            
            %Validate the input:
            validateattributes(value,{'double'}, {'nonempty', 'size', [1, 1], 'nonnegative'}); 
            obj.t2 = value;

        end %Tissue.set.t2 
        
%         function set.modelParams(obj, ~)
%       
%             ERRMSG = ['Error setting modelParams: Use the function obj.setModelParams function \n'...
%                         'to change the ehnhancement model parameters for a given tissue layer'];
%             error(ERRMSG, ' ');
% 
%         end
        
    end
    
    %---------------get Methods ----------------------%
    methods
        
%        function validValue = get.tissueLayerValid(obj)
%             
%             %A valid TissueLayer must have a mask and t10 value
%             if isempty(obj.mask) 
%                 validValue = false;
%                 return;
%             elseif isempty(obj.t10)
%                 validValue = false;
%                 return;
%             else
%                 % Either all pharmacokinetic values must be set or empty
%                 if isempty(obj.ktrans) && isempty(obj.ve) && isempty(obj.vp)
%                     validValue = true;
%                 elseif ~isempty(obj.ktrans) && ~isempty(obj.ve) && ~isempty(obj.vp)
%                     validValue = true;
%                 else
%                     validValue = false;
%                 end
%             end
%             
%         end %End obj.set.tissueLayerValid 
        
    end
    
    %---------------Other Methods ----------------------%
    methods (Access = {?Tissue})
        
        function setModelParams(obj, theField, theValue)
        % modelParam: parameters for the enhancement model chosen for 
        % simulation. Values must be a structure containing field names 
        % and values for simulation model.      
            
            %Validate the field name
            if ~isvarname(theField)
                ERRMSG = ['ERROR setting model parameters: the chosen field name is not valid'];
                error(ERRMSG);
            end
            
            obj.modelParams.(theField) = theValue;
            
        end
        
        function removeModelParams(obj, theField)
            obj.modelParams = rmfield(obj.modelParams.(theField));
        end
        
    end
    
    %---------------Restricted Methods ----------------------%
    methods(Access = 'protected')
        
        
        function obj = setDefaultTissue(obj, varargin) 
        % [ OBJ ] = setDefaultTissue(OBJ, STR) will set up the default
        % values for T1 and T2 relaxation times as well as ktrans, ve and
        % vp as appropriate.  The input OBJ is the Tissue object that is
        % being modified.  STR specifies the tissue type.  Options are:
        %
        %  'blood': Sets t10 = 1440 ms, ktrans = 0 mL/g/min,  ve = 1,  vp = 0.55.  
        %  'fat' :  Sets t10 = 296 ms,  ktrans = [], ve = [], vp = [].
        %  'fibroglandular': Sets t10 = 1266 ms, ktrans = 5e-4 mL/g/min; ve =0.3, and vp = 0 
        % 
        % Note: In the settings for the blood, vp setting chosen to account 
        % for assumptions made in the General Kinetic Modeling function in
        % the digital phantom.  Setting vp = 0.55 allows the vascular input
        % function to come through in the selected tissues unaltered.
        %
        % Note: For tissues that do not show Gad uptake like fat, ktrans,
        % ve and vp are left empty.
        %
        
            %Default Tissue Properties:
            %Blood:
            default(1).name = obj.defaultTissueTypes{1}; 
            default(1).t10 = 1440;
            default(1).t2 = 250;
            default(1).ktrans = 0;
            default(1).ve = 1;
            default(1).vp = 0.55;
            
            %Fat:
            default(2).name = obj.defaultTissueTypes{2}; 
            default(2).t10 = 296;
            default(2).t2 = 53;
            default(2).ktrans = [];
            default(2).ve = [];
            default(2).vp = [];
            
            %Fibroglandular
            default(3).name = obj.defaultTissueTypes{3}; 
            default(3).t10 = 1266;
            default(3).t2 = 57;
            default(3).ktrans  = 0.03/60;
            default(3).ve = 0.3;
            default(3).vp = 0;
            
        
            %Figure out what tissue type is being called:
            %if only 1 input, this is calling for setting all tissue
            %properties to the default tissue type:
            if length(varargin) == 2
                str = lower(varargin{2});
                validatestring(str, obj.defaultTissueTypes);
                
                for counter = 1:length(obj.defaultTissueTypes)
                    if strcmp(str, obj.defaultTissueTypes{counter})
                        obj.t10 = default(counter).t10;
                        obj.t2 = default(counter).t2;
                        %obj.ktrans = default(counter).ktrans;
                        %obj.ve = default(counter).ve;
                        %obj.vp = default(counter).vp;
                        break;
                    end
                end
            
            %we must have a single property or list of properties to set 
            %with the default value    
            else 
                the_properties = {'t10', 't2'};%{'t10', 't2', 'ktrans', 've', 'vp'};
                counter = 1;
                while counter <= length(varargin) - 1
                    for counter2 = 1:length(the_properties)
                        if strcmp(varargin{counter}, the_properties{counter2})

                            %Check the next variable for a matching
                            %default tissue type:
                            for counter3 = 1:length(obj.defaultTissueTypes)
                                if strcmp(varargin{counter + 1}, obj.defaultTissueTypes{counter3})
                                    obj.(varargin{counter}) = default(counter3).(varargin{counter});
                                    counter = counter + 1;
                                    break;
                                end
                            end %end defaultTissueTypes loop
                        end 
                    end %End properties loop
                    counter = counter + 1;
                end %End input loop 
            end   

        end %end setDefaultTissue
    
    end

end