classdef Tissue
% classdef for tissue object that contains the image space representation 
% of individual tissues for the the object to be simulated as well as 
% properties needed by the digital phantom simulator.  
%test change
    
    %Properties
    properties
        
        % tissueType defines the general class of the tissue.  Current
        % options are:
        %    fat
        %    water
        %    custom
        %
        tissueType = '';
        
        %Base T1 relaxation time for the tissue contained in s0.  If not
        %defined, value will default to 1266ms, the T1 of normal
        %fibroglandular tissue.  Regions with different T1s can be defined
        %by adding a tissue layer and assigning a different T1 value.
        defaultT10 = [];
        
        tissueName = '';
        
        xyzRes = [];

    end % End properties
    
    
    properties (SetAccess = 'protected')
        
        % noiseMask contains masks that allow the tissue to have different 
        % behaviors in different regions.  For example, the ‘Water’ type 
        % tissue may have regions defined as fibroglandular tissue and 
        % pectoralis muscle.
        tissueLayers = {};
        
        %Information for scignal precession:
        signalPrecession = struct;
        
        %customizable user parameter array:
        userParams = struct;
        
    end
    
    properties (SetAccess = 'protected', Dependent)
        % s0 is the base signal intensity for the tissue, ideally for a
        % proton density weighted image (independent of T1 and T2 decay).
        % Data may be magnitude data or complex.    This data is typically
        % obtained from fat/water separated images of a volunteer that have
        % been filtered and read into a matrix.
        %s0 = [];
        %
        s0 = [];
        
        % noiseMask contains values ranging from 0 to 1 that define the location of 
        % the Tissue.  This defines what fraction of each voxel belongs to this 
        % tissue type. This allows for exclusion of voxels in noise regions from 
        % dynamic processing and allows for partial voluming of different tissue 
        % types such as between fat and water.  If not set explicitly by the user 
        % values are assumed to be 1.  
        noiseMask = [];
        
    
    end
    
        
    properties(Dependent)
        
        %xyzDims indicates the size in the x, y and z directions of the s0
        % images.  All subsequently added masks must match this size.
        xyzDims=[];
        
        
        %layerNames: assigned names of the different tissue layers
        layerNames = {};
        
        t10 = [];
        t2 = [];
        
        modelParams = struct;
        
        
    end
    
    
    properties (SetAccess = 'protected', Hidden)

         noiseMask_qtimage = qt_image.empty(1,0);
         
         s0_median = [];
         s0_mean = [];
         s0_threshold = [];

    end  
    
    properties (SetAccess = {?Tissue, ?Phantom}, Hidden)
        
        s0_qtimage = qt_image.empty(1,0);
        s0_qtimageImag = qt_image.empty(1,0);
        
        complex_data = false;
        
    end
    
    properties (Access = 'protected', Hidden)
       
        %Allows for reversible adding and removing of inserted tissue
        %layers
        insertTissue_pipeline = {};
        
        linkedLayers = {};
        
    end
    
    properties (Dependent, Hidden)
        
        Dims =[];
        
    end
      
    
     % ----------------- Class Constructor ---------------------------%
    methods
        function obj = Tissue(theS0, varargin)
        % Tissue constructs an instance of the Tissue class
        % 
        %    OBJ = Tissue(IMAGE_DATA) creates a TISSUE object OBJ by reading in the 
        %    specified image data.  IMAGE_DATA can be:
        %       *String specifying a specific file or files to be read in.  
        %        Metadata will be obtained from the    header if possible.
        %       *Filepath to a directory containing image files.  All images in the 
        %        directory will be read in.  Images in subdirectories will be 
        %        ignored.
        %       *A single 2D or 3D numeric Matlab array containing the desired 
        %        image information.  
        % 
        %    OBJ = Tissue(IMAGE_DATA, MASK_DATA) creates a TISSUE object OBJ from 
        %    the information referenced by IMAGE_DATA and also creates the 
        %    corresponding tissue mask noiseMask. where the tissue mask defines 
        %    where the tissue is present (mask = 1) or not (mask = 0, i.e. noise or 
        %    other tissue types).  MASK_DATA can be of the same format as 
        %    TISSUE_DATA but values must range between 0 and 1.  If
        % 
        %    OBJ = Tissue(…, ‘PROP1’, VAL1, …) creates the TISSUE object as 
        %    described above but also sets the specified properties.
        
            % Do nothing with zero inputs...
            if ~nargin
                return
            end
        
            %Send the First input to qt_image:
            if isnumeric(theS0) && ~isreal(theS0)
                obj.complex_data = true;
                obj.s0_qtimage = qt_image(real(theS0));
                obj.s0_qtimageImag = qt_image(imag(theS0));
            elseif isa(theS0, 'qt_image')
                obj.s0_qtimage = theS0;
            else
                obj.s0_qtimage = qt_image(theS0);
            end
                
            %If the second input is not a string, set it the mask:
            if length(varargin) > 0 && ~ischar(varargin{1})
                obj.noiseMask = varargin{1};
            end
            
            %Set up default properties:
            thelayerName = 'Default nonenhancing Tissue';
            theT10 = [];
            theT2 = [];
            theModelParams = {};
            theppm = [];
            theRelAmp = [];
            
            %Assign other properties as appropriate:
            the_properties = properties(Tissue);
            counter = 1;
            while counter <= size(varargin, 2)-1
                thevar = varargin{counter};
                foundProp = false;   
                %If the input is a string, figure out what property it
                %corresponds to and call that set function:
                if ischar(thevar)
                    
                    %Catch setting of certain properties in constructor:
                    if strcmp (thevar, 'layerName')
                        foundProp = true;
                        if ischar(varargin{counter + 1})
                            thelayerName = varargin{counter + 1};
                        end
                        break;
                    elseif strcmp (thevar, 't10')
                        foundProp = true;
                        theT10 = varargin{counter + 1};
                        break;
                    elseif strcmp(thevar, 't2')
                        foundProp = true;
                        theT2 = varargin{counter + 1};
                        break;
                    elseif strcmp(thevar, 'modelParams')
                        foundProp = true;
                        paramNum = length(theModelParams + 1);
                        theModelParams(paramNum, 1) = varargin{counter + 1};
                        theModelParams(paramNum, 2) = varargin{counter + 2};
                        counter = counter + 1;
                        break;
                    elseif strcmp(thevar, 'ppm')
                        foundProp = true;
                        theppm = varargin{counter + 1};
                        break;
                    elseif strcmp(thevar, 'relAmp')
                        foundProp = tru; 
                        theRelAmp = varargin{counter + 1};
                        break;
                    else
                        %identify the property:
                        for prop_num = 1:length(the_properties)
                            if strcmp(the_properties{prop_num}, thevar)
                                obj.(thevar) = varargin{counter + 1};
                                foundProp = true;
                                break;
                            end
                        end
                    end
                end
                
                %If the second input variable is a string but not a string
                %corresponding to one of the properties, try using it to
                %set the noiseMask:
                if ~foundProp && counter == 1
                    try
                        obj.noiseMask = varargin{1};
                    catch
                    end   
                end
                
                
                if foundProp
                    counter = counter + 2;
                else
                    counter = counter + 1;
                end

            end
            
            %If no mask was set, set a default mask with all ones:
            if isempty(obj.noiseMask)
                obj.noiseMask = ones(obj.s0_qtimage.dimSize);
            end
            
            %If the tissue type was not set, set it to water: 
            if isempty(obj.tissueType)
                obj.tissueType = 'water';
            end
            
            %Set the species precession settings:
            if ~isempty(theppm) || ~isempty(theRelAmp)
                obj = obj.setPrecession(theppm, theRelAmp);
            elseif strcmp(obj.tissueType, 'fat')
                %theppm = -[3.2879,    2.4894,   -0.7359,    3.6950,    1.8318,    0.3601];
                %theRelAmp = [0.6200,    0.1500,    0.1000,    0.0600,    0.0300,    0.0400];
                
                theFatModel = fatModel('ModelType', 'SubFat');
                theppm = theFatModel.Model_Relative_ppm;
                theRelAmp = theFatModel.Model_Relative_Amplitude;
                obj = obj.setPrecession(theppm, theRelAmp);
                
            else
                theppm = [0];
                theRelAmp = [1];
                obj = obj.setPrecession(theppm, theRelAmp);
            end
            
            %If no T1 value was set by the user, set the default value:
            if isempty(obj.defaultT10) && isempty(theT10)
                if strcmp(obj.tissueType, 'water')
                    obj.defaultT10 = 1266; %T1 for Fibroglandular Tissue
                elseif strcmp(obj.tissueType, 'fat')
                    obj.defaultT10 = 296; %T1 for Fat
                else
                    ERRMSG = ['ERROR: No default t10 for the tissue type %s'];
                    error(ERRMSG, obj.tissueType);
                end
            elseif isempty(obj.defaultT10)
                obj.defaultT10 = theT10;
            end
            
            if isempty(theT10)
                theT10 = obj.defaultT10;
            end
            
            %Set up the first tissue layer as nonenhancing tissue:
            obj.tissueLayers(1) = {TissueLayer(obj.noiseMask_qtimage, 'layerName', thelayerName)};
            obj.insertTissue_pipeline(1) = {[]};
            
            %Add T1 and T2:
            obj = obj.setT10(1, theT10);
            obj = obj.setT2(1, theT2);
            
            %Add object resolution:
            if isempty(obj.xyzRes)
                obj.xyzRes = obj.s0_qtimage.elementSpacing;
            end
            
            %Add any assigned model parameters:
%             for counter = 1:length(theModelParams)
%                 obj = obj.setModelParams(theModelParams{counter, 1}, theModelParams{counter, 2});
%             end

            %Record the median value excluding the noise:
            if obj.complex_data == false
                %the_values = obj.s0(:,:,round(obj.xyzDims(3)/4):obj.xyzDims(3) - round(obj.xyzDims(3)/4));
                the_values = obj.s0(obj.noiseMask == 1);
                level = graythresh(1/obj.s0_qtimage.elementMax * the_values)*obj.s0_qtimage.elementMax;
            else
                %the_values = abs(obj.s0(:,:,round(obj.xyzDims(3)/4):obj.xyzDims(3) - round(obj.xyzDims(3)/4)));
                the_values = abs(obj.s0(obj.noiseMask == 1));
                the_max = max(the_values(:));
                level = graythresh(1/the_max * the_values)*the_max;
            end
            tissue_values = the_values(the_values > level);
            obj.s0_median = median(tissue_values);
            obj.s0_mean = mean(tissue_values);
            obj.s0_threshold = level;
                

        end % End Tissue constructor
    end    
    
    % ----------------- Set Methods ---------------------------%
    methods
        function obj = set.tissueName(obj, str)
        %tissueName: user defined tissue layer name.  Optional.
        %
            validateattributes(str, {'char'}, {'2d'});
            obj.tissueName = str;
        end %Tissue.set.layerName
        
        function obj = set.tissueType(obj, str)
%         tissueType defines the general class of the tissue.  Current
%         options are:
%            fat
%            water
%         
%             if isempty(str)%allows tissue type to be cleared
%                 obj.tissueType = '';
%             end



            
            str = lower(str);
            valid_strings = {'fat', 'water', 'custom'};
            
            validatestring(str, valid_strings);
            obj.tissueType = str;
            
        end %Tissue.set.tissueType
        
        function obj = set.noiseMask(obj, variable)
            
            %use the input to set the qt_image object:
            if islogical(variable)
               %temp_object = qt_image( double(variable), 'elementClass','logical' );
               temp_object = qt_image( double(variable));
            else
                temp_object = qt_image(variable);
            end
                
            %Verify values are in the proper range:
            if temp_object.elementMin < 0 || temp_object.elementMax > 1
                ERRMSG=(['ERROR setting noiseMask: The values in the mask object are expected to range \n'...
                        'from 0 to 1.  Values in the input range from %d to %d.']);
                
                error(ERRMSG, temp_object.elementMin, temp_object.elementMax);
            end
            
            %Verify the object is of the proper size:
            if max(temp_object.dimSize ~= obj.s0_qtimage.dimSize)
                ERRMSG=(['Error setting noiseMask: Expected the mask to be the same size \n'...
                         'as the s0 object (Tissue images).  The size of s0 is: [%d, %d, %d] \n' ...
                         'while the size of the current object is [%d, %d, %d]']);
                
                 error(ERRMSG, obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3), temp_object.dimSize(1), temp_object.dimSize(2), temp_object.dimSize(3));
            end

            obj.noiseMask_qtimage = temp_object;
            
        end  %End Tissue.set.noiseMask
        
        function obj = set.xyzDims(obj, thesize)
            
            %Validate the requested size:
            if ~isequal(ceil(thesize), floor(thesize))
                ERRMSG = ['Error using set.xyzDims.  Expected input values to be integers'];
                error(ERRMSG);
            end
            validateattributes(thesize, {'numeric'}, {'numel', 3, 'positive', 'nonnan', 'nonempty', 'vector', 'even'});
            
            %Check to see if the size is actually being changed:
            if isequal(thesize, obj.xyzDims)
                return;
            end
            
            original_dims = obj.xyzDims;
            
            %Call the resize function on all qt_image objects:
            obj.s0_qtimage.resize([thesize(2), thesize(1), thesize(3)]);
            obj.noiseMask_qtimage.resize([thesize(2), thesize(1), thesize(3)]);
            
            if obj.complex_data
                obj.s0_qtimageImag.resize([thesize(2), thesize(1), thesize(3)]);
            end
            
            for counter = 1:length(obj.tissueLayers)
                obj.tissueLayers{counter}.mask.resize([thesize(2), thesize(1), thesize(3)]);
            end
            
            obj.xyzRes = obj.xyzRes .* (original_dims)./obj.xyzDims;
            
        end %End Tissue.set.xyzDims
        
        function obj = set.xyzRes(obj, theres)
            
            validateattributes(theres, {'numeric'}, {'numel', 3, 'positive', 'nonnan', 'nonempty', 'vector'});
            
            %Set the element spacing for the S0 qt_image object:
            obj.xyzRes = theres;
            
        end
        
        function obj = set.layerNames(obj, ~)
            
            ERRMSG = (['Error trying to set tissue layerNames.  Use: \n'...
                       '   OBJ = OBJ.setLayerName(NUMBER, NAME)%s']);
                  
            error(ERRMSG, ' ');

        end %End Tissue.set.layerNames
        
        function obj = set.t10(obj, ~)
            
            ERRMSG = (['Error trying to set tissue layer t10.  Use: \n'...
                       '   OBJ = OBJ.setT10(NUMBER, T10)%s']);
                  
            error(ERRMSG, ' ');

        end %End Tissue.set.t10
        
        function obj = set.t2(obj, ~)
            
            ERRMSG = (['Error trying to set tissue layer t10.  Use: \n'...
                       '   OBJ = OBJ.setT2(NUMBER, T2)%s']);
                  
            error(ERRMSG, ' ');

        end %End Tissue.set.t2
                 
        function obj = set.defaultT10(obj, value)
        % defaultT10: Nominal T1 relaxation time for the tissue in ms
            
            %Validate the input:
            validateattributes(value,{'double'}, {'nonempty', 'size', [1, 1], 'nonnegative'}); 
            obj.defaultT10 = value;
           
            if ~isempty(obj.tissueLayers) && ~isempty(obj.tissueLayers{1})
                obj.tissueLayers{1}.t10 = value;
            end

        end %Tissue.set.defaultT10
        
        function obj = set.modelParams(obj, ~)
           
            ERRMSG = ['Error setting modelParams: Use the function obj.setModelParams function \n'...
                      'to change the ehnhancement model parameters for a given tissue layer'];
            error(ERRMSG, ' ');
            
        end


    end %End Set methods
    
    % ----------------- Get Methods ---------------------------%
    methods
        
        function the_size = get.xyzDims(obj)
            if isempty(obj.s0_qtimage)
                the_size = [0, 0, 0];
            else
                the_size = [obj.s0_qtimage.dimSize(2), obj.s0_qtimage.dimSize(1), obj.s0_qtimage.dimSize(3)];
            end
        end %end Tissue.get.xyzDims
        
        function theDims = get.Dims(obj)
           if isempty(obj.s0_qtimage)
               theDims = [0, 0, 0];
           else
               theDims = obj.s0_qtimage.dimSize; 
           end
        end
        
        function theLayerNames = get.layerNames(obj)
            
            theLayerNames = cell([1,length(obj.tissueLayers)]);
            
            for counter = 1:length(obj.tissueLayers)
               theLayerNames(counter) = {obj.tissueLayers{counter}.layerName};
            end

        end %End Tissue.get.layerNames
        
        function theS0 = get.s0(obj)
           
            if obj.complex_data
                theS0 = complex(obj.s0_qtimage.value, obj.s0_qtimageImag.value);
            else
                if isempty(obj.s0_qtimage)
                    theS0 = [];
                else
                    theS0 = obj.s0_qtimage.value;
                end
            end
            
        end
        
        function theNoiseMask = get.noiseMask(obj)
            
            if isempty(obj.noiseMask_qtimage)
                theNoiseMask = [];
            else
                theNoiseMask = obj.noiseMask_qtimage.value;
            end
            
        end
        
        function theT10s = get.t10(obj)
            
            theT10s = zeros([1,length(obj.tissueLayers)]);
            
            for counter = 1:length(obj.tissueLayers)
               theT10s(counter) = obj.tissueLayers{counter}.t10;
            end

        end %End Tissue.get.t10
        
        function theT2s = get.t2(obj)
            
            theT2s = cell([1,length(obj.tissueLayers)]);
            
            for counter = 1:length(obj.tissueLayers)
               theT2s(counter) = {obj.tissueLayers{counter}.t2};
            end

        end %End Tissue.get.t2
             
        function theModelParams = get.modelParams(obj)
                      
            theModelParams = struct;
            
            if isempty(obj.tissueLayers)
                return;
            end

            %Get the field names:
            fieldNames = fieldnames(obj.tissueLayers{1}.modelParams);
            
            %Initialize structure:
            for counter = 1:length(fieldNames)
                theModelParams.(fieldNames{counter}) = {};                
            end
            
            
            for counter = 1:length(obj.tissueLayers)
                for counter2 = 1:length(fieldNames)
                    theModelParams.(fieldNames{counter2})(counter) = {obj.tissueLayers{counter}.modelParams.(fieldNames{counter2})};
                end
            end
            
        end %End get.modelParams
                
    end
    
    
    % ----------------- Tissue Methods ---------------------------%
    methods
        
        function obj = insertTissueObject(obj, TissueObject)
            
            %Validate the inputs:
            if ~isa(TissueObject, 'Tissue')
                ERRMSG = ['ERROR: Input object is not a Tissue Object'];
                error(ERRMSG);
            end
            
            if ~isequal(obj.Dims, TissueObject.Dims)
                ERRMSG = ['ERROR: Input Tissue object is not the same size as this Tissue object'];
                error(ERRMSG);
            end
            
            %-----Replace the s0 in the Tissue object:  -----%
            tissueS0 = obj.s0;
            insertedS0 = TissueObject.s0;
            lesion_location = find(TissueObject.noiseMask > 0);
            noise_threshold = obj.s0_threshold;
            
            %If the tissue is made up of complex data and the lesion is
            %not, assign phase information to it:
            if obj.complex_data && ~TissueObject.complex_data 
            
               %If at least a third of the voxels in the tissue have temp
               %magnitudes greater than the noise floor, assign the lesion
               %phase to match that phase, else assign it to 0:
               all_values = tissueS0(lesion_location);
               phase_locations = lesion_location(abs(all_values) > noise_threshold);
               the_values = tissueS0(phase_locations);
               
               if length(the_values) > length(all_values)*1/3    
                   mean_vector = mean(the_values);
                   mean_phase = angle(mean_vector);
                   apply_phase = true;
               else
                   apply_phase = false;
               end
               
               insertedS0 = complex(insertedS0);
               if apply_phase == true;
                   
                   %For voxels in the tissue with enough signal, assign the
                   %lesion to have that phase:
                   insertedS0(phase_locations) = complex(insertedS0(phase_locations).*cos(angle(tissueS0(phase_locations))), insertedS0(phase_locations).*sin(angle(tissueS0(phase_locations))));
                   
                   %for voxels without enough signal, assign them to have
                   %the mean phase:
                   other_locations = lesion_location(abs(all_values) <= noise_threshold);
                   insertedS0(other_locations) = complex(insertedS0(other_locations)*cos(mean_phase), insertedS0(other_locations).*sin(angle(mean_phase)));
               end
               
            end
               
            %Now, insert the insertedS0 into the tissueS0:
            tissueS0 = tissueS0 .*(abs(1-TissueObject.noiseMask)) + insertedS0;

            %Check for complex values:
            if obj.complex_data
                obj.s0_qtimage = qt_image(real(tissueS0));
                obj.s0_qtimageImag = qt_image(imag(tissueS0));
            else
                obj.s0_qtimage = qt_image(tissueS0);
            end
            
            
            %-----Now, insert the Layers:  -----%
            linked_layers = [];
            starting_layer = length(obj.tissueLayers)+1;
            for counter = 1:length(TissueObject.tissueLayers)
                
                %Skip Default nonenhancing Tissue:
                if strcmp(TissueObject.tissueLayers{counter}.layerName, 'Default nonenhancing Tissue')
                    continue;
                end
                
                %Update model parameters as appropriate:
                [obj, TissueObject.tissueLayers{counter}] = modelParamCheck(obj, TissueObject.tissueLayers{counter});
                
                %Update the tissue name: 
                TissueObject = TissueObject.setLayerName(counter, [TissueObject.tissueName, ' ', TissueObject.tissueLayers{counter}.layerName]);
                
                %Next, add the tissue layer to the stack of tissue layers:
                inserted_layer = length(obj.tissueLayers)+1;
                linked_layers = [linked_layers, inserted_layer];
                obj.tissueLayers(inserted_layer) = TissueObject.tissueLayers(counter);
                obj.insertTissue_pipeline(inserted_layer) = {TissueObject.insertTissue_pipeline{counter}-counter + inserted_layer};
                
                current_layer = logical(TissueObject.tissueLayers{counter}.mask.value);
                
                %Identify other existing layers this layer intersects with:
                for tissueCounter = 1:starting_layer - 1
                    
                    existing_layer = logical(obj.tissueLayers{tissueCounter}.mask.value);
                    the_intersect = existing_layer & current_layer;
                    intersect_loc = find(the_intersect(:) == true, 1, 'first');
                    
                    if ~isempty(intersect_loc)
                       obj.insertTissue_pipeline(tissueCounter) = { [obj.insertTissue_pipeline{tissueCounter}, inserted_layer]};
                    end
            
                end
                
            end
            
        end
        
        function obj = insertTissueLayer(obj, varargin)
        % ~for comments - if passing pieces of tissuelayer, first input must be the tissue mask.  
        % ~Inserting modifies masks of other layers
            
            %Determine if any pre-processing of the data to be added needs
            %to occur:
            [TissueLayerAdds] = obj.preProcessAddTissue(varargin{1:end});
            
            
            %Reversibly "insert" the tissue:
            for counter = 1:length(TissueLayerAdds)

                %Check the model parameters before inserting the layer:
                %First, add any fields from this tissue to the
                %obj.modelParams:
                the_fields = fieldnames(TissueLayerAdds{counter}.modelParams);
                for counter2 = 1:length(the_fields)
                    
                    %If its not an existing field and its not empty, add it
                    %to all layers:
                    if ~isfield(obj.modelParams, the_fields{counter2}) && ~isempty(TissueLayerAdds{counter}.modelParams.(the_fields{counter2}))
                        obj = obj.setEmptyModelParams(the_fields{counter2});
                    %If it is not an existing field and it is empty, remove
                    %it from the modelparams structure
                    %elseif ~isfield(obj.modelParams, the_fields{counter2}) && isempty(TissueLayerAdds{counter}.modelParams.(the_fields{counter2}))
                    %    TissueLayerAdds{counter}.modelParams = rmfield(TissueLayerAdds{counter}.modelParams, the_fields{counter2});
                    end
                end
                
                %Next, add any fields that exist in the current obj but not
                %in the tissue layer:
                obj_fields = fieldnames(obj.modelParams);
                for counter2 = 1:length(obj_fields)
                    if ~isfield(TissueLayerAdds{counter}.modelParams,obj_fields{counter2})
                        TissueLayerAdds{counter}.modelParams.(obj_fields{counter2}) = [];
                    end
                end
                
                %Next, add the tissue layer to the stack of tissue layers:
                inserted_layer = length(obj.tissueLayers)+1;
                obj.tissueLayers(inserted_layer) = TissueLayerAdds(counter);
                obj.insertTissue_pipeline(inserted_layer) = {[]};
                
                current_layer = logical(TissueLayerAdds{counter}.mask.value);
                
                %Identify other existing layers this layer intersects with:
                for tissueCounter = 1:inserted_layer - 1
                    
                    existing_layer = logical(obj.tissueLayers{tissueCounter}.mask.value);
                    the_intersect = existing_layer & current_layer;
                    intersect_loc = find(the_intersect(:) == true, 1, 'first');
                    
                    if ~isempty(intersect_loc)
                       obj.insertTissue_pipeline(tissueCounter) = { [obj.insertTissue_pipeline{tissueCounter}, inserted_layer]};
                    end
            
                end
            end    
            
        end
        
        function obj = removeTissueLayer(obj, number)
            
            %First validate the input:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
            if number == 1
                ERRMSG = ['ERROR: Layer 1 contains the non-enhancing tissue and can not be removed'];
                error(ERRMSG);
            end
            
            %Remove the tissue from the tissueLayers array:
            removeLayer = obj.tissueLayers{number};
            obj.tissueLayers(number) = [];
            
            %Update the insertTissue_pipline array as appropriate:
            %Remove the indicaited number from the insertTissue_pipline
            %list and decrease all numbers larger than it by 1.
            for counter = 1:length(obj.tissueLayers)
                layerTracker = obj.insertTissue_pipeline{counter};
                
                layerloc = find(layerTracker == number);
                if ~isempty(layerloc)
                    layerTracker(layerloc) = [];
                end
                
                layerTracker(layerTracker > number) = layerTracker(layerTracker > number) - 1;
                
                obj.insertTissue_pipeline(counter) = {layerTracker};
                
            end
            
            %remove the layer from tracking:
            inserted_layers = obj.insertTissue_pipeline{number};
            obj.insertTissue_pipeline(number) = [];
            
            %ensure that the area contained in the "removed tissue" is
            %either covered by anothr tissue or becomes part of the
            %unenhancing tissue:
            if isempty(find(obj.insertTissue_pipeline{1} == number))
                removemask = removeLayer.mask.value;
                for counter = number-1:-1:1
                    if isempty(find(obj.insertTissue_pipeline{counter} == number))
                        continue;
                    else 
                        removemask = removemask - obj.getLayerMask(counter);
                    end
                    
                end
                
                removemask(removemask < 0) = 0;
                mask1 = qt_image(obj.getLayerMask(1) + removemask);
                
                obj.tissueLayers{1}.mask = mask1;
            end

       
            %See if layers that were inserted need to be inserted into any
            %other layers:
            for counter = 1:length(inserted_layers)
                
                current_layer = logical(obj.tissueLayers{inserted_layers(counter)}.mask.value);
                
                for tissueCounter = 1:number - 1
                    
                    %First, check to see if the insrted layer is already
                    %recorded as being inserted into the current layer:
                    if ~isempty(find(obj.insertTissue_pipeline{tissueCounter} == inserted_layers(counter)))
                        continue;
                    end
                    
                    existing_layer = logical(obj.tissueLayers{tissueCounter}.mask.value);
                    the_intersect = existing_layer & current_layer;
                    intersect_loc = find(the_intersect(:) == true, 1, 'first');

                    if ~isempty(intersect_loc)
                       obj.insertTissue_pipeline(tissueCounter) = { [obj.insertTissue_pipeline{tissueCounter}, inserted_layers(counter)]};
                    end

                end
            end

        end
        
        function obj = flattenTissueLayers(obj)
           
            %"Flattens" all tissue layers, meaning that tissue layers can
            %no longer be reversibly removed
            
           for counter = 1:length(obj.tissueLayers)
                
                base_mask = obj.getLayerMask(counter);
               
                %Save the model parameters
                temp_params = fieldnames(obj.modelParams);
                for counter2 = 1:length(temp_params)
                    temp_params(counter2, 2) = obj.modelParams.(temp_params{counter2})(counter);
                end
                
                if strcmp(obj.tissueLayers{counter}.mask.elementClass, 'logical');
                    %Create a new TissueLayer object to avoid overwriting the users
                    %original input
                    obj.tissueLayers(counter) = {TissueLayer(qt_image(base_mask, 'elementClass','logical'), 'layerName', obj.layerNames{counter}, 't10', obj.t10(counter), 't2', obj.t2{counter})};
                    
                else
                    %Create a new TissueLayer object to avoid overwriting the users
                    %original input
                    obj.tissueLayers(counter) = {TissueLayer(qt_image(base_mask), 'layerName', obj.layerNames{counter}, 't10', obj.t10(counter), 't2', obj.t2{counter})};
                end
            
                %Set the layer parameters:
                for counter2 = 1:length(temp_params)
                    obj.tissueLayers{counter}.setModelParams(temp_params{counter2, 1}, temp_params{counter2, 2});
                end
                
                %Resest the insert tissue pipeline as all layers have
                %essentially been added to this layer now:
                obj.insertTissue_pipeline(counter) = {[]};
           end
        end
        
        function layerMask = getLayerMask(obj, number)
            
            %Validat the input:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
            %Insert other layers as needed:
            insert_layers = obj.insertTissue_pipeline{number};
            layerMask = obj.tissueLayers{number}.mask.value;
            
            mask_layers = zeros(size(layerMask)); 
            for layerCounter = insert_layers(1:end)
               
                 %Get the layer to insert:
                 mask_layers = mask_layers + obj.tissueLayers{layerCounter}.mask.value;
                 
            end
            
            %Remove any values over 1:
            mask_layers(mask_layers > 1) = 1;

            %invert the layer:
            mask_layers = abs(mask_layers - 1);
            
            %The value in mask_layers now defines the maximum value
            %remaining in any given voxel.
            layerMask(layerMask > mask_layers) = mask_layers(layerMask > mask_layers);            
            
        end
        
        function theEnhancementMap = getEnhancementMap(obj)
            
            theEnhancementMap = zeros([obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3)]);
                        
            %Loop by tissue layer:
            slice_weights = zeros([obj.xyzDims(2), obj.xyzDims(1),obj.xyzDims(3),length(obj.tissueLayers)]);
            for counter = 2:length(obj.tissueLayers)

                %Get the weights for each voxel in each layer
                temp_mask = obj.getLayerMask(counter);
                slice_weights(:,:,:, counter) = temp_mask;
            end

            %Account for the object mask:
            summation_mask = sum(slice_weights, 4);
            slice_weights(:,:,:,1) = obj.noiseMask-summation_mask;

            %Determine the max weighting.  Voxel will be assigned to
            %this color in the enhancement map:
            [max_val, the_index] = max(slice_weights,[], 4);
            the_index(obj.noiseMask == 0) = 0;

            theEnhancementMap = the_index;
                              
        end
        
        function obj = setLayerName(obj, number, name)
            
            %Validate the inputs:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
            
            if ~ischar(name)
                ERRMSG = ['Error using obj.setLayerName: Expected ''name'' to be a string, instead it is a %s'];
                error(ERRMSG, class(name));
            end
            
            obj.tissueLayers{number}.layerName = name;
            
        end %end Tissue.setLayername
        
        function obj = setT10(obj, number, value)
            
            %Validate the inputs:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
                        
            if ~isempty(value)
                if ischar(value)
                    valid_strings = obj.tissueLayers{number}.defaultTissueTypes;
                    validatestring(value, valid_strings);
                elseif isequal(value, -1)
                    value = [];
                else
                    validateattributes(value, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan', 'positive'});
                end
            end

            
            obj.tissueLayers{number}.t10 = value;
            
        end %obj.setT10
        
        function obj = setT2(obj, number, value)
            
            %Validate the inputs:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
                        
            if ~isempty(value)
                if ischar(value)
                    valid_strings = obj.tissueLayers{number}.defaultTissueTypes;
                    validatestring(value, valid_strings);
                elseif isequal(value, -1)
                    value = [];
                else
                    validateattributes(value, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan', 'positive'});
                end
            end

            obj.tissueLayers{number}.t2 = value;
            
        end %obj.setT2
        
        function obj = setModelParams(obj, number, varargin) %theFieldname, value)
            
            %Validate the inputs:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.layerNames), '>', 1 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
            
            %Get a list of current fieldnames:
            current_fieldnames = fieldnames(obj.tissueLayers{1}.modelParams);
            
            %Pull out the properties:
            for fieldcounter = 1:2:length(varargin)
            
                theFieldname = varargin{fieldcounter};
                value = varargin{fieldcounter+1};
                
                if ~isvarname(theFieldname)
                    ERRMSG = ['ERROR: The selected fieldname is not valid'];
                    error(ERRMSG);
                end
            
                %If a valid fieldname has been supplied, go ahead and add the
                %the property.  validation of the value will have to be
                %performed by the simulation model function because we don't
                %know ahead of time what the model is and by extension what
                %valid values are.

                %If there is an existing, matching fieldname, go ahead and set
                %the value:
                exist_field = false;
                for counter = 1:length(current_fieldnames)
                   if strcmp(theFieldname, current_fieldnames{counter})
                       exist_field = true;
                       break;
                   end
                end

                %If there is not an existing field, add it to all layers.  Set
                %defult value to empty;
                if ~exist_field && ~isempty(value)
                    obj = obj.setEmptyModelParams(theFieldname);
                end

                %If the value is being set to empty, remove that particular
                %fieldname if the value is empty for all tissue layers:
                if exist_field && isempty(value)
                    remove_field = true;

                    for counter = 1:length(obj.tissueLayers)

                        if ~isempty(obj.tissueLayers{counter}.modelParams.(theFieldname))
                            remove_field = false;
                             break;
                        end
                    end

                    if remove_field
                        for counter = 1:length(obj.tissueLayers)
                            obj.tissueLayers{counter}.removeModelParams(theFieldname);
                        end
                        return;
                    end
                end

                %Set the model parameter for the specified layer:
                obj.tissueLayers{number}.setModelParams(theFieldname, value);
            end

            
        end %obj.setModelParams  
                
        function obj = combineLayers(obj, varargin)
        %combineLayers combines layers with matching tissue properties.
        %with no inputs, it will search for all layers that match and
        %combine them.  Alternatly, provide a vector of layers to combine 
        %or a list of properties to match. Only layers with matching values 
        %will be combined.  
        
        %   input options:
        %       *vector contining layer numbers
        %       *cell array containing model parameter names
        
        
            if isempty(varargin)

                untested_layers = [1:length(obj.tissueLayers)];
                theParams = fieldnames(obj.modelParams);
                
            else
                
                %sort out the input 
                counter = 1;
                untested_layers = [];
                theParams = {};
                while counter <= length(varargin)
                    
                    %Check for vector containing layer numbers to combine:
                    if isnumeric(varargin{counter})
                        
                        if ~isempty(untested_layers)
                            ERRMSG = ['ERROR: Additional numeric input provided to combineLayers.  Only one \n'...
                                        'numeric input allowed to set the layers to be combined. %s'];
                            error(ERRMSG, '');
                        end
                        
                        untested_layers = varargin{1};
                
                        validateattributes(untested_layers, {'numeric'}, {'vector', 'positive', '<=', length(obj.layerNames), 'nonempty'});
                        if ~isequal(ceil(untested_layers), floor(untested_layers))
                            ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                            error(ERRMSG);
                        end
                    elseif iscell(varargin{counter})
                        
                        if ~isempty(theParams)
                            ERRMSG = ['ERROR: Additional cell array input provided to combineLayers.  Only one \n'...
                                        'cell array input containing the names of relevent modelParams allowed. %s'];
                            error(ERRMSG, '');
                        end
                        
                        %validate the inputs in the cell array
                         theParams = varargin{counter};
                         
                         validparams = fieldnames(obj.modelParams);
                         
                         for counter2 = 1:length(theParams)
                            validatestring(theParams{counter2}, validparams, counter2);
                         end
                    else
                        ERRMSG = ['Unknown input of class %s.  Valid inputs include a \n' ...
                                  'cell array containing a list of model parameters or \n' ...
                                  'a numeric vector containing layers to combine'];
                        error(ERRMSG, class(varargin{counter}));
                    end
                    
                    counter = counter + 1;
                end
                 
                if isempty(untested_layers) 
                    untested_layers = [1:length(obj.tissueLayers)];
                end
                if isempty(theParams)
                    theParams = fieldnames(obj.modelParams);
                end     
            end
            
            unmatched_layers = [];
            matched_layers = {};

            %Find all matching tissue layers:
            while ~isempty(untested_layers)

                %compare layers against each other to see if they
                %match: 
                num = untested_layers(1);
                the_match = [];

                %compare the first untested layer with all
                %subsequent untested layers:
                remove_layers = [];
                %theParams = fieldnames(obj.modelParams);
                for counter = 2:length(untested_layers)
                    layerMatch = true;

                    %Compare all tissue and model properties:
                    if ~isequal(obj.t10(num), obj.t10(untested_layers(counter)))
                        continue;
                    elseif ~isequal(obj.t2(num), obj.t2(untested_layers(counter)))
                        continue;
                    end

                    for paramCounter = 1:length(theParams)
                        if ~isequal(obj.modelParams.(theParams{paramCounter}){untested_layers(1)}, obj.modelParams.(theParams{paramCounter}){untested_layers(counter)})
                            layerMatch = false;
                            break;
                        end
                    end

                    if ~layerMatch
                        continue;
                    end

                    the_match = [the_match, untested_layers(counter)];
                    remove_layers = [remove_layers, counter];
                end

                if ~isempty(the_match)
                    matched_layers(length(matched_layers) + 1) = {[num, the_match]};

                    %Remove matched layers from the list of untested
                    %layers:
                    untested_layers(remove_layers) = [];
                else
                    unmatched_layers(length(unmatched_layers) + 1) = untested_layers(1);
                end

                untested_layers = untested_layers(2:end);
            end

            %If no layers match, return to calling function
            if isempty(matched_layers)
                return;
            end
            
            %Flatten all layers to simplify the combining process:
            obj = obj.flattenTissueLayers;
            
            %Combine matching layers:
            for counter = 1:length(matched_layers)

                base_layerNum = matched_layers{counter}(1);
                base_mask = obj.getLayerMask(base_layerNum);
                newlayerName = obj.tissueLayers{base_layerNum}.layerName;

                comb_logical = strcmp(obj.tissueLayers{base_layerNum}.mask.elementClass, 'logical');
                for counter2 = 2:length(matched_layers{counter})
                    comb_layerNum = matched_layers{counter}(counter2);
                    
                    if base_layerNum > 1
                        base_mask = base_mask + obj.getLayerMask(comb_layerNum);
                        newlayerName = [newlayerName ' + ' obj.tissueLayers{comb_layerNum}.layerName];
                    end

                    if ~strcmp(obj.tissueLayers{comb_layerNum}.mask.elementClass, 'logical')
                        comb_logical =  false;
                    end
                end

                temp_params = fieldnames(obj.modelParams);
                for counter2 = 1:length(temp_params)
                    temp_params(counter2, 2) = obj.modelParams.(temp_params{counter2})(base_layerNum);
                end
                
                
                if comb_logical
                    %obj.tissueLayers{base_layerNum}.mask = qt_image(base_mask, 'elementClass','logical' );

                    %Create the 'combined' tissue layer by creating a new
                    %TissueLayer object to avoid overwriting the users
                    %original input
                    obj.tissueLayers(base_layerNum) = {TissueLayer(qt_image(base_mask, 'elementClass','logical'), 'layerName', newlayerName, 't10', obj.t10(base_layerNum), 't2', obj.t2{base_layerNum})};
                    
                else
                    %obj.tissueLayers{base_layerNum}.mask = qt_image(base_mask);
                    
                    %Create the 'combined' tissue layer by creating a new
                    %TissueLayer object to avoid overwriting the users
                    %original input
                    obj.tissueLayers(base_layerNum) = {TissueLayer(qt_image(base_mask), 'layerName', newlayerName, 't10', obj.t10(base_layerNum), 't2', obj.t2{base_layerNum})};
                end
            
                %Set the layer parameters:
                for counter2 = 1:length(temp_params)
                    obj.tissueLayers{base_layerNum}.setModelParams(temp_params{counter2, 1}, temp_params{counter2, 2});
                end
                
                %Resest the insert tissue pipeline as all layers have
                %essentially been added to this layer now:
                obj.insertTissue_pipeline(base_layerNum) = {[]};
            end
            
            %Remove the combined layers:
            for counter = 1:length(matched_layers)

                for counter2 = 2:length(matched_layers{counter})
                    comb_layerNum = matched_layers{counter}(counter2);

                    %remove the layer:
                    obj = obj.removeTissueLayer(comb_layerNum);
                    
                    %decrement all values in matched layers less than
                    %comb_layerNum to accounter for removing this layer
                    for counter3 = 1:length(matched_layers)
                        for counter4 = 1:length(matched_layers{counter3})
                            if matched_layers{counter3}(counter4) > comb_layerNum
                                matched_layers{counter3}(counter4) = matched_layers{counter3}(counter4) - 1;
                            end
                        end
                    end

                end
            end

        end
        
        function obj = generateBPE(obj, baseMask, numLayers, varargin)
                
            %Required inputs:
            % baseMask - can be qt_image object or an array with values
            %            between 0 and 1.  Defines the fibroglanular tissue
            %            extent.
            % numLayers - defines how many different enhancement profiles
            %             are desired for the BPE.  Value greater than 1.
            %Optional inputs:
            % BPE - set the level of BPE.  Options: 
            %      *'minimal', 'mild', 'moderate', 'marked' -> will
            %       default to 10%, 37%, 63% and 87% tissue enhancement
            %       respectivly
            %      * Can also enter a decimal value between 0 and 1 to
            %       explicitly set the percentage of fibroglandular tissue 
            %       enhancing.
            % seedPercent - defines the percentage of voxels in the
            %       fibroglandular tissue used as seed locations for 
            %       defining the areas of enhancement. Default value is
            %       .0005 or 0.05%.  Decreasing this value will give more
            %       isolated pockets of ehnancement, increasing will
            %       provide more continuous regions of enhancement.
            % center - center point between the breasts.  Default is center
            %       of imaging volume.
            % x_radius / y_radius / z-radius - distance from the
            %       centerpoint to the lateral edge, anterior/posterior edge and
            %       superior/inferior edges
            % rand_seed - seed for the random number generator
            % 
            
            if isempty(baseMask) || isempty(numLayers)
                ERRMSG ='Error in calling generateBPE: inputs baseMask and numLayers are required';
                error(ERRMSG);
            end
            

            %Validate the required inputs:
            if isa(baseMask, 'TissueLayer')
                the_mask = baseMask.mask.value;
            elseif isa(baseMask, 'qt_image')
                the_mask = baseMask.value;
            elseif ~ismatrix(baseMask)
                the_mask = qt_image(baseMask).value;
            end
            validateattributes(the_mask, {'logical', 'double'}, {'size', [obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3)], 'nonnegative', '<=', 1});
            clear baseMask;
            
            if ~isequal(floor(numLayers), ceil(numLayers))
                ERRMSG = ['Error in selecting the number of BPE layers: expected the number of layers to be a nondecimal number.  numLayers = %d'];
                error(ERRMSG, numLayers);
            end
            
            validateattributes(numLayers, {'numeric'}, {'scalar', 'positive'});
            
            %Set up defaults before sorting through the varargin variables:
            BPE_percent = .37; %Default setting for mild BPE
            seed_locs_percent = .0005;  %default number of seed locs to use: .05% of all valid indicies
            mid_point = round(size(the_mask)/2);
            x_rad = (2/3)*mid_point(2);
            y_rad = (2/3)*mid_point(1);
            z_rad = (2/3)*mid_point(3);
            layerName = 'BPE';
            nbins = 10;
            rand_seed = 0;
            
            counter = 1;
            while counter <= length(varargin)
                if isstr(varargin{counter})
                    if strcmp(varargin{counter}, 'BPE')
                        
                        BPE_value = varargin{counter + 1};
                        
                        if isstr(BPE_value)
                            BPE_value = lower(BPE_value);
                            valid_strings = {'minimal', 'mild', 'moderate', 'marked'};
                            validatestring(BPE_value, valid_strings);

                            if strcmp(BPE_value, 'minimal')
                                BPE_percent = 0.01;
                            elseif strcmp(BPE_value, 'mild')
                                BPE_percent = 0.37*.5;
                            elseif strcmp(BPE_value, 'moderate')
                                BPE_percent = 0.63*.5;
                            elseif strcmp(BPE_value, 'marked')
                                BPE_percent = 0.87*.5;
                            end
                        else
                            validateattributes(BPE_value, {'double'}, {'scalar', 'nonnegative', '<=' , 1});
                            BPE_percent = BPE_value;
                        end
                        clear BPE_value;
                    elseif strcmp(varargin{counter}, 'seedPercent')
                        validateattributes(varargin{counter + 1}, {'double'}, {'scalar', 'positive', '<=', 1});
                        
                        seed_locs_percent = varargin{counter + 1};
                        
                        if seed_locs_percent > 0.05
                            WRNMSG('Warning: selcting more the 5% of all indicies to be seed locs may slow performance')
                            warning(WRNMSG);
                        end
                        
                    elseif strcmp(varargin{counter}, 'center')
                        validateattributes(varargin{counter + 1}, {'numeric'}, {'vector', 'numel', 3, 'positive'});
                        
                        mid_point = varargin{counter};
                        
                        if ~isequal(ceil(mid_point), floor(mid_point))
                            ERRMSG('Error in setting matrix center point: Expected non decimal inputs');
                            error(ERRMSG);
                        elseif mid_point(1) > obj.xyzDims(2) || mid_point(2) > obj.xyzDims(1) || midpoint(3) > obj.xyzDims(3)
                            ERRMSG('Error in setting matrix center point: Expected center point to be less than size(baseMask).  Requested center is [%d, %d, %d], but matrix size is: [%d, %d, %d].');
                            error(ERRMSG, mid_point(1), mid_point(2), mid_point(3), obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3));
                        end
                        
                    elseif strcmp(varargin{counter}, 'xradius')
                        validateattributes(varargin{counter+1}, {'double'}, {'scalar', 'positive', '<=', obj.xyzDims(1)});
                        x_rad = varargin{counter + 1};
                    elseif strcmp(varargin{counter}, 'yradius')
                        validateattributes(varargin{counter+1}, {'double'}, {'scalar', 'positive', '<=', obj.xyzDims(2)});
                        y_rad = varargin{counter + 1};
                    elseif strcmp(varargin{counter}, 'zradius')
                        validateattributes(varargin{counter+1}, {'double'}, {'scalar', 'positive', '<=', obj.xyzDims(3)});
                        z_rad = varargin{counter + 1};
                    elseif strcmp(varargin{counter}, 'layerName')
                        if isstr(varargin{counter+1})
                            layerName = varargin{counter + 1};
                        end
                    elseif strcmp(varargin{counter+1}, 'nbins')
                        validateattributes(varargin{counter+1}, {'numeric'}, {'positive', 'scalar'});
                        nbins = varargin{counter+1};
                    elseif strcmp(varargin{counter+1}, 'rand_seed')
                        rand_seed = varargin{counter+1};
                    end
                    counter = counter + 1;
                end
                counter = counter + 1;
            end
               
            %Ok, now actually start performing calculations for BPE:
            volume = zeros(size(the_mask));
            volume(the_mask > 0) = 1;
            
            the_indicies = find(volume == 1);
            num_pts = length(the_indicies);
            
            %Check for values of seed_locs_percent that are low enough we
            %will esseintially have no seeds.
            if round(num_pts * seed_locs_percent) < 20 && num_pts > 40
                seed_locs_percent = 20/num_pts;
            end
                    
            %divide breast into regions:
            xlimits = 0:mid_point(2)/nbins:mid_point(2);
            ylimits = 0:mid_point(1)/nbins:mid_point(1);
            zlimits = 0:mid_point(3)/nbins:mid_point(3);

            mean_bin = find(xlimits > x_rad, 1, 'first'); %bin containing the indicaited radii
            region_struct = cell([nbins, 1]);

            [yval, xval, zval] = ind2sub(size(volume), the_indicies);
            remaining_indicies = the_indicies;
            for counter = 2:length(xlimits)

                ovoid =  ((xval-mid_point(2)).^2)/(xlimits(counter).^2) + ((yval-mid_point(1)).^2)/(ylimits(counter).^2) + ((zval-mid_point(3)).^2)/(zlimits(counter).^2);
                region_locs = find(ovoid <=1);
                region_struct(counter - 1) = {[remaining_indicies(region_locs)]};  

                xval(region_locs) = [];
                yval(region_locs) = [];
                zval(region_locs) = [];
                remaining_indicies(region_locs) = [];
                clear region_locs; 
            end

            %Randomly set the seed locations:
            rng(rand_seed);
            
            the_stdev = nbins/3;
            %cumulative distribution function to help determin seeds per
            %bin:
            cumulative_dist = 0.5*[1+erf(((0:nbins) - mean_bin)/(the_stdev*sqrt(2)))]; 
 
            %Calculate fraction of seeds per bin
            bin_fraction = diff(cumulative_dist/cumulative_dist(end));
            %There is no '0' bin, include those seeds in the '1' bin
            bin_fraction(1) = bin_fraction(1) + cumulative_dist(1); 
 
            num_seeds_bin = ones([length(region_struct),1]);
            seed_locs = zeros([round(length(the_indicies)*seed_locs_percent*5), 1]); 
            for counter = 1:nbins
                num_seeds_bin(counter) = round(seed_locs_percent * (bin_fraction(counter)/(1/nbins)) * length(region_struct{counter}));
                seed_locs(sum(num_seeds_bin(1:counter-1))+1:sum(num_seeds_bin(1:counter-1))+num_seeds_bin(counter)) = region_struct{counter}(randperm(length(region_struct{counter}), num_seeds_bin(counter)));
            end

            seed_locs(seed_locs == 0) = [];
 
            %Grow the seeds:
             x = -10:10;
             y = -10:10;
             z = -10:10;
             [X, Y, Z] = meshgrid(x, y, z);
             seed_values = (sqrt(X.^2 + Y.^2 + Z.^2)).^-1;
             seed_values(seed_values == Inf) = 1;
             
             [row, col, slice] = ind2sub(size(volume), seed_locs);
             value_mask = zeros(size(volume));
             for counter = 1:length(row)

                 if row(counter)-10 < 1 || row(counter) + 10 > size(volume, 1)
                     continue;
                 elseif col(counter)-10 < 1 || col(counter) + 10 > size(volume, 2)
                     continue;
                 elseif slice(counter)-10 < 1 || slice(counter) + 10 > size(volume, 3)
                     continue;
                 end
                 
                 %Here we are essentially adding a tile at each seed
                 %location:
                 value_mask(row(counter) - 10:row(counter) + 10, col(counter)-10: col(counter)+10, slice(counter)-10: slice(counter) + 10) = ...
                     value_mask(row(counter) - 10:row(counter) + 10, col(counter)-10: col(counter)+10, slice(counter)-10: slice(counter) + 10) + seed_values;
             end

            %Limit the values to only those that are valid locations in the
            %mask:
             value_mask = value_mask .* volume;

             % use this map to generate the desired masks:
             [value,ind] = sort(value_mask(:), 'descend');
             num_mask_locs = round(BPE_percent * length(the_indicies));
             
             locs_per_layer = zeros([numLayers, 1]);
             for counter = 1:numLayers
                 locs_per_layer(counter) = round(num_mask_locs * ((numLayers - counter + 1)/numLayers)^3);
             end
             
             for counter = 1:numLayers
                layerMask = zeros(size(volume));
                layerMask(ind(1:(locs_per_layer(counter)))) = 1;
                
                if islogical(the_mask)
                    layerMask = logical(layerMask);
                else
                    [layerMask] = obj.blendMask(layerMask);
                    layerMask = abs(layerMask).*the_mask;
                    layerMask(layerMask > 1) = 1;
                    layerMask(layerMask < 0) = 0;
                end
                
                obj = obj.insertTissueLayer(layerMask, 'layerName', [layerName num2str(counter)]);
             end

        end
        
        function obj = generateHeterogeneous(obj, baseMask, numLayers, varargin)
           
            if isempty(baseMask) || isempty(numLayers)
                ERRMSG ='Error in calling generateBPE: inputs baseMask and numLayers are required';
                error(ERRMSG);
            end
            
            %Validate the required inputs:
            if isa(baseMask, 'TissueLayer')
                the_mask = baseMask.mask.value;
            elseif isa(baseMask, 'qt_image')
                the_mask = baseMask.value;
            elseif ~ismatrix(baseMask)
                the_mask = qt_image(baseMask).value;
            end
            validateattributes(the_mask, {'logical', 'double'}, {'size', [obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3)], 'nonnegative', '<=', 1});
            clear baseMask;
            
            validateattributes(numLayers, {'numeric'}, {'scalar', 'positive'});
            
            %Set up defaults for allowed varaibles in varargin:
            layerName = 'Heterogeneous';
            percent_seed_indicies = 0.01;
            
            %Randomly assign locations to each layer:
            rng(0);
            
            %identify the locations in the mask:
            the_indicies = find(the_mask > 0);
            
            if length(the_indicies) < numLayers
                ERRMSG('ERROR: not enough locations in mask to assign to different layers');
                error(ERRMSG);
            end

            %assign seeds for determining the layer:
            num_seed_indicies = round(percent_seed_indicies*length(the_indicies));
            seeds_per_layer = floor(num_seed_indicies / numLayers);
            
            %check for too few seed locs:
            if seeds_per_layer < 5 && length(the_indicies)>=5*numLayers
                seeds_per_layer = 5;
            elseif seeds_per_layer < 5 && length(the_indicies) < 5*numLayers
                seeds_per_layer = floor(seed_per_layer / numLayers);
            end
            
            %randomly assign the seed locs:
            seed_indicies = randperm(length(the_indicies), num_seed_indicies);
            
            %Sort locs by layer:
            seed_locs = zeros(numLayers, seeds_per_layer); 
            seed_rcs = zeros(numLayers, seeds_per_layer, 3);
            for counter = 1:numLayers
                seed_locs(counter, :) = the_indicies(seed_indicies((seeds_per_layer * (counter-1)+1):(seeds_per_layer * counter)));
                [row, col, slice] = ind2sub(size(the_mask), seed_locs(counter, :));
                seed_rcs(counter, :, 1) = row;
                seed_rcs(counter, :, 2) = col; 
                seed_rcs(counter, :, 3) = slice;
            end
            
            %gaussian tile:
            sigmax = 2 * obj.xyzRes(1)/obj.xyzRes(1);
            sigmay = 2 * obj.xyzRes(1)/obj.xyzRes(2);
            sigmaz = 2 * obj.xyzRes(1)/obj.xyzRes(3);
            
            x = -10:10;
            y = -10:10; 
            z = -10:10;
 
            [X, Y, Z] = meshgrid(x, y, z); 
            
            the_tile = exp(-(X.^2/(2*sigmax^2) + Y.^2/(2*sigmay^2) + Z.^2/(2*sigmaz^2)));
            the_tile = the_tile * 1/sum(the_tile(:));
            
            %Build the new layers:
            the_layers(1:numLayers) = {zeros(size(the_mask))};
            
            for counter = 1:numLayers
                for counter2 = 1:seeds_per_layer
                    row = seed_rcs(counter, counter2, 1);
                    col = seed_rcs(counter, counter2, 2);
                    slice = seed_rcs(counter, counter2, 3);
                    
                    if row < 11 || row > obj.xyzDims(2)-10
                        continue;
                    elseif col < 11 || col > obj.xyzDims(1)-10
                        continue;
                    elseif slice < 11 || slice > obj.xyzDims(3)-10
                        continue;
                    end
                    
                    the_layers{counter}(row-10:row+10, col-10:col+10, slice-10:slice+10) = the_layers{counter}(row-10:row+10, col-10:col+10, slice-10:slice+10) + the_tile;
                end
            end
           
            %Assign pixels based on which layer has the max value for that voxel:
            for counter = 1:numLayers
                for counter2 = counter+1:numLayers
                    test = the_layers{counter} - the_layers{counter2};
                    
                    %If another layer has a larger value, set voxel to 0:
                    the_layers{counter}(test < 0) = 0;
                    %If this layer is larger, set other layer voxel to 0:
                    the_layers{counter2}(test > 0) = 0;
                    %If they are the same and not zero, set other layer
                    %voxel value to 0:
                    the_layers{counter2}(test == 0 & the_mask > 0) = 0;
                end
                
                the_layers{counter}(the_layers{counter} > 0) = 1;
                the_layers{counter} = the_layers{counter} .* ceil(the_mask);
            end
                
            %Create tissue layer obejects and add them to the Tissue:
            for counter = 1:numLayers
%                 if ~islogical(the_mask)
%                     [the_layers{counter}] = obj.blendMask(the_layers{counter});
%                 end
                obj = obj.insertTissueLayer(the_layers{counter});
                obj.setLayerName(length(obj.tissueLayers), [layerName ' ' num2str(counter)]);
            end
            
        end
        
        function [filtered_mask] = blendMask(obj, the_mask, varargin)
           
            if ~isempty(varargin)
               
                counter = 1;
                while counter <=length(varargin)
                    if ischar(varargin{counter})
                        if strcmp('sigmax', varargin{counter})
                            sigmax = varargin{counter + 1};
                        elseif strcmp('sigmay', varargin{counter})
                            sigmay = varargin{counter + 1};
                        elseif strcmp('sigmaz', varargin{counter})
                            sigmaz = varargin{counter + 1};
                        else
                            ERRMSG = 'Error in blending mask, unrecognized input';
                            error(ERRMSG);
                        end
                        counter = counter + 1;
                    end
                    counter = counter + 1;
                end
            end
            
            if ~exist('sigmax', 'var')
                %sigmax = 1 / round(obj.xyzDims(1)) * obj.xyzRes(1)/obj.xyzRes(1);
                sigmax = 0.5 * obj.xyzRes(1)/obj.xyzRes(1);
            end
            if ~exist('sigmay', 'var')
                %sigmay = 1 / round(obj.xyzDims(2))* obj.xyzRes(1)/obj.xyzRes(2);
                sigmay = 0.5 * obj.xyzRes(1)/obj.xyzRes(2);
            end
            if ~exist('sigmaz', 'var')
                %sigmaz = 1 / round(obj.xyzDims(3))* obj.xyzRes(1)/obj.xyzRes(3);
                sigmaz = 0.5 *  obj.xyzRes(1)/obj.xyzRes(3);
            end
            
            %x = -floor(obj.xyzDims(1)/2):floor(obj.xyzDims(1)/2-1);
            %y = -floor(obj.xyzDims(2)/2):floor(obj.xyzDims(2)/2-1);
            %z = -floor(obj.xyzDims(3)/2):floor(obj.xyzDims(3)/2-1);
            
            x = -1:1;
            y = -1:1; 
            z = -1:1;
 
            [X, Y, Z] = meshgrid(x, y, z);
            
            %the_gaussian = exp(-(X.^2*(2*sigmax^2) + Y.^2*(2*sigmay^2) + Z.^2*(2*sigmaz^2)));
            
            the_gaussian = exp(-(X.^2/(2*sigmax^2) + Y.^2/(2*sigmay^2) + Z.^2/(2*sigmaz^2)));
            the_gaussian = the_gaussian * 1/sum(the_gaussian(:));
            
            filtered_mask = imfilter(the_mask, the_gaussian);
            
            %fftmask = fftshift(fftshift(fftshift(fft(fft(fft(fftshift(fftshift(fftshift(the_mask, 1), 2), 3), [], 1), [], 2), [], 3), 1), 2), 3).*the_gaussian;
            
            %filtered_mask = abs(ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(fftmask, 1), 2), 3), [], 1), [], 2), [], 3), 1), 2), 3));
        end
        
        function obj = setPrecession(obj, relppm, relAmp)
            %Set up the precession model.  ppm is relative to water
            
            %Validate the inputs:
            validateattributes(relppm, {'numeric'}, {'vector', '<=', 7, '>=', -7, 'nonempty'});
            validateattributes(relAmp, {'numeric'}, {'vector', '<=', 1, 'nonnegative', 'nonempty', 'numel', length(relppm)});
            
            if ~isequal(sum(relAmp(:)), 1)
                ERRMSG = ['ERROR: relative amplitude of all spectral components must sum to 1'];
                error(ERRMSG);
            end
            
            obj.signalPrecession.ppm = relppm;
            obj.signalPrecession.relAmp = relAmp;
                        
        end
        
        function obj = setUserParams(obj, fieldname, value)
            
            %Validate the fieldname:
            if ~isvarname(fieldname)
                ERR_MSG = ['Error in setUserParams: the supplied field name is not a valid matlab variable name'];
                error(ERR_MSG);
            end
            
            obj.userParams.(fieldname) = value;

        end
            
            
    end %End Tissue Methods
    
    methods (Access = 'protected', Hidden = true) 
        
        function [TissueLayerAdds] = preProcessAddTissue(obj,varargin)
                
            
            %Determine if we have been passed TissueLayer objects or
            %components to build a tissueLayer object:
            TissueLayerAdds = {};
            other_inputs = {};
            for counter = 1:length(varargin)
                if strcmp(class(varargin{counter}), 'TissueLayer')
                    TissueLayerAdds(length(TissueLayerAdds)+1) = varargin(counter);
                else
                    other_inputs(length(other_inputs)+ 1) = varargin(counter);
                end
            end
            
            %construct additional TissueLayer objects if needed:
            if ~isempty(other_inputs)
                %The first input should still be the information to make
                %the mask object:
                %TO DO: Add try/catch statements arount setting of qt_image
                %to explain error to the user

                if islogical(other_inputs{1})
                    the_mask = qt_image(double(other_inputs{1}));
                else
                    the_mask = qt_image(other_inputs{1});
                end

                TissueLayerAdds(length(TissueLayerAdds) +1) ={TissueLayer(the_mask, other_inputs{2:end})};
            end
            
            for counter = 1:length(TissueLayerAdds)
                                            
                %Verify that the mask is the proper size:
                if ~isequal(TissueLayerAdds{counter}.mask.dimSize, [obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3)]);
                    ERRMSG = ['Error in adding new mask layer: the provided mask does not match ...\n'
                        'the tissue size of [ %d, %d, %d]'];
                    error(ERRMSG, obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3));
                end
                
                if isempty(TissueLayerAdds{counter}.t10)
                    TissueLayerAdds{counter}.t10 = obj.defaultT10;
                end
                
            end
            
        end
                
        function [TissueLayerAdds] = OldpreProcessAddTissue(obj, varargin)
        %This function sorts out the inputs used to add new tissues   
        
        
            %If there is just one input, the choices are limited:
            if length(varargin) < 2
                %If this is a tissue layer object, set it to the temporary
                %variable:
                if isa(varargin{1}, 'TissueLayer')
                    TissueLayerAdds = varargin(1);

                %If this is not a tissueLayer object, it must be the input
                %required to make one.  Use it to make a tissueLayer
                %object:
                else
                    TissueLayerAdds = {TissueLayer(qt_image(varargin{1}))};

                end
            else
            %If there is more than 1 input, we need to identify the inputs
            %because we may need to call some additional functions:
            
                %The first input should still be the information to make
                %the mask object:
                %TO DO: Add try/catch statements arount setting of qt_image
                %to explain error to the user
                the_mask = qt_image(varargin{1});
                
                %Now look for ktrans, ve and vp.  If more than 1 value is
                %sent in, call function for heterogeneous assignement of
                %enhancement values:
                ktrans = [];
                ve = []; 
                vp = [];
                t10 = [];
                other_input = {};
                counter = 2;
                while counter <= length(varargin)-1
                    if ischar(varargin{counter})
                        if strcmp(varargin{counter}, 'ktrans')
                            ktrans = varargin{counter + 1};
                            counter = counter + 1;
                        elseif strcmp(varargin{counter}, 've')
                            ve = varargin{counter + 1};
                            counter = counter + 1;
                        elseif strcmp(varargin{counter}, 'vp')
                            vp = varargin{counter + 1};
                            counter = counter + 1;
                        elseif strcmp(varargin{counter}, 't10');
                            t10 = varargin{counter + 1};
                            counter = counter + 1;
                        else
                            other_input(length(other_input)+ 1) = varargin(counter);
                        end
                    else
                        other_input(length(other_input)+ 1) = varargin(counter);
                    end
                    counter = counter + 1;
                end

                %If they are all either scalar, empty or a string 
                %specifying a default tissue type, this is just a
                %simple add: 
                if (isempty(ktrans) || isscalar(ktrans) || ischar(ktrans)) && ...
                   (isempty(ve)     || isscalar(ve)     || ischar(ve))     && ...
                   (isempty(vp)     || isscalar(vp)     || ischar(vp))     && ...
                   (isempty(t10)    || isscalar(t10)    || ischar(t10)) 
                   
                    TissueLayerAdds = {TissueLayer(the_mask, varargin{2:end})};
                    
                else
                
                    %If this is not a simple add, verify inputs are appropriate
                    %for heterogeneous add:

                    %determine the size of each variable:
                    prop_size = cell([1,4]);
                    prop_size(1,1) = {size(ktrans)};
                    prop_size(1,2) = {size(ve)};
                    prop_size(1,3) = {size(vp)};
                    prop_size(1,4) = {size(t10)};

                    alternate_size = [1,1];
                    for counter = 1:length(prop_size)
                       temp_compare = prop_size{1,counter};

                       %Rule out a few cases:
                       %Can not have more than 3 dimensions:
                       if length(temp_compare) > 3
                            ERRMSG = ['Error attempting to add new tissue: when assigning T10, ktrans, ve and vp values, \n'...
                                      'the number of dimensions of the variables can not exceed 3.  \n'...
                                      'The size of the inputs are: \n'...
                                      'size(ktrans)=[%d, %d, %d] \n'... 
                                      'size(ve)    =[%d, %d, %d] \n'...
                                      'size(vp)    =[%d, %d, %d] \n'...
                                      'size(t10)   =[%d, %d, %d]']; 


                            error(ERRMSG, size(ktrans, 1), size(ktrans, 2), size(ktrans, 3),...
                                          size(ve, 1),     size(ve, 2),     size(ve, 3),...
                                          size(vp, 1),     size(vp, 2),     size(vp, 3),...
                                          size(t10, 1),    size(t10, 2),    size(t10, 3));
                       %If 3 dimensions, must be a map the same size as the
                       %images:
                       elseif length(temp_compare) == 3 && ~isequal(temp_compare, obj.xyzDims)
                           ERRMSG = ['Error attempting to add new tissue: when assigning T10, ktrans, ve and vp values, \n'...
                                      'if the number of dimensions are 3, the size must match the size of the tissue images.  \n'...
                                      'The size of the inputs are: \n'...
                                      'size(ktrans)=[%d, %d, %d] \n'... 
                                      'size(ve)    =[%d, %d, %d] \n'...
                                      'size(vp)    =[%d, %d, %d] \n'...
                                      'size(t10)   =[%d, %d, %d]']; 

                           error(ERRMSG, size(ktrans, 1), size(ktrans, 2), size(ktrans, 3),...
                                          size(ve, 1),     size(ve, 2),     size(ve, 3),...
                                          size(vp, 1),     size(vp, 2),     size(vp, 3),...
                                          size(t10, 1),    size(t10, 2),    size(t10, 3));

                       %If there are 2 dimensions, must be a vector of values:
                       elseif length(temp_compare) == 1 && (~isequal(temp_compare(1), 1) ||~isequal(temp_compare(2), 1))
                            ERRMSG = ['Error attempting to add new tissue: when assigning T10, ktrans, ve and vp values, \n'...
                                      'if the number of dimensions are 2, input must be a vector.  \n'...
                                      'The size of the inputs are: \n'...
                                      'size(ktrans)=[%d, %d, %d] \n'... 
                                      'size(ve)    =[%d, %d, %d] \n'...
                                      'size(vp)    =[%d, %d, %d] \n'...
                                      'size(t10)   =[%d, %d, %d]']; 

                           error(ERRMSG, size(ktrans, 1), size(ktrans, 2), size(ktrans, 3),...
                                          size(ve, 1),     size(ve, 2),     size(ve, 3),...
                                          size(vp, 1),     size(vp, 2),     size(vp, 3),...
                                          size(t10, 1),    size(t10, 2),    size(t10, 3));
                       end

                       %If we make it this far, we must have a value of a valid
                       %size.  Compare it to previous values to make sure they
                       %match.
                       if isempty(temp_compare) || isequal(temp_compare, alternate_size)
                           continue;
                       elseif isequal(alternate_size, [1,1])
                           alternate_size = temp_compare;
                       else
                           ERRMSG = ['Error attempting to add new tissue: when assigning T10, ktrans, ve and vp values, \n'...
                                      'all variables must either contain a single value or must be the same size.  \n'...
                                      'The size of the inputs are: \n'...
                                      'size(ktrans)=[%d, %d, %d] \n'... 
                                      'size(ve)    =[%d, %d, %d] \n'...
                                      'size(vp)    =[%d, %d, %d] \n'...
                                      'size(t10)   =[%d, %d, %d]']; 

                           error(ERRMSG, size(ktrans, 1), size(ktrans, 2), size(ktrans, 3),...
                                          size(ve, 1),     size(ve, 2),     size(ve, 3),...
                                          size(vp, 1),     size(vp, 2),     size(vp, 3),...
                                          size(t10, 1),    size(t10, 2),    size(t10, 3));


                       end

                    end %End for loop comparing property sizes
                
            
                    %If we have made it this far, we have a valid heterogeneous
                    %tissue to add:
                
                    TissueLayerAdds = obj.heterogeneousTissue(the_mask, ktrans, ve, vp, t10, other_input);
                end
            end
            
            
            %Verify that all layers have a T1 assigned before proceeding:
            for counter = 1:length(TissueLayerAdds)
                if isempty(TissueLayerAdds{counter}.t10)
                    TissueLayerAdds{counter}.t10 = obj.defaultT10;
                end
            end
            
            
                    
        end %End preProcessAddTissue

        function [obj] = setEmptyModelParams(obj, theFieldname)
            %If there is not an existing field, add it to all layers.  Set
            %defult value to empty;
            for counter = 1:length(obj.tissueLayers)
                obj.tissueLayers{counter}.setModelParams(theFieldname, []);
            end
        end
        
        function [obj, TissueLayer] = modelParamCheck(obj, TissueLayer)
            
                %First, add any fields from this tissue to the
                %obj.modelParams:
                the_fields = fieldnames(TissueLayer.modelParams);
                for counter2 = 1:length(the_fields)
                    
                    %If its not an existing field and its not empty, add it
                    %to all layers:
                    if ~isfield(obj.modelParams, the_fields{counter2}) && ~isempty(TissueLayer.modelParams.(the_fields{counter2}))
                        obj = obj.setEmptyModelParams(the_fields{counter2});
                    %If it is not an existing field and it is empty, remove
                    %it from the modelparams structure
                    %elseif ~isfield(obj.modelParams, the_fields{counter2}) && isempty(TissueLayerAdds{counter}.modelParams.(the_fields{counter2}))
                     %   TissueLayer.modelParams = rmfield(TissueLayer.modelParams, the_fields{counter2});
                    end
                end
                
                %Next, add any fields that exist in the current obj but not
                %in the tissue layer:
                obj_fields = fieldnames(obj.modelParams);
                for counter2 = 1:length(obj_fields)
                    if ~isfield(TissueLayer.modelParams,obj_fields{counter2})
                        TissueLayer.modelParams.(obj_fields{counter2}) = [];
                    end
                end
        end
        
    end %End protected methods
    
end
