classdef Phantom < handle
% classdef for Phantom object that contains the image space representation 
% of the object to be simulated as well as properties needed by the digital
% phantom simulator.  

    
    %Properties
    properties 
        
        %Vascular Input Function - the vascular input function 
        %  This can either be a matrix with 2 columns (time and
        %  concentration of gadolinium in Moles ([Gd+], M) or a function
        %  handle to a function that generates these values
        VIF = [];
        
        
        xyzRes = [];

        %  patientPosition defines the orientation of the images
        %  that make up this phantom object.  Valid options for
        %  patientPosition include:
        %
        %       'HFS' - Head First Supine
        %       'HFP' - Head First Prone
        %       'FFS' - Feet First Supine
        %       'FFP' - Feet First Prone
        %       'HFDR' - Head First-Decubitus Right
        %       'HFDL' - Head First-Decubitus Left
        %       'FFDR' - Feet First-Decubitus Right
        %       'FFDL' - Feet First-Decubitus Left
        %
        %  If both patientPosition and plane are set, the orientation will 
        %  also be set.
        %  
        patientPosition = '';
        
         
        % plane sets the imaging plane of the images that make up the
        % phantom object.  Valid choices are: 'axial', 'coronal' and
        % 'sagittal'.  If the orientation property has not been set yet, 
        % setting plane will automatically set orientation following the
        % LPS convention:
        %     Axial:    obj.orientation = [1, 0,  0; 
        %                                  0, 1,  0; 
        %                                  0, 0,  1];
        %     Coronal:  obj.orientation = [1, 0,  0; 
        %                                  0, 0,  1; 
        %                                  0, 1,  1];
        %     Sagittal: obj.orientation = [0, 1,  0; 
        %                                  0, 0,  1; 
        %                                  1, 0,  1];
        %
        %plane = '';
        
        % orientation is a matrix defining the 'x', 'y' and 'z' of the 
        % images making up the phantom according to the LPS convention
        % (origin is the top left corner of the image):
        %
        %  [R] - Right - The direction in which X decreases. 
        %  [L] - Left - The direction in which X increases. 
        %  [A] - Anterior - The direction in which Y decreases. 
        %  [P] - Posterior - The direction in which Y increases.
        %  [I] - Inferior - The direction in which Z decreases. 
        %  [S] - Superior - The direction in which Z increases.
        %
        %   If the "Patient" is in the Feet-First-Prone position:
        %   For an Axial scan, you could define to following:
        %      orientation =  [-1,  0, 0;     -> x direction L to R
        %                       0, -1, 0;     -> y direction is P to A
        %                       0   0  1];    -> z direction is I to S
        %    For a breast patient, the would mean the images are oriented
        %    such that the breasts are pointint downwards.    
        %
        % This value iwll be set automatically if obj.plane is 
        % set and obj.orientation is still empty. However, changing 
        % obj.plane after obj.orientation has been set will not overwrite
        % obj.orientation.  
        %
        %orientation=[];
        
        partialVolume_Layer = true;
        
        partialVolume_Tissue = true;      
        
    end
    
    properties (Dependent)

        %xyzDims indicates the size in the x, y and z directions of the s0
        %images.  All subsequently added masks must match this size.
        xyzDims=[];
        
        
                
    end
    
    properties (SetAccess = 'protected')
       
       %Array to hold different tissues making up the phantom:
       theTissues = {};
       
       %customizable user parameter array:
       userParams = struct;
       
    end
    
    properties (SetAccess = 'protected', Dependent)
        
       xyzFOV = [];
       
       tissueNames = {};
       
    end
    
    properties (Dependent, Hidden)
       
        Dims = [];
        
    end
    
    properties (SetAccess = 'protected', Hidden)
       
        tissueTypes = struct;
        
    end
    
    
     % ----------------- Class Constructor ---------------------------%
    methods
        function obj = Phantom(varargin)
            
            %Input notes: first input must be initialTissue, a Tissue
            %object or a cell array of inputs used to make the object;
            
            % Do nothing with zero inputs...
            if ~nargin
                return
            end

            initialTissue = varargin{1};

            varargin(1) = [];
            
            if ~isa(initialTissue, 'Tissue')
                 if iscell(initialTissue)
                    obj.theTissues(1) = {Tissue(initialTissue{1:length(initialTissue)})};
                 else
                    obj.theTissues(1) = {Tissue(initialTissue)};
                 end
            else
                obj.theTissues(1) = {initialTissue};
            end
            
            %Record the Tissue type for later use:
            theFieldName = obj.theTissues{1}.tissueType;
            obj.tissueTypes.(theFieldName) = [1];
            
            %Loop throup properties and set as appropriate:
            the_properties = properties(Phantom);
            counter = 1;
            while counter <= length(varargin)
                thevar = varargin{counter};
                foundProp = false;   
                %If the input is a string, figure out what property it
                %corresponds to and call that set function:
                if ischar(thevar) 
                    for prop_num = 1:length(the_properties)
                        if strcmp(the_properties{prop_num}, thevar)
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
            
            %If no patientPosition was set, attempt to set from qt_image
            %object:
            if isempty(obj.patientPosition)
                
                %Check to see if it is an exisiting field in the qt_image
                %metadata ( i.e. if the qt_object was made from a set of DICOM
                %images)
                if isfield(obj.theTissues{1}.s0_qtimage.metaData, 'patientPosition')
                
                    the_string = obj.theTissues{1}.s0_qtimage.metaData.patientPosition;
                    
                    %Valid positions:
                    %       'HFS' - Head First Supine
                    %       'HFP' - Head First Prone
                    %       'FFS' - Feet First Supine
                    %       'FFP' - Feet First Prone
                    %       'HFDR' - Head First-Decubitus Right
                    %       'HFDL' - Head First-Decubitus Left
                    %       'FFDR' - Feet First-Decubitus Right
                    %       'FFDL' - Feet First-Decubitus Left
                    
                    valid_string = {'HFS', 'HFP', 'FFS', 'FFP', 'HFDR', 'HFDL', 'FFDR', 'FFDL'};
                    obj.patientPosition = the_string;
                    
                    try
                        validatestring(the_string, valid_string);
                    catch 
                        warning('Unrecognized patient position in phantom metadata.  Setting patient position to FFP');
                        obj.patientPosition = 'FFP';
                    end
                else
                    obj.patientPosition = 'FFP';
                end
            end
            
            %Initialize the Resolution
            if isempty(obj.xyzRes)
                obj.xyzRes = obj.theTissues{1}.xyzRes;
            end
            
            
        end
    end
    
     % ----------------- Set Methods ---------------------------%
    methods
        
        function set.xyzDims(obj, thesize)
            
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
            
            %Resize all Tissue objects::
            for counter = 1:length(obj.theTissues)
                obj.theTissues{counter}.xyzDims = thesize;
            end
            
        end %End set.xyzDims
        
        function set.xyzRes(obj, theres)
            
            validateattributes(theres, {'numeric'}, {'numel', 3, 'positive', 'nonnan', 'nonempty', 'vector'});
            
            %if value is unchanged, return
            if isequal(theres, obj.xyzRes)
                return;
            end
            
            %Set the element spacing for the S0 qt_image object:
            obj.xyzRes = theres;
            
            %Reset the resolution of the Tissue objects residing in the
            %phantom:
            for counter = 1:length(obj.theTissues)
                obj.theTissues{counter}.xyzRes = obj.xyzRes;
            end

        end

        function set.patientPosition(obj, str)
         %  patientPosition defines the orientation of the 
         %  images that make up this phantom object.  Valid options for STR
         %  include:
         %
         %       'HFS' - Head First Supine
         %       'HFP' - Head First Prone
         %       'FFS' - Feet First Supine
         %       'FFP' - Feet First Prone
         %       'HFDR' - Head First-Decubitus Right
         %       'HFDL' - Head First-Decubitus Left
         %       'FFDR' - Feet First-Decubitus Right
         %       'FFDL' - Feet First-Decubitus Left
         %

            
            if isempty(str) 
                obj.patientPosition = 'FFP';
                return;
            end
         
            %Validate the input:
            valid_strings = {'HFS', 'HFP', 'FFS', 'FFP', 'HFDR', 'HFDL',...
                'FFDR', 'FFDL'};
            
            validatestring(str, valid_strings);
            obj.patientPosition = str;
            
%             %If both patientPosition and Plane are set, set orientation:
%             obj = setOrientation(obj);

        end %Phantom.set.patientPosition
        
        function set.VIF(obj, vif)
            
            if isa(vif,'function_handle') 
                %get the values from the function:
                thevif = vif;
            else
                thevif = vif;
            end

            %Initial validation of general matrix properties:
            validateattributes(thevif, {'numeric'}, {'ncols', 2, 'nonnan', 'nonempty'});
            
            %Validate the time and concentration columns individually:
            vif_time = vif(:,1);
            vif_conc = vif(:,2);

            validateattributes(vif_time, {'numeric'}, {'vector', 'nonnan', 'nonempty', 'finite', 'nondecreasing'});
            validateattributes(vif_conc, {'numeric'}, {'vector', 'nonnan', 'nonempty', 'finite', 'nonnegative'});    
            
            obj.VIF = vif;
        end
            
        
    end %End Set Methods
    
    % ----------------- Get Methods ---------------------------%
    methods
        
        function thexyzDims = get.xyzDims(obj)
            
            thexyzDims = obj.theTissues{1}.xyzDims;
            
        end
        
        function theDims = get.Dims(obj)
           theDims = [obj.xyzDims(2), obj.xyzDims(1), obj.xyzDims(3)]; 
        end
        
        function thexyzFOV = get.xyzFOV(obj)
           
            thexyzFOV = obj.xyzDims .* obj.xyzRes * 1/10;
            
        end
        
        function theTissueNames = get.tissueNames(obj)
            
            theTissueNames = cell([1,length(obj.theTissues)]);
            
            for counter = 1:length(obj.theTissues)
               theTissueNames(counter) = {obj.theTissues{counter}.tissueName};
            end

        end %End get.layerNames
        
    end %End Get Methods
    
    % ----------------- Other Methods ---------------------------%
    methods
        
        function obj = addTissue(obj, addTissue)
            %Add new Tissue. input theTissue must be a Tissue object or a cell array
            %containing the inputs necessary to make one:
            
            %Validate the inputs.  The addTissue must already  be a
            %Tissue object a cell array containing parameters to make one:
            if ~isa(addTissue, 'Tissue')
                 if iscell(addTissue)
                    newTissue = Tissue(addTissue{1:length(addTissue)});
                 else
                    newTissue = Tissue(addTissue);
                 end
            else
                newTissue = addTissue;
            end
            
            %Verify that the matrix size matches existing tissue objects:
            if ~isequal(obj.xyzDims, newTissue.xyzDims)
                ERRMSG = ['ERROR:  The Tissue object being added to the Phantom object has \n' ...
                            ' different xyz dimensions than the those in the Phantom object. %s'];
                        
                error(ERRMSG, '');
            end
            
            %Set the resolution to match the resolution of the phantom:
            newTissue.xyzRes = obj.xyzRes;            
            
            %Add the tissue to the array of tissue objects:
            thePosition = length(obj.theTissues) + 1;
            obj.theTissues(thePosition) = {newTissue};
            
            %Record the tissue type for later use:
            theFieldName = obj.theTissues{thePosition}.tissueType;
            if ~isfield(obj.tissueTypes, theFieldName)
                obj.tissueTypes.(theFieldName) = thePosition;
            else
                obj.tissueTypes.(theFieldName) = [obj.tissueTypes.(theFieldName), thePosition];
            end
            
        end
        
        function obj = removeTissue(obj, number)
          
            %First validate the input:
            validateattributes(number, {'numeric'}, {'scalar', 'positive', '<=', length(obj.theTissues), 'nonempty'});
            if ~isequal(ceil(number), floor(number))
                ERRMSG = ['ERROR: Expected the layer number to be non-decminal number'];
                error(ERRMSG);
            end
            
            if length(obj.theTissues) == 1
                ERRMSG = ['ERROR: Can not remove the only Tissue object from the Phantom'];
                error(ERRMSG);
            end
            
            %Remove that tissue from tracking:
            theFieldname = obj.theTissues{number}.tissueType;
            fieldValues = obj.tissueTypes.(theFieldname);
            fieldValues(fieldValues == number) = [];
            
            if isempty(fieldValues)
                obj.tissueTypes = rmfield(obj.tissueTypes, theFieldname);
            else
                fieldValues(fieldValues > number) = fieldValues(fieldValues > number) - 1;
                obj.tissueTypes.(theFieldname) = fieldValues;
            end
            
            %Remove the desired element from the cell array:
            obj.theTissues(number) = [];
            
        end

        function tissueLesion = generateLesion(obj, tissueType, lesionShape, varargin)
            % Inputs:
            %   tissueType - the tissueType that the lesion should be
            %      inserted into
            %   lesionShape - this is either a matrix providing a mask to
            %      define the lesion shape or it is a string specifying which
            %      shape lesion to generate
            %   lesionLoc - relative location of lesion in xyz.  values range from -1 to 1.  
            %   varargin - additional optional inputs used to customze the
            %      lesion shape
          
            
            %Validate tissueType:
            valid_tissueTypes = fieldnames(obj.tissueTypes);
            try
                validatestring(tissueType, valid_tissueTypes);
            catch
                if length(valid_tissueTypes) > 1
                    format_string = [];
                    for counter = 1:length(valid_tissueTypes)
                        format_string = [format_string, '   ''%s'' \n'];
                    end
                    ERRMSG = ['ERROR: unrecognized lesion tissue type.  Options for lesion tissue type those contained in the current Phantom object: \n', format_string];
                    error(ERRMSG, valid_tissueTypes{1:length(valid_tissueTypes)});
                else
                    ERRMSG = ['ERROR: unrecognized tissue type.'];
                    error(ERRMSG);
                end
            end
            
            %Create the lesion object: 
             tissueLesion = Lesion(lesionShape, 'tissueType', tissueType, varargin{1:end});
            
        end
        
        function [obj, lesionLoc_relxyz] = insertLesion(obj, theLesion, varargin)
            %Inserts a lesion (contained in a Tissue object) into the
            %Digital Phantom.  Use 'obj.generateLesion' to generate an
            %appropriate tissue object.
            %
            % theLesion  - Tissue object containing the lesion for
            % insertion
            %
            % tissueNumber - optional argument specifying which tissue
            % object in the phantom the lesion should be inserted into.
            %
            % location - optional argument specifying the location in the
            % phantom the lesion is to be inserted at.  Location is the
            % relative location with the center point of the phantom
            % located at [0, 0, 0]  and the edges at +/- 1.
            
            %Validate the inputs:
            if ~isa(theLesion, 'Lesion')
                ERRMSG = ['ERROR: the lesion object is not a Lesion object. Use the \n'...
                            'Phantom.generateLesion method to generate an appropriate lesion \n'... 
                            'or create your own using the Lesion class. %s'];
                error(ERRMSG, '');
            end
            
            %Sort throught the varargin to look for parameters defining
            %lesion parameters.  First set up defaults, then overwrite with 
            %contents of varargin as appropriate:
            valid_options = {'location', 'tissueNumber'};

            lesionLoc = [];
            lesionSize = theLesion.lesionSize; 
            lesionMeanIntensity = theLesion.lesionMeanIntensity;
            lesionStdev = theLesion.lesionStdev;
            tissueNumber = [];

            counter = 1;
            while counter <= length(varargin)
                if ischar(varargin{counter})
                    
                    %Skip it if its not one of the options required by this function:
                    try 
                        validatestring(varargin{counter}, valid_options);
                    catch
                        counter = counter + 1;
                        continue;
                    end
                    
                    
                    if strcmp('location', varargin{counter})
                        lesionLoc = varargin{counter + 1};
                        validateattributes(lesionLoc, {'numeric'}, {'numel', 3, 'vector', '<=', 1, '>=', -1});
                        lesionLoc = lesionLoc(:)';
                        counter = counter + 1;
                    elseif strcmp('tissueNumber', varargin{counter})
                        validateattributes(tissueNumber, {'numeric'}, {'scalar', 'positive', '<=', length(obj.theTissues)});
                
                        if ~strcmp(obj.theTissues{tissueNumber}.tissueType, theLesion.tissueType)
                            ERRMSG = ['ERROR: The lesion object is not the same tissue type as the \n'...
                              'tissue it is being inserted into: \n'...
                              'Lesion tissueType: %s \n'...
                              'Existing tissueType: %s \n'];
                            error(ERRMSG, theLesion.tissueType, obj.theTissues{tissueNumber}.tissueType);
                        end
                    end
                end
                counter = counter + 1;
            end
            
            %Determine the tissue the lesion is to be inserted into:
            if isempty(tissueNumber)
                tissueNumber = obj.tissueTypes.(theLesion.tissueType)(1);
            end
            
            % Determine the desired lesion location.  If lesionLoc is 
            % empty, call GUI to pick location:
            if isempty(lesionLoc)
                lesionLoc = Lesion_Placement_GUI(obj.theTissues{tissueNumber}.s0);
                
                lesionLoc_relxyz = 2*lesionLoc./obj.xyzDims - 1;
                
                the_fig = gcf;
                close(the_fig);
            else
                lesionLoc_relxyz = lesionLoc;
                lesionLoc = round(lesionLoc .* [obj.xyzDims]*1/2 + [obj.xyzDims]*1/2); 
            end
            
            %Resize the lesion to the desired size:
            %desired_xyzDims = round(theLesion.xyzRes .* theLesion.xyzDims .* obj.xyzRes.^-1);

            theLesion = theLesion.partialVolumeLesion(obj.xyzRes);
            
            %Place the lesion at the desired location:
            theLesion = theLesion.placeAtLocation(lesionLoc, obj.xyzDims);
            
            %If the Lesion s0 has not been set, add the fill:          
            %First, obtain desired mean intensity and standard deviation
            %for the voxels making up the lesion:
            if ~theLesion.s0_set
                if isempty(lesionMeanIntensity) || isempty(lesionStdev)

                    %Base the mean and standard deviation off the center slice 
                    %the lesion is to be inserted in:
                    lesion_locations = find(theLesion.noiseMask > 0);
                    the_values = zeros([length(lesion_locations), 1]);
                    noise_level = 0;
                    central_slice = zeros([obj.Dims(1), obj.Dims(2)]);
                    for counter = 1:length(obj.tissueTypes.(theLesion.tissueType))
                        central_slice = central_slice + obj.theTissues{obj.tissueTypes.(theLesion.tissueType)(counter)}.s0(:,:, lesionLoc(3));
                        the_values = the_values + abs(obj.theTissues{obj.tissueTypes.(theLesion.tissueType)(counter)}.s0(lesion_locations));
                        noise_level = noise_level + obj.theTissues{obj.tissueTypes.(theLesion.tissueType)(counter)}.s0_threshold;
                        central_slice = central_slice + abs(obj.theTissues{obj.tissueTypes.(theLesion.tissueType)(counter)}.s0(:,:, lesionLoc(3)));
                    end
                    noise_level = noise_level / length(obj.tissueTypes.(theLesion.tissueType));

                    %assign the mean and std based on the values greater than the
                    %nosise floor contined in the mask:
                    values_for_calc = the_values(the_values > noise_level);

                    if isempty(lesionMeanIntensity)
                        %If at least a third of the values are above the noise
                        %thershold, use them to determine the mean:
                        if length(values_for_calc) > length(the_values)/3
                            lesionMeanIntensity = mean(values_for_calc(:));
                        else 
                            values_for_calc = central_slice(central_slice > noise_level);
                            lesionMeanIntensity = mean(values_for_calc(:));
                        end
                    end

                    if isempty(lesionStdev)

                        %If at least a third of the values contained in the
                        %mask are greater than the noise floor, use the
                        %standard deviation of those voxels:
                         if length(values_for_calc) > length(the_values)/3
                             lesionStdev = std(values_for_calc);
                         else

                             %Look at 10x10 tiles in the central slice:
                             the_stdev = zeros([floor(obj.Dims(1)/10)*floor(obj.Dims(2)/10),1]);
                             counter = 1;
                             for rowi = 0:floor(obj.Dims(1)/10)-1
                                 for coli = 0:floor(obj.Dims(2)/10)-1
                                     the_tile = central_slice(rowi*10+1:rowi*10+10, coli*10+1:coli*10+10);
                                     tile_values = the_tile(the_tile > noise_level);

                                     if length(tile_values)> numel(the_tile)/3
                                         the_stdev(counter) = std(tile_values);
                                     else
                                         the_stdev(counter) = 0;
                                     end

                                     counter = counter + 1; 
                                 end
                             end

                             the_stdev(the_stdev == 0) = [];
                             lesionStdev = mean(the_stdev);
                         end
                    end
                end

                %Call function to apply lesion fill:
                theLesion = theLesion.fillS0(lesionMeanIntensity, lesionStdev);
            end
            
            orig_layers = length(obj.theTissues{tissueNumber}.tissueLayers);
            
            obj.theTissues{tissueNumber} = obj.theTissues{tissueNumber}.insertTissueObject(theLesion);
            
            %combine "default non enhancing tissue" layers as appropriate:
            obj.theTissues{tissueNumber} = obj.theTissues{tissueNumber}.combineLayers([1, orig_layers + 1]);
            
%             %Replace the s0 in the Tissue object:
%             tissueS0 = obj.theTissues{tissueNumber}.s0;
%             lesionS0 = theLesion.s0;
%             lesion_location = find(theLesion.noiseMask > 0);
%             noise_threshold = obj.theTissues{tissueNumber}.s0_threshold;
%             
%             %If the tissue is made up of complex data and the lesion is
%             %not, assign phase information to it:
%             if obj.theTissues{tissueNumber}.complex_data && ~theLesion.complex_data 
%             
%                %If at least a third of the voxels in the tissue have temp
%                %magnitudes greater than the noise floor, assign the lesion
%                %phase to match that phase, else assign it to 0:
%                all_values = tissueS0(lesion_location);
%                phase_locations = lesion_location(abs(all_values) > noise_threshold);
%                the_values = tissueS0(phase_locations);
%                
%                if length(the_values) > length(all_values)*1/3    
%                    mean_vector = mean(the_values);
%                    mean_phase = angle(mean_vector);
%                    apply_phase = true;
%                else
%                    apply_phase = false;
%                end
%                
%                lesionS0 = complex(lesionS0);
%                if apply_phase == true;
%                    
%                    %For voxels in the tissue with enough signal, assign the
%                    %lesion to have that phase:
%                    lesionS0(phase_locations) = complex(lesionS0(phase_locations).*cos(angle(tissueS0(phase_locations))), lesionS0(phase_locations).*sin(angle(tissueS0(phase_locations))));
%                    
%                    %for voxels without enough signal, assign them to have
%                    %the mean phase:
%                    other_locations = lesion_location(abs(all_values) <= noise_threshold);
%                    lesionS0(other_locations) = complex(lesionS0(other_locations)*cos(mean_phase), lesionS0(other_locations).*sin(angle(mean_phase)));
%                end
%                
%             end
%                
%             %Now, insert the lesionS0 into the tissueS0:
%             tissueS0 = tissueS0 .*(abs(1-theLesion.noiseMask)) + lesionS0;
%             obj.theTissues{tissueNumber}.s0_qtimage = qt_image(tissueS0);
%             
%             %Add the lesion enhancement information:
%             for counter = 2:length(theLesion.tissueLayers)
%                 theLesion = theLesion.setLayerName(counter, [theLesion.tissueName, ' ', theLesion.tissueLayers{counter}.layerName]);
%                 obj.theTissues{tissueNumber} = obj.theTissues{tissueNumber}.insertTissueLayer(theLesion.tissueLayers{counter});
%             end
            

            %Finally, exclude the region of the lesion from other tissues:
            for counter = 1:length(obj.theTissues)
                
                %Skip the Tissue that was already processed:
                if counter == tissueNumber
                    continue;
                end
                
                %Remove the lesion area from the s0 images:
                tissueS0 = obj.theTissues{counter}.s0;
                tissueS0 = tissueS0 .*(abs(1-theLesion.noiseMask));
                
                if obj.theTissues{counter}.complex_data
                    obj.theTissues{counter}.s0_qtimage = qt_image(real(tissueS0));
                    obj.theTissues{counter}.s0_qtimageImag = qt_image(imag(tissueS0));
                else
                    obj.theTissues{counter}.s0_qtimage = qt_image(tissueS0);
                end
                
                %Remove the lesion area from the enhancement maps:
                if length(obj.theTissues{counter}.tissueLayers) > 1
                    obj.theTissues{counter} = obj.theTissues{tissueNumber}.insertTissueLayer(theLesion.noiseMask);
                    obj.theTissues{counter} = obj.theTissues{counter}.combineLayers([1, length(obj.theTissues{counter}.tissueLayers)]);
                end
            end 

        end
          
        function obj = flattenTissue(obj)
            
            for counter = 1:length(obj.theTissues)
                obj.theTissues(counter) = {obj.theTissues{counter}.combineLayers};
            end
        end
        
        function obj = setUserParams(obj, fieldname, value)
            
            %Validate the fieldname:
            if ~isvarname(fieldname)
                ERR_MSG = ['Error in setUserParams: the supplied field name is not a valid matlab variable name'];
                error(ERR_MSG);
            end
            
            obj.userParams.(fieldname) = value;

        end
        
    end
    
    % ----------------- Private Methods ---------------------------%
    methods (Access = 'private')

        function filtered_image = filterFunction(obj, input_images, sigma_scalar)
            
            if nargin < 3
                sigma_scalar = 0.5;
            end
            
            sigma = sigma_scalar * [obj.xyzRes(1), obj.xyzRes(1), obj.xyzRes(1)]./[obj.xyzRes];

            x = -1:1;
            y = -1:1; 
            z = -1:1;

            [X, Y, Z] = meshgrid(x, y, z);

            the_gaussian = exp(-(X.^2/(2*sigma(1)^2) + Y.^2/(2*sigma(2)^2) + Z.^2/(2*sigma(3)^2)));
            the_gaussian = the_gaussian * 1/sum(the_gaussian(:));
            
            filtered_image = imfilter(input_images, the_gaussian);
            
        end
        
        
    end %End Private Methods
    
end
