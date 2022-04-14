classdef Lesion < Tissue
% classdef for Lesion objects.  Lesions are a subclass of Tissue.

    
    %Properties
    properties (SetAccess = 'protected')
        %properties related to building lesion morphology.  Must be
        %provided in constructor otherwise they will be assigned default
        %values
        lesionShape = '';
        lesionFeatures =[];
        
        randSeed = 0;
        lesionFeatureSize = [];
        
    end
    
    properties
       %Properties that are changable after initial construction
        
        lesionSize = 3;  %lesion size in cm
        lesionMeanIntensity = [];
        lesionStdev = [];
        
    end
    
    properties(SetAccess = {?Lesion, ?Phantom}, Hidden)
       
        s0_set = false;
        
    end
    
    properties (SetAccess = 'protected', Hidden)
        
        lesionBaseSize = 5;
        buildRequired = true;
    end
    
 % ----------------- Class Constructor ---------------------------%
    methods
        function obj = Lesion(thelesionShape, varargin)
            
            %Recognized optional inputs:
            % *All valid properties of Tissue class
            % 

            isS0_set = false;    
            
            % Do nothing with zero inputs...
            if nargin == 0
               TissueInputs = {};
               
            else
                
                TissueInputs = {};
                
                %First, create the underlying Tissue object.  Search the
                %inputs for anything labled 's0'
                counter = 1;
                while counter <= length(varargin)
                    if ischar(varargin{counter})
                        if strcmp(varargin{counter}, 's0')
                        
                            TissueInputs(1) = {varargin{counter + 1}};
                                        
                            %Remove s0 matrix/object from the list of variable
                            %inputs:
                            varargin(counter + 1) = [];
                        
                            %Remove s0 string from the list of variable inputs:
                            varargin(counter) = [];
                            
                            counter = counter - 1;
                        
                            break;
                        end
                    end
                    counter = counter + 1;
                end
                
                %Determine if a mask was sent in for thelesionShape
                %variable:
                if ~ischar(thelesionShape) && ismatrix(thelesionShape)
                    %If the user sent in a mask defining the lesion location,
                    %validate its 
                    validateattributes(thelesionShape, {'numeric'}, {'nonnegative', '<=', 1});
                    lesion_mask = thelesionShape;
                    
                    %crop the lesion mask such that there is a single row of
                    %zeros around the edges:
                    [row, col, slice] = ind2sub(size(lesion_mask), find(lesion_mask > 0));
                
                    newDims = [max(row)-min(row) + 3, max(col)-min(col) + 3, max(slice)-min(slice)+3];
                    temp_lesionMask = lesion_mask;
                    lesion_mask = zeros(newDims);
                
                    lesion_mask(2:end-1, 2:end-1, 2:end-1) = temp_lesionMask(min(row):max(row), min(col):max(col), min(slice):max(slice));
                    clear temp_lesionMask;
                    
                    %resize s0 to match if appropriate:
                    if exist('TissueInputs') && ~isempty(TissueInputs{1})
                       news0 = zeros(newDims);
                       news0(2:end-1, 2:end-1, 2:end-1) = TissueInputs{1}(min(row):max(row), min(col):max(col), min(slice):max(slice));
                       TissueInputs{1} = news0;
                       clear news0;
                    end
                    
                    %Put mask into appropriate location for the Tissue
                    %class input:
                    TissueInputs(2) = lesion_mask;
                    thelesionShape = 'custom';
                    
                    clear lesion_mask;

                end
                
                if ~isempty(TissueInputs)
                    if isempty(TissueInputs{1}) && ~isempty(TissueInputs{2})
                        %If no s0 was set by the user but a desired lesion mask
                        %was provided, set s0 to a matrix of 0s.
                        TissueInputs(1) = zeros(size(TissueInputs{2}));
                        
                        isS0_set = false;
                    else
                        isS0_set = true;
                    end
                end
                
                TissueInputs(3:length(varargin)+2) = varargin;
                
            end
            
            %If a lesion shape was specified, set up inputs as temporary
            %placeholders in initial tissue object input as appropriate:
            if exist('thelesionShape', 'var') && (isempty(TissueInputs) || isempty(TissueInputs{1}))
                TissueInputs{1} = zeros([1, 1, 1]);
            end
            
            if exist('thelesionShape', 'var') && length(TissueInputs) >1 && isempty(TissueInputs{2})
                TissueInputs(2) = [];
            end
            
            %Call the Tissue constructor:
            obj@Tissue(TissueInputs{:});
                        
            if nargin == 0
                return;
            end
            
            obj.s0_set = isS0_set;
            
            %Now customize the tissue object into the desired lesion
            
            %Set the lesion shape:           
            obj.lesionShape = thelesionShape;
            
            %Call the function to generate the lesion:              
            switch thelesionShape
                case 'round'
                    obj.lesionFeatures = 0;
                    obj.lesionFeatureSize = 0;
                case 'lobulated'
                    obj.lesionFeatures = 3;
                    obj.lesionFeatureSize = obj.lesionBaseSize;
                case 'irregular'
                    obj.lesionFeatures = 15;
                    obj.lesionFeatureSize = obj.lesionBaseSize*.4;
                case 'spiculated'
                    obj.lesionFeatures = 15;
                    obj.lesionFeatureSize = obj.lesionBaseSize;
            end
            
            %Finish sorting through the variable inputs:
                
            %Assemble a list of the lesion specific parameters:
            full_prop = properties(Lesion);
            tissue_prop = properties(Tissue);

            lesion_prop = {};

            for counter1 = 1:length(full_prop)
                found_prop = 0;
                for counter2 = 1:length(tissue_prop)

                    if strcmp(full_prop{counter1}, tissue_prop{counter2})
                        found_prop = 1;
                        break;
                    end
                end

                if found_prop == 0
                    %Remove lesionShape from the list as it has already been
                    %checked
                    if ~strcmp(full_prop(counter1), 'lesionShape')
                        lesion_prop(length(lesion_prop) + 1) = full_prop(counter1);
                    end
                end
            end
                
            %Identify and set lesion properties as appropriate:
            counter = 1;
            while counter <= length(varargin)

                if ischar(varargin{counter})
                    found_prop = 0;

                    for counter2 = 1:length(lesion_prop)
                        if strcmp(varargin{counter}, lesion_prop{counter2})
                            found_prop = 1;
                            
                            obj.(lesion_prop{counter2}) = varargin{counter + 1};

                            varargin(counter + 1) = [];
                            varargin(counter) = [];

                            break;

                        end
                    end

                    if found_prop == 0
                        counter= counter + 1;
                    end

                else
                    counter = counter + 1;
                end
            end

            %Call the function to customize the lesion:
            obj = obj.buildLesion;
            

        end
    end
    
  % ----------------- Set methods ---------------------------%  
    methods

        function obj = set.lesionShape(obj, thestr)
           
            if ~ischar(thestr)
                ERRMSG=(['ERROR setting lesionShape: The lesion shape must be a string.']);
                error(ERRMSG);
            end
                
            %If lesion shape was set to custom in the constructor, do not
            %allow for changes:
            if strcmp(obj.lesionShape, 'custom')
                ERRMSG=(['ERROR setting lesionShape: The lesion shape was set to custom and can not be changed']);
                error(ERRMSG);
            elseif strcmp(obj.lesionShape, thestr)
                %if the lesion shape is not changing, do noting
                return;
            else
                % Validate lesion shape: 
                valid_lesionShapes = {'round', 'lobulated', 'irregular', 'spiculated', 'custom'};
                try
                    validatestring(thestr, valid_lesionShapes);
                catch
                    if length(valid_lesionShapes) > 1
                        format_string = [];
                        for counter = 1:length(valid_lesionShapes)
                            format_string = [format_string, '   ''%s'' \n'];
                        end
                        ERRMSG = ['ERROR: unrecognized lesion type.  Options for lesion type are: \n', format_string];
                        error(ERRMSG, valid_lesionShapes{1:length(valid_lesionShapes)});
                    else
                        ERRMSG = ['ERROR: unrecognized lesion shape.'];
                        error(ERRMSG);
                    end
                end
                
                obj.lesionShape = thestr;
                %obj.buildRequired = true;
            end
        end
        
        function obj = set.lesionSize(obj, thesize)
        %Lesion size in cm
            
            if isequal(obj.lesionSize, thesize);
                return;
            end
        
            %Validate the input:
            validateattributes(thesize,{'double'}, {'nonempty', 'size', [1, 1], 'positive', 'nonnan', 'finite'}); 
            obj.lesionSize = thesize;
            
            %Set obj.xyzRes if appropriate
            if isequal(obj.xyzDims, [1, 1, 1])
                obj.xyzRes = [thesize, thesize, thesize];
            elseif isequal(obj.xyzDims, [1, 1, 1])
                [max_val, dim] = max(obj.xyzDims);
                scalar = obj.lesionSize/(max_val-2) * 10;
                obj.xyzRes = obj.xyzRes * scalar;
            else
                [max_val, dim] = max(obj.xyzDims);
                scalar1 = obj.lesionSize/(max_val-2) * 10;
                scalar2 = scalar1/obj.xyzRes(dim);
                obj.xyzRes = obj.xyzRes * scalar2;
            end
            
            
            %Set flag requiring lesion build
            %obj.buildRequired = true;
            
        end
        
        function obj = set.lesionFeatures(obj, numFeatures)
        %Number of features (ie lobules, bumps or spicules) on the lesion
        
            if isequal(obj.lesionFeatures, numFeatures);
                return;
            end
        
            %Validate the input:
            validateattributes(numFeatures,{'double'}, {'size', [1, 1], 'nonnegative'}); 
            obj.lesionFeatures = numFeatures;
  
            %Set flag requiring lesion build
            %obj.buildRequired = true;
            
        end
        
        function obj = set.randSeed(obj, theseed)
        %Set the seed for the random number generator
            rng(theseed);
            obj.randSeed = theseed;
        end
        
        function obj = set.lesionFeatureSize(obj, thesize)
        %the lesion size in cm
        
            if isequal(obj.lesionFeatureSize, thesize);
                return;
            end
                   
            %Validate the input:
            validateattributes(thesize,{'double'}, {'nonnegative', 'size', [1, 1], 'finite', 'nonnan'}); 
            obj.lesionFeatureSize = thesize;
  
            %Set flag requiring lesion build
            %obj.buildRequired = true;
            
        end
              
        function obj = set.lesionMeanIntensity(obj, meanIntensity)
            validateattributes(meanIntensity,{'numeric'}, {'size', [1, 1], 'finite', 'nonnan'}); 
            obj.lesionMeanIntensity = meanIntensity;
            
        end
        
        function obj = set.lesionStdev(obj, thestdev)
            validateattributes(thestdev,{'double'}, {'size', [1, 1], 'finite', 'nonnan', 'nonnegative'}); 
            obj.lesionStdev = thestdev;
            
        end
        
        
    end
    
 % ----------------- Other Methods ---------------------------%     
    methods
       
        function obj = buildLesion(obj)
        %Function to build the lesion using built in functions 
        
            if ~strcmp(obj.lesionShape, 'custom')
        
                lesionMask = LesionGenerator(obj);
                
                %crop the lesion mask such that there is a single row of
                %zeros around the edges:
                [row, col, slice] = ind2sub(size(lesionMask), find(lesionMask > 0));
                
                newDims = [max(row)-min(row) + 3, max(col)-min(col) + 3, max(slice)-min(slice)+3];
                temp_lesionMask = lesionMask;
                lesionMask = zeros(newDims);
                
                lesionMask(2:end-1, 2:end-1, 2:end-1) = temp_lesionMask(min(row):max(row), min(col):max(col), min(slice):max(slice));
                
                clear temp_lesionMask;
                
                %set S0:
                obj.s0_qtimage = qt_image(zeros(size(lesionMask)));
                
                %Set the lesion mask:
                obj.noiseMask = lesionMask;
                tempT1 = obj.tissueLayers{1}.t10;
                tempT2 = obj.tissueLayers{1}.t2;
                obj.tissueLayers(1) = {TissueLayer(obj.noiseMask_qtimage, 'layerName', obj.layerNames{1}, 't10', tempT1, 't2', tempT2)};

                %Apply size criteria:
                [max_val, dim] = max(obj.xyzDims);
                
                new_resolution = (obj.lesionSize*10)/(max_val-2); %Multiply by 10 to change units from cm to mm
                obj.xyzRes = [new_resolution, new_resolution, new_resolution];

            else
                
                if isequal(obj.xyzRes, [1, 1, 1])
                    WRG_MSG = ['Lesion resolution still at default value.  Set intended resolution using obj.xyzRes'];
                    warning(WRG_MSG);
                    
                end
            
            end
            
        end
        
        function obj = partialVolumeLesion(obj, newXyzRes)
        %Function to shrink the lesion down to the appropriate size to fit 
        %in the phantom.  
        %
        % newXyzRes - the desired xyzSize
        
            %First, Validate the input:
            %Validate the requested resolution:
            validateattributes(newXyzRes, {'numeric'}, {'numel', 3, 'positive', 'nonnan', 'nonempty', 'vector'});
            
            %If the new size is the same as the existing size, return:
            if isequal(obj.xyzRes, newXyzRes)
                return;
            end
            
                        
            %Determine the matrix size that will get us closes to the
            %desire resolution:
            newMatrixSize = round(obj.xyzRes ./ newXyzRes .* obj.xyzDims);
            
            %Reset the resolution:
            obj.xyzRes = newXyzRes;
            
            %If any of the inputs are larger than the current matrix size,
            %call the resize function:
            if (newMatrixSize > obj.xyzDims)
                desiredSize = max([thesize; obj.xyzDims], [], 1);
                obj.xyzDims = desiredSize;
            end
            
            
            %Resize s0:
            if obj.s0_set
                obj.s0_qtimage = parVolResize(obj.s0_qtimage, newMatrixSize);
            else
                obj.s0_qtimage = qt_image(zeros([newMatrixSize(2), newMatrixSize(1), newMatrixSize(3)]));
            end
            
            %switch to matlab coordinates
            thesize_mat = [newMatrixSize(2), newMatrixSize(1), newMatrixSize(3)];
            
            %Resize the noise mask:
            newMask = obj.parVolResize(obj.noiseMask_qtimage.value, thesize_mat);
            newMask(newMask < 0) = 0;
            newMask(newMask > 1) = 1;
            obj.noiseMask_qtimage =  qt_image(newMask);
            
            %Resize all masks.  Start from upper most mask.  Ensure the sum
            %of succesesive layers do not exceed the value in the noise
            %along the edges of the noise mask or values will not come out
            %correctly in phantom:
            sumMask = zeros(size(obj.noiseMask));
            store_masks = cell([1, length(obj.tissueLayers)]);
            for counter = length(obj.tissueLayers):-1:1
                newMask = obj.parVolResize(obj.tissueLayers{counter}.mask.value, thesize_mat);
                newMask(newMask < 0) = 0;
                newMask(newMask > 1) = 1;
                
                %Check to see that values are not exceeding those of the
                %noise mask:
                temp_sumMask = sumMask + newMask;               
                
                %Calculate the difference between the sumMask and the noise
                %Mask:
                differenceMask = temp_sumMask - obj.noiseMask;
                
                %Ignore any areas where the noiseMask is 1 as behaviour is
                %correct in these areas:
                differenceMask(obj.noiseMask>=1) = 0;
                
                %Remove values where the the value of the noise mask is not
                %exceeded:
                differenceMask(temp_sumMask <= obj.noiseMask) = 0;
                
                %Generate a new mask that is reduced by the amount the
                %noise mask is exceeded by:
                adjustedMask = newMask - differenceMask;
                sumMask = sumMask + adjustedMask;
                
                %Store mask temporarily so more processing can be done:
                store_mask(counter) = {adjustedMask};
            end
                
            %We are also seeing some areas where we wind up with values
            %less than the noise mask especially in areas around the rim
            %mask.  address this:
            differenceMask = sumMask - obj.noiseMask;
            differenceMask(differenceMask > 0) = 0;
            
            scale_mask = 1./(1-abs(differenceMask));
            
            
            %Make the new layer objects:
            for counter = length(obj.tissueLayers):-1:1
                newLayer = TissueLayer(qt_image(store_mask{counter}.*scale_mask)); 
                thelayerprops = properties(TissueLayer);
                
                for counter2 = 1:length(thelayerprops)
                    if strcmp(thelayerprops{counter2}, 'mask') || strcmp(thelayerprops{counter2}, 'defaultTissueTypes')
                        continue;
                    end
                    newLayer.(thelayerprops{counter2}) = obj.tissueLayers{counter}.(thelayerprops{counter2});
                end
                obj.tissueLayers(counter) = {newLayer};
            end
                        
        end

        function obj = placeAtLocation(obj, desiredxyzLoc, desiredxyzDims)
            
            %If the lesion is larger than the desired new matrix size,
            %return
            if desiredxyzDims < obj.xyzDims
                return;
            end
            
            desiredLoc = [desiredxyzLoc(2), desiredxyzLoc(1), desiredxyzLoc(3)];
            desiredDims = [desiredxyzDims(2), desiredxyzDims(1), desiredxyzDims(3)];
            
            %Re-map s0:
            if obj.s0_set
                obj.s0_qtimage = qt_image(obj.mapLesion(obj.s0, desiredLoc, desiredDims));
            else
                obj.s0_qtimage = qt_image(zeros(desiredDims));
                obj.xyzDims = desiredxyzDims;
            end
            
            %Re-map the noise mask:
            obj.noiseMask_qtimage = qt_image(obj.mapLesion(obj.noiseMask, desiredLoc, desiredDims));
            
            for counter = 1:length(obj.tissueLayers)
                obj.tissueLayers{counter}.mask = qt_image(obj.mapLesion(obj.tissueLayers{counter}.mask.value, desiredLoc, desiredDims));
            end
            
        end
        
        function obj = fillS0(obj, meanIntensity, stdev)
            
            %Determine the locations to add signal:
            the_locations = find(obj.noiseMask>0);
            
            %Initialize the rng:
            rng(0);
            
            %Assign the S0 values:
            lesionS0 = zeros(obj.Dims);
            lesionS0(the_locations) = meanIntensity + stdev * randn([length(the_locations),1]);
            lesionS0(~obj.noiseMask) = meanIntensity;
            
            %Blend so it doesn't look so much like static:
            lesionS0 = obj.filterFunction(lesionS0, 2);
            
            %Apply the noise mask to it:
            lesionS0 = lesionS0 .* obj.noiseMask;
            %lesionS0(obj.noiseMask==0) = 0;
            
            
            %Finally, reset s0:
            obj.s0_qtimage = qt_image(lesionS0);
        end
        
        function obj = generateRim(obj, basemask, thickness)
        %base mask - mask to make the rim around or layer number
        %Thickness - thickness of the rim in mm.  However, if the
        %obj.xyzRes has not been set (i.e. values is [1, 1, 1]), this is 
        %this will be the rim thickness in pixels at the current size.  Rim
        %thickness is the same number of pixels in all directions.
            
        
            %Valiedate the inputs:
            if isscalar(basemask) %The baseMask input indicates a particular tissue layer to use as a template:
                validateattributes(basemask, {'numeric'}, {'scalar', 'positive', 'nonnan', 'nonempty', '<=', length(obj.tissueLayers)});
                original_layer = basemask;
                basemask = obj.getLayerMask(basemask);
            else
                validateattributes(basemask, {'numeric'}, {'size', obj.Dims, 'nonempty', '>=', 0 '<=', 1});
            end

            validateattributes(thickness, {'numeric'}, {'scalar', 'positive', 'nonnan', 'nonempty'});

            thickness = thickness * obj.xyzRes .^-1;

            se_xy = strel('disk', round((thickness(1) + thickness(2))/2));
            se_z = strel('line', round(thickness(3)), 0);

            logical_mask = ceil(basemask);

            %Erode in the xy direction:
            rim_mask = (logical_mask - imerode(logical_mask, se_xy));

            %erode in the z direction:
            the_mask_z = permute(logical_mask, [1, 3, 2]);
            the_mask_z = the_mask_z - imerode(the_mask_z, se_z);
            the_mask_z = permute(the_mask_z, [1, 3, 2]);

            %Put it all together:
            rim_mask = rim_mask | the_mask_z;

            %apply the scaling from the basemask:
            rim_mask = rim_mask.*basemask;

            %Add this as a new layer:
            obj = obj.insertTissueLayer(rim_mask, 'layerName', 'Rim Enhancement');

        end
        
        function displayLesionSurface(obj)

            figure;
            p = patch(isosurface(obj.noiseMask));
            set(p,'FaceColor','red');
            set(p,'EdgeColor','none');
            daspect([1,1,1]);
            axis tight
            view(-38,30)
            camlight
            lighting gouraud
            title('IsoSurface Representation')

        end
        
    end
 
   methods (Access = 'protected')
      
       lesionMask = LesionGenerator(obj)
       
       function newMask = parVolResize(obj, existingMask, thesize)
          
           newMask = zeros(thesize);
           theDims = size(existingMask);
           
           row_lines = linspace(0, theDims(1), thesize(1)+1);
           col_lines = linspace(0, theDims(2), thesize(2)+1);
           slice_lines = linspace(0, theDims(3), thesize(3)+1);

           voxel_vol = row_lines(2)*col_lines(2)*slice_lines(2);
                       
            for row_c = 1:thesize(1)
                for col_c = 1:thesize(2)
                    for slice_c = 1:thesize(3)

                        clear the_tile;

                        the_tile = existingMask(floor(row_lines(row_c))+1:ceil(row_lines(row_c +1)),...
                                              floor(col_lines(col_c))+1:ceil(col_lines(col_c +1)),...
                                              floor(slice_lines(slice_c))+1:ceil(slice_lines(slice_c+1)));

                        %Scale the edges:

                        %row edges:
                        if ~isequal(floor(row_lines(row_c)), ceil(row_lines(row_c)))
                            the_tile(1,:,:) = the_tile(1,:,:) * (1-(row_lines(row_c)-floor(row_lines(row_c))));
                        end
                        if ~isequal(floor(row_lines(row_c+1)), ceil(row_lines(row_c+1)))
                            the_tile(end,:,:) = the_tile(end,:,:) * (1-(ceil(row_lines(row_c+1))-row_lines(row_c+1)));
                        end

                        %col edges:
                        if ~isequal(floor(col_lines(col_c)), ceil(col_lines(col_c)))
                            the_tile(:,1,:) = the_tile(:,1,:) * (1-(col_lines(col_c)-floor(col_lines(col_c))));
                        end
                        if ~isequal(floor(col_lines(col_c+1)), ceil(col_lines(col_c+1)))
                            the_tile(:,end,:) = the_tile(:,end,:) * (1-(ceil(col_lines(col_c+1))-col_lines(col_c+1)));
                        end

                        %slice edges:
                        if ~isequal(floor(slice_lines(slice_c)), ceil(slice_lines(slice_c)))
                            the_tile(:,:,1) = the_tile(:,:,1) * (1-(slice_lines(slice_c)-floor(slice_lines(slice_c))));
                        end
                        if ~isequal(floor(slice_lines(slice_c+1)), ceil(slice_lines(slice_c+1)))
                            the_tile(:,:,end) = the_tile(:,:,end) * (1-(ceil(slice_lines(slice_c+1))-slice_lines(slice_c+1)));
                        end
    
                        newMask(row_c, col_c, slice_c) = sum(the_tile(:))*1/voxel_vol;
                    end
                end
            end
       end
       
       function largeMatrix = mapLesion(obj, smallMatrix, desiredLoc, desiredDims)
           
            %Place the lesion mask at the desired location:
            largeMatrix = zeros([desiredDims]);
            lesionDims = size(smallMatrix);
            
            if numel(lesionDims) < 3
                lesionDims = [lesionDims, 1];
            end
                
            %Determine bounds of lesion:
            bounds = [desiredLoc - floor(lesionDims/2); desiredLoc + ceil(lesionDims/2)-1];
            %Validate bounds
            if min(bounds(1,:)) < 1
                locs = find(bounds(1,:) < 1);
                for counter = 1:length(locs)
                    bounds(:,locs(counter)) = bounds(:,locs(counter))+ (1- bounds(1,locs(counter)));
                end
            end
            if sum(bounds(2,:) > desiredDims)
                locs = find(bounds(2,:) > desiredDims );
                for counter = 1:length(locs)
                    bounds(:,locs(counter)) = bounds(:,locs(counter))- (bounds(2,locs(counter)) - desiredDims(locs(counter)));
                end
            end
        
            largeMatrix(bounds(1,1):bounds(2,1), bounds(1,2):bounds(2,2), bounds(1,3):bounds(2,3)) = smallMatrix;
       end
       
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
       
   end
    
end