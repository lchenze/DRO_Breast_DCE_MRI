classdef MRSystem < handle
% Classdef for the MRSystem object 

    properties (SetAccess = 'protected')
        
        MRSystem_Validated = false;
        
        %MR System: 
        kspacePoints = [];
        
        %customizable user parameter array:
        userParams = struct;
              
    end

    properties 
                        
        %fieldstrength, in Tsla
        fieldStrength = 1.5;
        
        %Trajectory:
        isCartesian = 1;
        is2D = 0; %By default the simulation assumes a 3D slab excitation
        dk = [];
        kmax = [];
        kmin = [];
        kcenter = [];
        
        %Settings:
        bandwidth = 100; %in kHz
        flipAngle = 10; %in degrees
        TR = 10; %in ms
        time = []; %in s
        imageContrastFunc = {};
        
        %Coil:
        sMaps = [];
        
        %Contrast Agent Relaxivity:
        relaxivity = 4.9;
        
        %Settings for ChemSat:
        chemSat_On = 0; %ChemSat turned on - 1 or off - 0
        chemSat_flipAngle = 100; %Flip angle of inversion pulse, in degrees
        chemSat_TR = 100; %in ms
        chemSat_T1Fat = 300; %in ms.  T1 used for MR system calculations.  Actual tissue T1 used in phantom calculaitons.
        chemSat_delayedStart = 10; %in ms.  Delay after chemsat pulse before standard sequence TRs begom
        chemSat_delayedEnd = 0; %in ms.  Dead time at the end of the chemSat TR where no data is acquired.
        
    end
    
    properties (Dependent)
       
       samplingTime = []; %seconds per sample
       nCoil = [];
       finalMatrixSize = [];
       TE = []; %in ms
       
       DAQoff = []; %in ms
       DAQon = []; %in ms
       
       %ChemSat variables:
       chemSat_TI = []; %in ms
       
    end
    
    properties (Hidden, Dependent)
       DAQtime = []; %in ms
       
    end
    
    properties(Access = 'protected')
       
        minkspaceloc = [];
        
        DAQ_start = [];
        DAQ_end = [];
        DAQpercent = .3; %percentage of the TR that is allocated for data aqusion
        
        original_time = [];
        
        override_timeBackUp = 0;
        
    end
    
    events
       timeVectorSet 
       chemSatOn
       chemSatOff
       chemSatValueChange
    end
        
     % ----------- Class Constructor -------------%
    methods
        function obj = MRSystem(kspacePoints, varargin)
            
            % Do nothing with zero inputs...this is required by MATLAB for
            % smooth operation
            if ~nargin
                return
            end
            
            %Add listeners:
            addlistener(obj, 'chemSatOn', @obj.setChemSatTime);
            addlistener(obj, 'timeVectorSet', @obj.setChemSatTime);
            addlistener(obj, 'chemSatOff', @obj.setChemSatTime);
            addlistener(obj, 'chemSatValueChange', @obj.setChemSatTime);
%             addlistener(obj, 'checkChemSatParams', @obj.checkChemSat);
            

            %Validate the input:
            obj.kspacePoints = kspacePoints;
                       
            %Loop through the inputs and set other properties as
            %appropriate:
            the_properties = properties(MRSystem);
            counter = 1;
            while counter <= length(varargin)  
                if ischar(varargin{counter})
                    for prop_counter = 1:length(the_properties)
                        if strcmp(varargin{counter}, the_properties{prop_counter})
                            obj.(the_properties{prop_counter}) = varargin{counter + 1};
                            the_properties(prop_counter) = [];
                            counter = counter + 1;
                            break;
                        end
                    end
                end
                counter = counter + 1;
            end
            
            %Set the image contrast function:
            if isempty(obj.imageContrastFunc)
                obj.imageContrastFunc = @obj.SPGRT1Weighting;
            end
            
            %set kmax based on input kdata:
            if isempty(obj.kmax)
                obj.autoSetKmax;
            end
            
            %set kmin based on input kdata:
            if isempty(obj.kmin)
                obj.autoSetKmin;
            end
            
            if isempty(obj.kcenter)
                obj.autoSetKcenter;
            end
            
            %Set the time vector if it was not set by the user:
            if isempty(obj.time)
                obj.time = [0:obj.TR:(size(obj.kspacePoints, 3)-1)*obj.TR]*10^-3;
            end
                 
           
        end
    end
    
    % ----------- set methods -------------% 
    methods
       
        function set.kspacePoints(obj, thekspacePoints)
            validateattributes(thekspacePoints, {'numeric'}, {'real', 'ndims', 3, 'ncols', 3, 'nonnan', 'nonempty', 'finite'});
            
            obj.kspacePoints = thekspacePoints;
            
            obj.MRSystem_Validated = false;
        end
            
        function set.fieldStrength(obj, thefieldstrength)
            validateattributes(thefieldstrength, {'numeric'}, {'real', 'scalar', 'positive', '<=', 10, 'nonempty'});
            
            obj.fieldStrength = thefieldstrength;
            
        end
        
        function set.isCartesian(obj, thevalue)
            validateattributes(thevalue, {'logical'}, {'scalar', 'nonempty'});
            
            obj.isCartesian = thevalue;
            
        end
        
        function set.is2D(obj, thevalue)
            validateattributes(thevalue, {'logical'}, {'scalar', 'nonempty'});
            
            obj.is2D = thevalue;
            
        end
           
        function set.dk(obj, thedk)
        
            validateattributes(thedk, {'numeric'}, {'real', 'numel', 3, 'nonnan', 'nonempty', 'finite', 'positive'});
            
            obj.dk = thedk;
        end
        
        function set.kmax(obj, thekmax)
            validateattributes(thekmax, {'numeric'}, {'real', 'numel', 3, 'nonnan', 'nonempty', 'finite'});
            
            obj.kmax = thekmax;
        end
        
        function set.kmin(obj, thekmin)
            validateattributes(thekmin, {'numeric'}, {'real', 'numel', 3, 'nonnan', 'nonempty', 'finite'});
            
            obj.kmin = thekmin;
        end
        
        function set.kcenter(obj, thekcenter)
            validateattributes(thekcenter, {'numeric'}, {'real', 'numel', 3, 'nonnan', 'nonempty', 'finite'});
            
            obj.kcenter = thekcenter;
        end
            
        function set.bandwidth(obj, thebandwidth)
           
            validateattributes(thebandwidth, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
            
            obj.bandwidth = thebandwidth;
        end
        
        function set.flipAngle(obj, theflipAngle)
            
            validateattributes(theflipAngle, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', '<=', 180, '>=', -180});
           
            obj.flipAngle = theflipAngle;
            
        end
        
        function set.TR(obj, theTR)
           
            validateattributes(theTR, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
            obj.TR = theTR;
            
            %Check other DAQ settings:
            if obj.DAQon > obj.TR
                obj.DAQon_start = [];
            end
            
            if obj.DAQoff > obj.TR
                obj.DAQ_end = [];
            end

        end
                
        function set.time(obj, thetime)
            
            validateattributes(thetime, {'numeric'}, {'real', 'vector', 'nonnan', 'nonempty', 'finite', 'nonnegative', 'nondecreasing'});
            
            %Check the size of the time vector.  It must be the same length
            %as the third dimension of the kspacePoints property
            if ~isequal(length(thetime), size(obj.kspacePoints, 3))
                ERRMSG = 'ERROR: the length of the time vector must be the same length as the third dimension of the kspacePoints array';
                error(ERRMSG);
            end
            
            obj.time = thetime;
            
            notify(obj, 'timeVectorSet'); 
        end
        
        function set.sMaps(obj, thesMaps)
            
            validateattributes(thesMaps, {'numeric'}, {'nonnan', 'nonempty', 'finite'});
            
            if ndims(thesMaps) < 3 || ndims(thesMaps) > 4
                ERRMSG = 'ERROR: Incorrect number of dimensions for the sensitivity maps.  ndims of sensitivity maps should be nx by ny by nz by ncoil';
                error(ERRMSG);
            end
            
            obj.sMaps = thesMaps;
            
        end

        function set.relaxivity(obj, theRelaxivity)
            
            validateattributes(theRelaxivity, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
            
            obj.relaxivity = theRelaxivity;
            
        end
        
        function set.samplingTime(obj, theSamplingTime) 
            obj.bandwidth = (theSamplingTime.^-1) * 10^-3*1/2;
        end
        
        function set.TE(obj, ~)

            ERR_MSG = ['Error: TE is a dependent variable.  The value of TE is determined both by the DAQ time and \n'...
                       'the kspace sampling pattern selected.  If the echo time is not correct, check the DAQ time \n'...
                       'settings (obj.setDAQ) and verify that obj.kcenter contains the correct location for the \n'...
                       'center of kspace for your sampling pattern.%s'];
           error(ERR_MSG, '');
        end
        
        function set.DAQon(obj, onTime)

             validateattributes(onTime, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
             
             if onTime >= obj.TR
                 ERR_MSG = 'Error: The selected time is greater than the TR time';
                 error(ERR_MSG);
             end
             
             obj.DAQ_start = onTime;
             
             %Adjust other parameters as necessary:
             if isempty(obj.DAQ_end) && obj.DAQ_start + obj.DAQpercent * obj.TR < obj.TR
                 return;
             elseif isempty(obj.DAQ_end)
                 obj.DAQpercent = (obj.TR-obj.DAQ_start)/obj.TR;
             elseif obj.DAQ_end < obj.DAQ_start 
                 obj.DAQ_end = [];

                 if obj.DAQ_start + obj.DAQpercent * obj.TR > obj.TR
                    obj.DAQpercent = (obj.TR-obj.DAQ_start)/obj.TR;
                 end
             end
        end
        
        function set.DAQoff(obj, offTime)
            
             validateattributes(offTime, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
             
             if offTime >= obj.TR
                 ERR_MSG = 'Error: The selected time is greater than the TR time';
                 error(ERR_MSG);
             end
            
             obj.DAQ_end = offTime;
             
             %Adjust other parameters as necessary:
             if isempty(obj.DAQ_start) && obj.TR - obj.DAQ_end - (obj.DAQpercent * obj.TR) > 0
                 return;
             elseif isempty(obj.DAQ_start)
                 obj.DAQpercent = (obj.DAQ_end)/obj.TR;
             elseif obj.DAQ_start > obj.DAQ_end 
                 obj.DAQ_start = [];

                 if obj.TR - obj.DAQ_end - obj.DAQpercent * obj.TR < 0
                    obj.DAQpercent = (obj.DAQ_end)/obj.TR;
                 end
             end
        end 
        
        function set.imageContrastFunc(obj, fhandle)
            
            if isa(fhandle, 'function_handle')
                obj.imageContrastFunc = fhandle;
            else
                ERRMSG = 'ERROR: imageContrastFunc must be set with a function handle defining the desired Image Contrast Function';
                error(ERRMSG);
            end
            
            
        end
        
        function set.chemSat_On(obj, value)
            
            validateattributes(value, {'numeric'}, {'scalar','nonempty', 'nonnegative' '<=', 1});
            
            if obj.chemSat_On == round(value)
                return;
            end
            
            obj.chemSat_On = round(value);
            
            if obj.chemSat_On == 1
                notify(obj, 'chemSatOn');
            else
                notify(obj, 'chemSatOff');
            end
            
        end
        
        function set.chemSat_flipAngle(obj, theflipAngle)
            
            validateattributes(theflipAngle, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', '<=', 180, '>=', 90});
           
            obj.chemSat_flipAngle = theflipAngle;
            
        end
        
        function set.chemSat_TR(obj, theTR)
           
            validateattributes(theTR, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
            
            if theTR < obj.TR + obj.chemSat_delayedStart + obj.chemSat_delayedEnd
                ERR_MSG = ['The desired Chem Sat TR of %d ms is too short. \n' ...
                            'The chemSat TR must be larger than the sequence TR + any dead time at the begining + any dead time at the end of the sequence:\n\n' ...
                            '    obj.chemSat_TR > obj.TR + obj.chemSat_delayedStart + obj.chemSat_delayedEnd \n' ...
                            '    obj.chemSat_TR > %d + %d + %d'];
                        
                error(ERR_MSG, theTR, obj.TR, obj.chemSat_delayedStart, obj.chemSat_delayedEnd);
            end
            
            obj.chemSat_TR = theTR;
            
            notify(obj, 'chemSatValueChange');
        end
        
        function set.chemSat_delayedStart(obj, theTime)
           
            validateattributes(theTime, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'nonnegative'});
            
            if theTime > obj.chemSat_TR - obj.TR - obj.chemSat_delayedEnd
                ERR_MSG = ['The desired Chem Sat delayed start time of %d ms is too long. \n' ...
                            'The chemSat delayed start time must be smaller than the chemSat TR - the sequence TR - any dead time at the end of the sequence:\n\n' ...
                            '    obj.chemSat_delayedStart < obj.chemSat_TR - obj.TR - obj.chemSat_delayedEnd \n' ...
                            '    obj.chemSat_delayedStart < %d + %d + %d'];
                        
                error(ERR_MSG, theTime, obj.chemSat_TR, obj.TR, obj.chemSat_delayedEnd);
            end
            
            obj.chemSat_delayedStart = theTime;
            
            notify('chemSatValueChange');
        end
        
        function set.chemSat_delayedEnd(obj, theTime)
           
            validateattributes(theTime, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'nonnegative'});
            
            if theTime > obj.chemSat_TR - obj.TR - obj.chemSat_delayedStart
                ERR_MSG = ['The desired Chem Sat delayed end time of %d ms is too long. \n' ...
                            'The chemSat delayed end time must be smaller than the chemSat TR - the sequence TR - any dead time at the begining of the sequence:\n\n' ...
                            '    obj.chemSat_delayedEnd < obj.chemSat_TR - obj.TR - obj.chemSat_delayedStart \n' ...
                            '    obj.chemSat_delayedEnd < %d + %d + %d'];
                        
                error(ERR_MSG, theTime, obj.chemSat_TR, obj.TR, obj.chemSat_delayedStart);
            end
            
            obj.chemSat_delayedEnd = theTime;
            
            notify('chemSatValueChange');
        end
        
        function set.chemSat_T1Fat(obj, fatT1)
           
            validateattributes(fatT1, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'});
            
            obj.chemSat_T1Fat = fatT1;
            
        end
        
        function set.chemSat_TI(obj, theTI)
            
           validateattributes(theTI, {'numeric'}, {'real', 'scalar', 'nonnan', 'nonempty', 'finite', 'positive'}); 
           
            E1 = exp(-obj.chemSat_TR/obj.chemSat_T1Fat);
            E2 = exp(-theTI/obj.chemSat_T1Fat);
            
            a = (-1 + E2)./(E2 * (1-E1));
            
            obj.chemSat_flipAngle = acosd(a/(1+a*E1));
                       
        end
                      
    end
    
    % ----------- Get methods -------------%
    methods
    
        function theSampleTime = get.samplingTime(obj)

            theSampleTime = 1./(2*obj.bandwidth * 10^3);
            
        end
            
        function thenCoil = get.nCoil(obj)
           
            if ~isempty(obj.sMaps)
                thenCoil = size(obj.sMaps, 4);
            else
                thenCoil = 1;
            end
            
        end
        
        function thefinalMatrixSize = get.finalMatrixSize(obj)
           
            if isempty(obj.kmax) || isempty(obj.kmin) || isempty(obj.dk)
                thefinalMatrixSize = [];
                return;
            end
            
            for dim = 1:3
                if min(obj.kmin(dim)) <= 0
                    thefinalMatrixSize(dim) = obj.kmax(dim)./obj.dk(dim) + abs(obj.kmin(dim))./obj.dk(dim);
                else
                    thefinalMatrixSize(dim) = obj.kmax(dim)./obj.dk(dim) - obj.kmin(dim)./obj.dk(dim) + 1;
                end
            end
            
            
            thefinalMatrixSize = floor(thefinalMatrixSize);
            
            %Ensure final matrix size is even:
            thefinalMatrixSize = thefinalMatrixSize + rem(thefinalMatrixSize, 2);
            
        end
        
        function theTE = get.TE(obj)
           
            if isempty(obj.DAQon) 
                theTE = [];
            else
                
                kspace_line = obj.kspacePoints(:,:,1);
                thedist = sum((kspace_line - repmat(obj.kcenter, [length(kspace_line), 1])).^2,2).^(1/2);
                
                [peak, locs] = findpeaks(-1*thedist);
                
                theDAQtime = obj.DAQtime;
                
                theTE = theDAQtime(locs);
                
            end

        end
        
        function theDAQtime = get.DAQtime(obj)

            theDAQtime = linspace(obj.DAQon, obj.DAQoff, size(obj.kspacePoints,1));

        end
        
        function theDAQon = get.DAQon(obj)
            
            if ~isempty(obj.DAQ_start)
                theDAQon = obj.DAQ_start;
            elseif ~isempty(obj.DAQ_end)
                theDAQon = obj.DAQ_end - obj.DAQpercent*obj.TR;
            elseif obj.DAQpercent <= .75
                theDAQon = 0.25 * obj.TR;
            else
                theDAQon = (1-obj.DAQpercent)/2;
            end
            
        end
        
        function theDAQoff = get.DAQoff(obj)
            if ~isempty(obj.DAQ_end)
                theDAQoff = obj.DAQ_end;
            elseif ~isempty(obj.DAQ_start)
                theDAQoff = obj.DAQ_start + obj.DAQpercent*obj.TR;
            elseif obj.DAQpercent <= .75
                theDAQoff =(0.25 + obj.DAQpercent) * obj.TR ;
            else
                theDAQoff = obj.TR - (1-obj.DAQpercent)/2;
            end
        end
        
        function theTI = get.chemSat_TI(obj)
            
           MzTR = (1-exp(-obj.chemSat_TR/obj.chemSat_T1Fat))/(1-cosd(obj.chemSat_flipAngle)*exp(-obj.chemSat_TR/obj.chemSat_T1Fat));
            
           theTI = -log(1/(1-MzTR*cosd(obj.chemSat_flipAngle)))*obj.chemSat_T1Fat;
            
        end
        
    end
    
    % ----------- Other methods -------------%
    
    methods
       
        function validateMRSystem(obj)
            
            %If this MRSystem has already been validated, return
            if obj.MRSystem_Validated
                return;
            end
            
            %Verify that all properties have values:
            theProperties = properties(MRSystem);
            
            %Remove special cases from the list of properties:
            specialProperties = {'sMaps', 'MRSystem_Validated',...
                                 'nCoil', 'samplingTime', 'finalMatrixSize', 'kmax', ...
                                 'kmin', 'kcenter', 'relaxivity', 'userParams', 'TE', 'DAQon', 'DAQoff'};
            thelocs = [];
            
            for counter = 1:length(theProperties)
                
                try 
                    validatestring(theProperties{counter}, specialProperties);
                catch
                    thelocs = [thelocs, counter];
                end
            end
            
            theProperties = theProperties(thelocs);
            
            for counter = 1:length(theProperties)
                if isempty(obj.(theProperties{counter}))
                    ERRMSG = 'ERROR: The property %s has not been set';
                    error(ERRMSG, theProperties{counter});
                end
            end
            
            %Check the special cases as appropriate:
            %   Dependent properties, no check needed:
            %       nCoil
            %       samplingRate
            %       finalMatrix
            %   Optional properties, no check needed:
            %       sMaps
            %       relaxivity
            %       userParams
            %   MRSystem_validated: setAccess protected, no check needed
            

            if isempty(obj.kmax)
                WRNGMSG = 'Warning: no kmax set by user.  Determining kmax based on kspacePoints';
                warning(WRNGMSG);
                obj.autoSetKmax;
            end
            
            if isempty(obj.kmin)
                WRNGMSG = 'Warning: no kmin set by user.  Determining kmin based on kspacePoints';
                warning(WRNGMSG);
                obj.autoSetKmin;
            end
            
            if isempty(obj.kcenter)
                WRNGMSG = 'Warning: no kcenter set by user.  Determining kcenter based on kspacePoints';
                warning(WRNGMSG);
                obj.autoSetKcenter;
            end
            
            if isempty(obj.DAQon)
                WRNGMSG = 'Warning: no DAQ times set by user.  Setting start time equal to 0 and endtime equal to TR';
                warning(WRNGMSG);
                obj.setDAQ(0, obj.TR);
            end

            
            %Verify that the time vector and the 3rd dimension of kspace
            %points are still the same size (i.e. kspacePoints was not
            %reset w/o updateing time)
            if ~isequal(size(obj.kspacePoints, 3), length(obj.time))
                ERRMSG = 'ERROR: the length of the time vector must be the same length as the third dimension of the kspacePoints array';
                error(ERRMSG);
            end
            
            obj.MRSystem_Validated = true;
            
        end
        
        function autoSetKmax(obj)
           
           kxpoints = obj.kspacePoints(:,1,:);
           kxmax = max(kxpoints(:));
           
           kypoints = obj.kspacePoints(:,2,:);
           kymax = max(kypoints(:));
           
           kzpoints = obj.kspacePoints(:,3,:);
           kzmax = max(kzpoints(:));
           
           obj.kmax = [kxmax, kymax, kzmax];
           
        end
        
        function autoSetKmin(obj)
           
           kxpoints = obj.kspacePoints(:,1,:);
           kxmin = min(kxpoints(:));
           
           kypoints = obj.kspacePoints(:,2,:);
           kymin = min(kypoints(:));
           
           kzpoints = obj.kspacePoints(:,3,:);
           kzmin = min(kzpoints(:));
           
           obj.kmin = [kxmin, kymin, kzmin];
           
        end
        
        function autoSetKcenter(obj)
           
            if isempty(obj.kmax)
                obj.autoSetKmax;
            end
            
            if isempty(obj.kmin)
                obj.autoSetKmin;
            end
            
            %Determine if values are given on a standard kspace grid:
            if min(obj.kmin(:)) <=0
                obj.kcenter = [0, 0, 0];
            else
                %If none of the values are below 0, assume values are
                %specified in matlab matrix notation.  Throw a warning to
                %notify the user of the assumption:
                
                obj.kcenter = (obj.kmax - obj.kmin)/2;
                
                WRNG_MSG = 'Warning: No value for center of kspace (obj.kcenter) set, setting  center location to: [%d, %d, %d].';
                warning(WRNG_MSG, obj.kcenter(1), obj.kcenter(2), obj.kcenter(3));
            end
           
        end
        
        function autoSetTime(obj)
           
            if isempty(obj.TR)
                WRNGMSG = 'Warning: value for TR must be set before auto generating a time vector.  value for time not set';
                warning(WRNGMSG);
                return;
            end
            
            obj.time = 0:obj.TR*10^-3:obj.TR*10^-3 * (size(obj.kspacePoints, 3)-1);
            
        end
        
        function setUserParams(obj, fieldname, value)
            
            %Validate the fieldname:
            if ~isvarname(fieldname)
                ERR_MSG = 'Error in setUserParams: the supplied field name is not a valid matlab variable name';
                error(ERR_MSG);
            end
            
            obj.userParams.(fieldname) = value;

        end
        
        function resizesMaps(obj)
           
           origsMaps = obj.sMaps;
           newsMaps = complex(zeros([obj.finalMatrixSize, size(obj.sMaps, 4)]), zeros([obj.finalMatrixSize, size(obj.sMaps, 4)]));
           
           %Prepare new grid:
           iorig = linspace(0, 1, size(origsMaps, 1));
           inew = linspace(0, 1, obj.finalMatrixSize(1));

           jorig = linspace(0, 1, size(origsMaps, 2));
           jnew = linspace(0, 1, obj.finalMatrixSize(2));

           korig = linspace(0, 1, size(origsMaps, 3));
           knew = linspace(0, 1, obj.finalMatrixSize(3));

           [Iorig, Jorig, Korig] = meshgrid(iorig, jorig, korig);
           [Inew, Jnew, Knew] = meshgrid(inew, jnew, knew);
           
           
           %Loop for each coil:
           for counter = 1:size(origsMaps, 4)
                theMag = abs(origsMaps(:,:,:,counter));
                theAngle = angle(origsMaps(:,:,:,counter));
           
                %Make all phase information positive:
                minAngle = min(theAngle(:));
                theAngle = theAngle - minAngle;
           
                %Resize:
                % ** imresize3 only available in 2017 or later **
                %theMag = imresize3(theMag, obj.finalMatrixSize);
                %theAngle = imresize3(theAngle, obj.finalMatrixSize);
                     
                theMag = interp3(Iorig, Jorig, Korig, theMag, Inew, Jnew, Knew, 'cubic');
                theAngle = interp3(Iorig, Jorig, Korig, theAngle, Inew, Jnew, Knew, 'cubic');
           
                %Put phase back to original range: 
                theAngle = theAngle + minAngle;
           
                %Put back into complex array:
                newsMaps(:,:,:,counter) = theMag .* complex(cos(theAngle), sin(theAngle));
           end
           
           %put the updated smaps into the MRsystem object:
           obj.sMaps = newsMaps;

        end
        
        
%         function setDAQ(obj, DAQstart, DAQend, numpts)
%         % this function sets up the DAQ timing vector.  
%         % Inputs:
%         %      DAQstart - a single value indicaiting when the DAQ turns on.
%         %                 This value can also be a vector of start times if
%         %                 the DAQ turns on and off during the aqusition
%         %      DAQend - a single value indicaiting when the DAQ turns off.
%         %                 This value can also be a vector of end times if
%         %                 the DAQ turns on and off during the aqusition
%         %      numpts - Optinal parameter.  If used, it must be the same 
%         %               length as DAQstart and end.  Indicaites how many point during 
%         %               the TR occur between the DAQ start and DAQ end time
%         %               at the corresponding position.  This value must sum
%         %               to the number of points in a single TR.  If not set
%         %               by the user, it will be calculated based on the
%         %               relative amount of time the DAQ is on each time it
%         %               is turned on.
%            
%             %The TR must be set before setting the DAQ:
%             if isempty(obj.TR)
%                 ERR_MSG = ['ERROR: Please set the TR before attempting to set the DAQ'];
%                 error(ERR_MSG);
%             end
%         
%             %Validate the inputs:
%             validateattributes(DAQstart, {'numeric'}, {'real', 'vector', 'nonnan', 'nonempty', 'finite', 'nonnegative', 'nondecreasing'});
%             validateattributes(DAQend, {'numeric'}, {'real', 'vector', 'nonnan', 'nonempty', 'finite', 'nonnegative', 'nondecreasing', 'numel', numel(DAQstart), '<=', obj.TR});
%             
%             %Verify that start and stop times listed do not overlap
%             %eachother:
%             if length(DAQstart) > 1
%                 temparray = zeros([2, length(DAQstart)]);
%                 temparray(1,:) = DAQstart;
%                 temparray(2,:) = DAQend;
%                 temparray = temparray(:);
% 
%                 peaks = findpeaks(temparray);
%                 if ~isempty(peaks)
%                     ERR_MSG = ['ERROR setting DAQ Start/end times: start/end times can not overlap and must be nondecreasing'];
%                     error(ERR_MSG);
%                 end
%             else
%                 if DAQstart > DAQend
%             
%                     ERR_MSG = ['ERROR setting DAQ Start/end times: start/end times can not overlap and must be nondecreasing'];
%                     error(ERR_MSG);
%                 end
%             end
%             
%             obj.DAQ_start = DAQstart;
%             obj.DAQ_end = DAQend;
%             
%             %Validate the numpts variable if it exists:
%             if exist('numpts', 'var')
%                 validateattributes(DAQstart, {'numeric'}, {'real', 'vector', 'nonnan', 'nonempty', 'positive', '<=', size(obj.kspacePoints, 1), 'numel', numel(DAQstart)});
%                 
%                 %ensure that the number of points sum to the number of
%                 %points in a TR:
%                 try
%                     isequal(sum(numpts(:)),size(obj.kspacePoints, 1))
%                 catch
%                     ERR_MSG = ['ERROR: the specified number of points for the DAQ is not equal to the number of points in the TR'];
%                     error(ERR_MSG);
%                 end
%             else
%                 %If numpts was not set by the user, assign it:
%                 segment_time = DAQend - DAQstart;
%                 total_time = sum(segment_time(:));
%                 
%                 numpts = round(size(obj.kspacePoints, 1) * segment_time./total_time);
%                 
%                 %Verify the correct number of points were assigned:
%                 if ~isequal(sum(numpts(:)), size(obj.kspacePoints, 1))
%                    
%                     the_diff = sum(numpts(:)) - size(obj.kspacePoints, 1);
%                     
%                     numpts(1:abs(the_diff)) = numpts(1:abs(the_diff)) + sign(the_diff); 
%                     
%                 end
%             end    
%             obj.DAQ_numpts = numpts;
%         end
        
        function [s0] = SPGRT1Weighting(obj, s0,  t10)
        % SPGRT1Weighting: Apply T1 weighting
        %
        %  Required inputs:
        %     obj - the MRSystem object
        %     s0 - the image space data for weighting to be applied to
        %     t10 - the T1 value of the tissue in the image space data
        %
        % Optional inputs:
        % 
            
            
            %Adjust signal contrast using the model for a perfectly spoiled
            %gradient echo, ignoring T2* effects:
            s0 = s0 * (sind(obj.flipAngle) *...
                (1-exp(-obj.TR/(t10)))/...
                (1-cosd(obj.flipAngle)*exp(-obj.TR/(t10))));
            
%             display(['SPGR T1 Scaling: ' num2str((sind(obj.theMRSystem.flipAngle) *...
%                 (1-exp(-obj.theMRSystem.TR/(t10)))/...
%                 (1-cosd(obj.theMRSystem.flipAngle)*exp(-obj.theMRSystem.TR/(t10)))))]);

        end
        
        function [kdata] = SpecialChemFatSat(obj, kdata,  t10, tissueType)
        % SpecialChemFatSat: Apply a simplified version of the SPECIAL fat 
        % saturation as described in:
        %
        %       Foo TK, Sawyer AM, Faulkner WH, Mills DG. Inversion in the 
        %       steady state: contrast optimization and reduced imaging time 
        %       with fast three-dimensional inversion-recovery-prepared GRE 
        %       pulse sequences. Radiology 1994;191(1):85-90.
        %
        %  For small flip angles, calculating the longitudinal fat recovery 
        %  following the spectrally selective inversion pulse provides a
        %  reasonable approximation of the amount of fat leakage into the
        %  transverse plane as the imaging sequence progresses.
        %
        %  Required inputs:
        %     obj - the MRSystem object
        %     kdata - the kdata for the fat tissue to be supresed
        %     t10 - the T1 value of the tissue in the image space data
        %     tissueType - the tissue type -> verifys this method is only
        %     applied to fat data.
        %
        % Optional inputs:
        % 
            
            %Only apply if obj.chemSat_On has been turned on
            if ~obj.chemSat_On
                return;
            end

            %Only apply this method to fat tissues
            if ~strcmp(tissueType, 'fat')
                return;
            end
            
            %Determine the TR timing in relation to the most recent chemsat
            %pulse
            pulseTime = obj.time - floor(obj.time/(obj.chemSat_TR*10^-3))*obj.chemSat_TR*10^-3;

            %Calculate the longitudinal fat magnetization prior to the
            %sequence RF pulse:
            MzTR = (1-exp(-obj.chemSat_TR / t10))/(1-cosd(obj.chemSat_flipAngle) * exp(-obj.chemSat_TR / t10));            
            Mzt = 1-(1-MzTR*cosd(obj.chemSat_flipAngle))*exp(-pulseTime / (t10*10^-3));
            
            %The weighting for the sequence flip angle was already applied
            %in the image space domain so only apply addiitonal weighting
            %based on longitudinal T1 regrowth following inversion pulse:
            kdata = kdata .* repmat(Mzt, [size(kdata, 1),1]);
            
        end
        
    end
    
    
    % ----------- Protected methods -------------%
    methods (Access = 'protected')
        
        function DAQvalidate(obj, theprop, value)
            
            if strcmp(theprop, 'TR')
                if isempty(obj.DAQ_start)
                    obj.MRSystem_Validated = false;
                elseif obj.DAQ_start > value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested TR is less than the set data aquisition start time. \n  DAQ start = %d \n Requested TR = %d \n';
                    error(ERRMSG, obj.DAQ_start, value);
                end
                
                if isempty(obj.DAQ_end)
                    obj.MRSystem_Validated = false;
                elseif obj.DAQ_end > value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested TR is less than the set data aquisition end time. \n  DAQ end = %d \n Requested TR = %d \n';
                    error(ERRMSG, obj.DAQ_end, value);
                end
            elseif strcmp(theprop, 'DAQ_start')
                
                if isempty(obj.TR)
                    obj.MRSystem_Validated = false;
                elseif obj.TR < value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested data acqusition start time is greater than the set TR time. \n  TR = %d \n Requested DAQ start = %d \n';
                    error(ERRMSG, obj.TR, value);
                end
                
                if isempty(obj.DAQ_end)
                  obj.MRSystem_Validated = false;
                elseif obj.DAQ_end < value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested data acqusition start time is greater than the set data acquisition end time. \n  DAQ end = %d \n Requested DAQ start = %d \n';
                    error(ERRMSG, obj.DAQ_end, value);
                end
                
            elseif strcmp(theprop, 'DAQ_end')
                
                if isempty(obj.TR)
                    obj.MRSystem_Validated = false;
                elseif obj.TR < value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested data acqusition end time is greater than the set TR time. \n  TR = %d \n Requested DAQ end = %d \n';
                    error(ERRMSG, obj.TR, value);
                end
                
                if isempty(obj.DAQ_start)
                  obj.MRSystem_Validated = false;
                elseif obj.DAQ_start > value
                    obj.MRSystem_Validated = false;
                    ERRMSG = 'the Requested data acqusition end time is less than the set data acquisition start time. \n  DAQ start = %d \n Requested DAQ end = %d \n';
                    error(ERRMSG, obj.DAQ_start, value);
                end
            else
                ERRMSG = 'ERROR: Unrecognized property';
                error(ERRMSG);
            end

        end
        
        
    end
    
    
    
    % ----------- Event and Listener Methods -------------%
    methods(Static, Access = 'protected')
        
        %Listener for changes to chem sat and time vector
        function setChemSatTime(obj, eventData)
            
            if strcmp(eventData.EventName, 'timeVectorSet') && ~obj.chemSat_On
                return;
            end
            
            if strcmp(eventData.EventName, 'chemSatOff')
                obj.time = obj.original_time;
                return;
            end
            
            if ~obj.chemSat_On
                return;
            end
            
            if obj.override_timeBackUp
                obj.override_timeBackUp = 0;
                return;
            end
            
            %Save a back up copy of the time vector:
            if ~strcmp(eventData.EventName, 'chemSatValueChange')
                obj.original_time = obj.time;  
            end
            
            
            %Adjust the time vector to accomodate the inversion pulse
            %parameters:
            theTime = obj.original_time;
            
            %Adjust for the first TR seperatly:
            if theTime(1) < obj.chemSat_delayedStart*10^-3
                theTime = theTime + obj.chemSat_delayedStart*10^-3;
            end
            
            %loop through time vector and adjust TR segments as necessary
            %to fit in the interval between inversion pulses:
            dead_time = (obj.chemSat_delayedStart + obj.chemSat_delayedEnd)*10^-3;
            interval_time = (obj.chemSat_TR - obj.chemSat_delayedStart - obj.chemSat_delayedEnd)*10^-3;
            
            start_time = obj.chemSat_delayedStart*10^-3+interval_time;
            while start_time < theTime(end)
                
                %Find the first location that is greater than start time:
                [loc] = find(theTime+obj.TR*10^-3 > start_time+eps, 1, 'first');
                
                if theTime(loc) > start_time + dead_time
                    start_time = start_time + dead_time + interval_time;
                    continue;
                end
                
                theTime(loc:end) = theTime(loc:end) +(round((start_time - theTime(loc))*10^6))*10^-6+dead_time;
                
                start_time = start_time + dead_time + interval_time;
            end

            obj.override_timeBackUp = 1;
            obj.time = theTime;
            
        end
        
    end
    
end