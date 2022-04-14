classdef MRSimulator < handle
    % Classdef for the MRSimulator object 

    properties (SetAccess = 'protected')
    
        %Simulation objects:
        thePhantom = Phantom.empty(1,0);
        theMRSystem = MRSystem.empty(1,0);
        
        %Kdata output:
        kdata = [];
        imageData = []; %Image data befor kspace processing pipeline
        
        %Pipeline function names:
        ImageFuncs = {};
        KspaceFuncs = {};   
        
        %Applied Enhancment Curves:
        enhancementCurves = {};
        
        %Simulation SNR Settings:
        simSNR = [];
    end
    
    properties
                      
        %Simulation settings:
        useGPU = false;
        
        %Output Information:
        outputFilename = [];
        
    end
    
    properties (Dependent)
        %FFT handel
        Simfft
    end

    
    properties (Access = 'protected')
        
        %Cell arrays containing function calls for simulator:
        SimImagePipe = {};
        SimKspacePipe = {};
        
        %Cell array containing information to call user defined DCE kspace
        %functions
        userDCEFuncs = {};
             
        simKspacePoints = []; %Used to re-map kspace points for cartesian sampling  
        
        pipelineFuncs = {}; %List of built in pipeline functions, initialized by constructor
        
    end
    
    properties (Access = 'private')
       
        %Simulation counter variables:
        TissueCounter = 0;
        LayerCounter = 0;
        CoilCounter = 0;
        
        %FFT handel:
        UserDefined_FFT = [];
        
        %Signal Mean for SNR Calculations:
        signalMean = [];
        
    end
    
    
    % ----------- Class Constructor -------------%
    methods
        function obj = MRSimulator(phantom, mrsystem)
            
            % Do nothing with zero inputs...this is required by MATLAB for
            % smooth operation
            if ~nargin
                return
            end
            
            %Set the Phantom and MRSystem properties:
            obj.thePhantom = phantom;
            obj.theMRSystem = mrsystem;
                        
            %Allocate output array:
            obj.kdata = zeros([size(obj.theMRSystem.kspacePoints, 1), length(obj.theMRSystem.time), obj.theMRSystem.nCoil]);
            obj.enhancementCurves = {};     
                        
            
            %% Initialize the list of built in functions:
            objmethods = methods('MRSimulator'); %list of all built in functions
            
            %Built in functions that are not pipeline functions:
            nonpipeline_functions = {'MRSimulator', 'addPipeFunc',...
                                     'addUserDCEModel', 'removePipeFunc',...
                                     'reorderPipeFunc', 'simulate'};
            
            %functions inherited from handle class:                     
            handle_methods = methods('handle');
            
            %Remove the methods inherited from handle from the objmethods list
            for counter = 1:length(handle_methods)
                counter2 = 1;
                while counter2 <= length(objmethods)
                    if strcmp(handle_methods{counter}, objmethods{counter2})
                        objmethods(counter2) = [];
                        continue;
                    else
                        counter2 = counter2 + 1;
                    end
                end
            end
            
            %Remove nonpipeline functions:
            for counter = 1:length(nonpipeline_functions)
                counter2 = 1;
                while counter2 <= length(objmethods)
                    if strcmp(nonpipeline_functions{counter}, objmethods{counter2})
                        objmethods(counter2) = [];
                        continue;
                    else
                        counter2 = counter2 + 1;
                    end
                end
            end
            
            %The functions remaining in objmethods should be the built in
            %pipeline functions:
            obj.pipelineFuncs = objmethods;
            
            %% Add the default methods to the pipeline:
            
            %Apply appropriate signal weighting based on T1, TR and flip
            %angle:
            obj.addPipeFunc('image',  obj.theMRSystem.imageContrastFunc, {'t10'})
                        
            %Apply chemical fat saturation if chemFatSat is turned on:
            if obj.theMRSystem.chemSat_On
                obj.addPipeFunc('kspace',  @(varargin)obj.theMRSystem.SpecialChemFatSat(varargin{:}), {'t10', 'tissueType'});
            end
            
            %Add the signal precession method to the SimkspaceMod pipeline:
            obj.addPipeFunc('kspace',  @obj.signalPrecession);         
            
            %Add the GKM model to the pipeline at the Layer level:
            obj.addPipeFunc('kspace',  @obj.GeneralKineticModel); 

        end
    end
    
    % ----------- set methods -------------% 
    methods
       
        function set.thePhantom(obj, phantom)
           
           %Verify this is a Phantom object:
           if ~isa(phantom, 'Phantom') || isempty(phantom)
               ERRMSG = 'Error: input does not contain a valid Phantom object.  Use Phantom.m to create a valid Phantom object';
               error(ERRMSG);
           end
           
           %Set the Phantom object:
           obj.thePhantom = phantom;
            
        end
        
        function set.theMRSystem(obj, mrsystem)
           
            %Verify this is an MRSystem object:
           if ~isa(mrsystem, 'MRSystem') || isempty(mrsystem)
               ERRMSG = 'Error: input does not contain a valid MRSystem object.  Use MRSystem.m to create a valid MRSystem object';
               error(ERRMSG);
           end
           
           %Check to be sure the MRSystem object has been validated:
           if ~mrsystem.MRSystem_Validated
              ERRMSG = 'Error: the MRSystem object has not been validated.  Run theMRSystem.validateMRSystem';
               error(ERRMSG);
           end
           
           %Set the MRSystem object:
           obj.theMRSystem = mrsystem;
            
        end
        
        function set.Simfft(obj, fhandle)
            
            if isempty(fhandle)
                obj.UserDefined_FFT = [];
            elseif isa(fhandle, 'function_handle')
                obj.UserDefined_FFT = fhandle;
            else
                ERRMSG = 'ERROR: Simfft must be set with a function handle defining the desired Fourier Transform';
                error(ERRMSG);
            end
        end
        
        function set.useGPU(obj, thevalue)
            validateattributes(thevalue, {'logical'}, {'scalar', 'nonempty'});
            obj.useGPU = thevalue;
            
            if ~exist('mex_gpuNUFFT_forw_f.mexw64', 'file')
                WRNG_MSG = ['The executable for the 3D noncartesian fourier transform has not been compiled. \n'...
                'Disable the GPU fuctionality or compiled the function for the GPU NUFFT \n'];
                warning(WRNG_MSG, '');
            end
            
        end
        
        function set.simSNR(obj , thevalue)

            %Check to see if WhiteNoise is in the processing pipeline:
            found_whiteNoise = 0;
            for counter = 1:size(obj.KspaceFuncs, 1)
                if strcmp(obj.KspaceFuncs{counter}, 'WhiteNoise')
                    found_whiteNoise = 1;
                end
            end
            
            if found_whiteNoise == 0
                ERR_MSG = ['Error: Use obj.addWhiteNoise (where obj is the name of your MRSimulator object)\n'...
                           'to add the white noise function to the processing pipeline before setting the SNR'];
                error(ERR_MSG);
            end
            
            validateattributes(thevalue, {'numeric'}, {'real', 'scalar', 'positive', 'nonempty'});
            
            obj.simSNR = thevalue;
            
        end
            
    end
    
    % ----------- get methods -------------%
    methods 
        
        function fhandle=get.Simfft(obj)

            %If there is a user defined FFT, use that one:
            if ~isempty(obj.UserDefined_FFT)
                fhandle = obj.UserDefined_FFT;
                return;
            
            %If this is a Cartesian simulation:
            elseif obj.theMRSystem.isCartesian
                %For Cartesian objects, use the standard Matlab function:
                if obj.theMRSystem.is2D
                    fhandle = @fft2;
                else
                    fhandle = @fftn;
                end
                return;
            
            %If this is a non Cartesian 2D simulation, use the 2D NUFFT:
            elseif obj.theMRSystem.is2D
                fhandle = @iFGG_2d_type2;
                return;
                
            %If the user has selected the GPU option, return that function:    
            elseif obj.useGPU
                fhandle = @gpuNUFFT;

            %Else this is a non Cartesian 3D simulation w/o GPU.  Use the 3D NUFFT
            else
                fhandle = @iFGG_3d_type2;
            end
            
        end
    end
        
    
    
    % ----------- Simulator function -------------% 
    methods
       
        function simulate(obj)
            
            disp('Starting the MRSimulator..... ');
            
            %Run some pre-simulation checks:
            
            %Verify that the VIF is longer than the requested simulation
            %time:
            if obj.thePhantom.VIF(end, 1) < obj.theMRSystem.time(end)
                err_msg=['Error: The requested simulation time defined by the time vector in the MR System object is longer than the input VIF.  Unable to simulate.'];
                error(error_msg);
            end
            

            %Resize the phantom to the appropriate size:
            if ~isequal(obj.thePhantom.xyzDims, obj.theMRSystem.finalMatrixSize)
                disp('Resizing the Phantom object to match the final matrix size given in the MRSystem object ..... ');
                obj.thePhantom.xyzDims = obj.theMRSystem.finalMatrixSize;
                             
                
                
                disp('...Resizing Completed on the Phantom object ..... ');
            end
            
            if ~isempty(obj.theMRSystem.sMaps) 
                if ~isequal(size(obj.theMRSystem.sMaps), obj.theMRSystem.finalMatrixSize)
                    disp('Resizing the Sensitivity Maps to match the final matrix size given in the MRSystem object ..... ');
                    
                    obj.theMRSystem.resizesMaps;
                    disp('...Resizing Completed on the sensitivity maps object ..... ');
                end
            end
            
            %Clear any data already existing in obj.kdata:
            obj.kdata = zeros([size(obj.theMRSystem.kspacePoints, 1), size(obj.theMRSystem.kspacePoints, 3), obj.theMRSystem.nCoil]);
            obj.imageData = zeros([obj.thePhantom.xyzDims(2), obj.thePhantom.xyzDims(1), obj.thePhantom.xyzDims(3)]);
            
            %Initialize the output array holding the percent enhancement
            %for each tissue type:
            obj.enhancementCurves = cell([1, size(obj.thePhantom.theTissues, 2)]);
            for counter = 1:size(obj.thePhantom.theTissues, 2)
                obj.enhancementCurves(1,counter) = {zeros([length(obj.theMRSystem.time), length(obj.thePhantom.theTissues{counter}.tissueLayers)])};
            end
            
            %Perform any required pre-processing before loops:
            
            %If this is a Cartesian object using fftn, speed operations by transforming kspace points to represent matrix locations 
            if obj.theMRSystem.isCartesian && (strcmp(func2str(obj.Simfft), 'fftn') || strcmp(func2str(obj.Simfft), 'fft2'))
               
               if obj.theMRSystem.is2D
                   ERRMSG = 'ERROR: Not set up for 2D simulations yet';
                    error(ERRMSG);
               end
                
                
               disp('Mapping Cartesian kspace points to desired matrix index locations...');
               obj.simKspacePoints = round(obj.theMRSystem.kspacePoints...
                                     - repmat(obj.theMRSystem.kmin, [size(obj.theMRSystem.kspacePoints, 1), 1, size(obj.theMRSystem.kspacePoints, 3)])...
                                     ./repmat(obj.theMRSystem.dk, [size(obj.theMRSystem.kspacePoints, 1), 1, size(obj.theMRSystem.kspacePoints, 3)])) + 1;
               
               %Calculate new kmax values:                          
               kmax = round(obj.theMRSystem.kmax - obj.theMRSystem.kmin)./obj.theMRSystem.dk + 1;
               
               %To further speed operations by avoiding use of fftshift,
               %calculate the index of the desired points before fftshift
               %is applied:
               obj.simKspacePoints = obj.fftshift_points(obj.theMRSystem.kspacePoints, kmax);
               
               disp('... Cartesian Mapping completed'); 
            end
            
            
            
            %Loop by coil:
            for nCoil = 1:obj.theMRSystem.nCoil
                obj.CoilCounter = nCoil;
                display(['Looping by coil.  Coil number: ' num2str(nCoil) ' of ' num2str(obj.theMRSystem.nCoil)]);
                
                %Loop by Tissue:
                for nTissue = 1:length(obj.thePhantom.theTissues)
                    obj.TissueCounter = nTissue;
                    display(['... Looping by Tissue.  Tissue number: ' num2str(nTissue) ' of ' num2str(length(obj.thePhantom.theTissues))]);
                    display(['......The Tissue Name is: ' obj.thePhantom.theTissues{nTissue}.tissueName]); 
                    
                    %Obtain the pre-contrast images for the tissue:
                    s0 = obj.thePhantom.theTissues{nTissue}.s0;
                    
                    %Apply the noise mask:
                    s0 = s0.*obj.thePhantom.theTissues{nTissue}.noiseMask;
                    
                    %Apply any image based operations:
                    %Loop through each layer of the tissue:
                    if ~isempty(obj.SimImagePipe)
                        news0 = zeros(size(s0));
                        for nLayer = 1:length(obj.thePhantom.theTissues{nTissue}.tissueLayers)
                            obj.LayerCounter = nLayer;
                            display(['......Performing image based operations on layer ' num2str(nLayer) ' of ' num2str(length(obj.thePhantom.theTissues{nTissue}.tissueLayers))]);
                            layers0 = s0 .* obj.thePhantom.theTissues{nTissue}.getLayerMask(nLayer);
                            tissueProps = obj.grabTissueParams(nTissue, nLayer);
                            
                            %Loop through image pipeline:
                            for pipelineCounter = 1:size(obj.SimImagePipe, 1)

                                %apply pipeline operations:
                                fhandle = obj.SimImagePipe{pipelineCounter, 1};
                                params = obj.SimImagePipe{pipelineCounter, 2};
                                %obtain actual values for params:
                                for param_counter = 1:length(params)
                                    try
                                        params(param_counter) = {eval(params{param_counter})};
                                    catch
                                        params(param_counter) = {[]};
                                    end
                                end
                                
                                layers0 = fhandle(layers0, params{1:length(params)});
                            end
                            news0 = news0 + layers0;
                        end
                        obj.LayerCounter = 0;
                        s0 = news0;
                        clear news0;
                        clear layers0;
                    end
                            
                    obj.imageData = obj.imageData + s0;

                    %Apply the coil map:
                    if ~isempty(obj.theMRSystem.sMaps)
                        s0 = s0.*obj.theMRSystem.sMaps(:,:,:,nCoil);
                    end  
                    
                    
                    %Loop through each layer of the tissue:
                    if isempty(obj.SimKspacePipe)
                        %Perform Fourier Transform on S0:
                        disp('......Performing FT of S0');
                        obj.kdata(:,:,nCoil) = obj.kdata(:,:,nCoil) + obj.obtainKspace(s0);
                        disp('......... FT of S0 complete');
                    else
                        %Loop through each layer, perform kspace based
                        %processes on a per layer basis:
                        for nLayer = 1:length(obj.thePhantom.theTissues{nTissue}.tissueLayers)
                            obj.LayerCounter = nLayer;
                            display(['......Performing FT of layer ' num2str(nLayer) ' of ' num2str(length(obj.thePhantom.theTissues{nTissue}.tissueLayers))]);
                            thekdata = obj.obtainKspace(s0.*obj.thePhantom.theTissues{nTissue}.getLayerMask(nLayer));
                            disp(['......... FT of layer ' num2str(nLayer) ' complete']);

                            %obtain tissue properties for the current
                            %layer:
                            tissueProps = obj.grabTissueParams(nTissue, nLayer);
                            
                            %Loop through operations pipeline:
                            for pipelineCounter = 1:size(obj.SimKspacePipe,1)

                                %apply pipeline operations:
                                fhandle = obj.SimKspacePipe{pipelineCounter, 1};
                                params = obj.SimKspacePipe{pipelineCounter, 2};
                                %obtain actual values for params:
                                for param_counter = 1:length(params)
                                    try
                                        params(param_counter) = {eval(params{param_counter})};
                                    catch
                                        params(param_counter) = {[]};
                                    end
                                end
                                
                                %apply pipeline processing:
                                [thekdata] = fhandle(thekdata, params{1:length(params)});
                                
                                %Grab the enhancement information returned
                                %by the General 
                            end
                            
                            obj.kdata(:,:,nCoil) = obj.kdata(:,:,nCoil) + thekdata;
                            
                        end
                        obj.LayerCounter = 0;
                        
                    end
                    
                end %End tissue looop
                obj.TissueCounter = 0;
                
            end %end coil loop
            obj.CoilCounter = 0;
            
        end
        
    end
    
    
    % ----------- Functions for simulation Pipeline -------------% 
    methods
        
        %Image pipeline functions:

        
        %Kspace pipeline functions:
        function [thekdata, output] = signalPrecession(obj, thekdata, nCoil, nTissue, nLayer, pipelineCounter, varargin)
            
            output = []; %no output to return from this function

            %If There is no precession to apply, return:
            if isequal(obj.thePhantom.theTissues{nTissue}.signalPrecession.ppm, 0) && isequal(obj.thePhantom.theTissues{nTissue}.signalPrecession.relAmp, 1)
                return;
            end
            
            disp('......Applying signal precession');
            
            %Set up the time vector:
            time = obj.theMRSystem.DAQtime * 10^-3;
            time = time(:);
            
            %Calculate the frequency precession:
            gamma = 42.58 * 10^6; %Hz/tesla
            ppm = obj.thePhantom.theTissues{nTissue}.signalPrecession.ppm; %ppm shift relative to water
            relAmp = obj.thePhantom.theTissues{nTissue}.signalPrecession.relAmp; %relative amplitude of each component
            fieldStrength = obj.theMRSystem.fieldStrength; %field strength in T
            
            %Get the frequency shift in Hz:
            freq_shift = fieldStrength * gamma * ppm * 10^-6;
            
            %Calculate the changes in magnitude / phase over the course of
            %the TR:
            scaling_vector = repmat(relAmp, [length(time), 1]) .*...
                             complex(cos(repmat(freq_shift, [length(time), 1]) .* repmat(time, [1, length(relAmp)])*2*pi()),...
                                     sin(repmat(freq_shift, [length(time), 1]) .* repmat(time, [1, length(relAmp)])*2*pi()));
            scaling_vector = sum(scaling_vector, 2);
            
            %Apply the signal precession to the kdata:
            anglekdata = angle(thekdata);
            anglekdata = anglekdata + repmat(unwrap(angle(scaling_vector)), [1, size(thekdata, 2)]);
            
            thekdata = abs(thekdata).* repmat(abs(scaling_vector), [1, size(thekdata, 2)]) .* complex(cos(anglekdata), sin(anglekdata));
            clear anglekdata;

            disp('.........Signal precession complete');
        end
        
        function [thekdata, PercentEnhancement] = GeneralKineticModel(obj, thekdata, nCoil, nTissue, nLayer, pipelineCounter, varargin)
                        
            %Verify GKM function exists:
            if ~exist('gkm_ve')
                ERRMSG = 'ERROR: The function gkm_ve does not exist on the search path';
                error(ERRMSG);
            end
            
            %Verify other necessary inputs have been set:
            if isempty(obj.thePhantom.VIF)
                ERRMSG = ['ERROR: the vascular input function  (VIF) of the Phantom \n'...
                          'object is unset. VIF is required for pharmacokinetic modeling \n'...
                          'using the general kinetic model function.  Either set the VIF or \n'...
                          'remove the GeneralKineticModel function from your processing \n'...
                          'pipeline. %s'];
                error(ERRMSG, '');
            end
            
            if isempty(obj.theMRSystem.relaxivity)
                ERRMSG = ['ERROR: the contrast relaxivity value has not been set.  \n'...
                          'Please set this value the MRSystem object or remove the \n'...
                          'GeneralKineticModel function from your processing pipeline. %s'];
                error(ERRMSG, '');
            end
            
            if ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 'ktrans') || ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 've') || ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 'vp')
                WRNGMSG = ['Warning: variables ktrans, ve, and/or vp required by \n'...
                           'GeneralKineticModel function are unset.  Model not applied. %s'];
                warning(WRNGMSG);
                return;
            end
                           
            if isempty(obj.thePhantom.theTissues{nTissue}.modelParams.ktrans{nLayer}) ||...
               isempty(obj.thePhantom.theTissues{nTissue}.modelParams.ve{nLayer})||...
               isempty(obj.thePhantom.theTissues{nTissue}.modelParams.vp{nLayer})
                %No model parameters for this layer, do not process
                return;
            end
            
            disp('......Applying Pharmacokinetic Modeling');
            
            %split the vif apart into time and concentration vectors:
            vif = obj.thePhantom.VIF;
            vif_time = vif(:,1);
            vif = vif(:,2);
            
            %Convert vif from M of Gd to mM
            %vif = 1000*vif; %convert from M of Gd to mM
            
            %Get the model parameters:
            ktrans = obj.thePhantom.theTissues{nTissue}.modelParams.ktrans{nLayer};
            ve = obj.thePhantom.theTissues{nTissue}.modelParams.ve{nLayer};
            vp = obj.thePhantom.theTissues{nTissue}.modelParams.vp{nLayer};
            t10 = obj.thePhantom.theTissues{nTissue}.t10(nLayer);
            
            %Get the simulated curve:
            thecurve = gkm_ve([ktrans, ve, vp], vif_time, vif, 0.45);
            
            %interpolate the curve to match the simulation time curve:
            thecurve = interp1(vif_time, thecurve, obj.theMRSystem.time);
            
            %Convert from Gd concentration to relative signal change:
            T1 = (1/t10  + obj.theMRSystem.relaxivity*thecurve).^-1;
            PercentEnhancement = fspgr_full(1, T1, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)./fspgr_full(1, t10, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)-1;

            %Record the percent enhancement information:
            obj.enhancementCurves{1, nTissue}(:,nLayer) = PercentEnhancement;
            
            %Apply the curve to the kdata:
            thekdata = thekdata + thekdata .* repmat(PercentEnhancement, [size(thekdata, 1), 1]);
            
        end
        
        function [thekdata] = WhiteNoise(obj, thekdata, nCoil, nTissue, nLayer, pipelineCounter)
        %This function applies white gaussian noise to the kdata corresponding to the desired SNR level
        %
        %  Required inputs:
        %  Most of these inputs are used to locate the approprite, tissue
        %  specific data in the phantom.
        %     obj - the MRSim object
        %     thekdata - the kspace data noise is to be added to
        %     nCoil - the current coil number, used to locate relevent data
        %               in the phantom
        %     nTissue - the current tissue number, used to locate relevent 
        %                tissue specific data in the phantom object           
        %     nLayer - the current layer number, used to locate relevent
        %                tissue specific data in the phantom object.
        %
        % To add the function to the processing pipeline use the following
        % option:
        %
        %       obj.addWhiteNoise
        %
        % Where obj is the name of your MRSimulator object
            
            %This is really a global function, not a tissue or layer
            %specific function so only apply it once:
            if ~isequal(nTissue, length(obj.thePhantom.theTissues)) || ~isequal(nLayer, 1)
                return;
            end
            
            disp('......Adding White Noise');
            
            %Add the noise to the kdata
            variance = abs(obj.signalMean)/obj.simSNR;
            
            %Scale the variance by the matrix size:
            variance = variance * obj.theMRSystem.finalMatrixSize(1)*obj.theMRSystem.finalMatrixSize(2)*obj.theMRSystem.finalMatrixSize(3);
            
            %Generate normally distributed white gaussian noise:
            noise_matrix = complex(randn(size(thekdata)), randn(size(thekdata)))*1/sqrt(2)*sqrt(variance);
            thekdata = thekdata + noise_matrix;

        end
        
    end
    
    
    % ----------- Functions for controling the pipeline -------------% 
    methods
        
        function addPipeFunc(obj, pipeID, fhandle, varargin)
        %addPipeFunc(obj, pipeID, fhandle, parameterNames, optionalInputs) 
        %adds a function to the processing pipeline.
        %
        %    addPipeFunc(obj, PIPEID, FHANDLE, parameterNames) adds the 
        %    function specified by the function handle FHANDLE to the 
        %    processing pipeline indicaited by PIPEID.  Options for PIPEID 
        %    are 'image' for an image space processing function and 'kspace' 
        %    for a kspace based processing function.  ParameterNames is a 
        %    cell array contaitaining a list of the parameters required by
        %    the function specified by fhandle.  The first parameter passed
        %    to the function FHANDLE is the image space or kspace data as 
        %    indicaited by PIPEID and does no need to be specified in
        %    parameterNames.  ParameterNames is only required for user 
        %    defined functions, pipeline functions that 
        %    are part of the MRSimulator class will set this variable 
        %    automatically.  
        %
        %    A list of possible parameters can by obtained by calling the function 
        %    'displayPossibleInputs' of the MRSimulator object.  When the
        %    function specified by FHANDLE is called during simulation
        %    processing, the parameters corresponding to the names in
        %    parameterNames will be passed in order to the function 
        %    specified by FHANDLE.
        %
        %    addPipeFunc(..., optionalInputs) the optional input parameter
        %    'optionalInputs' allows the user to specify settings for the
        %    function specified by FHANDLE.  For a lsit of optional input
        %    parameters for built in pipeline functions, type 
        %    help MRSimulator.FUNCTIONNAME 
            
            %validate the pipeID input:
            valid_pipelines = {'image', 'kspace'};
            pipeID = lower(pipeID);
            
            validatestring(pipeID, valid_pipelines);
            
            if strcmp(pipeID, 'image')
                pipename = 'SimImagePipe';
                listname = 'ImageFuncs';
            else
                pipename = 'SimKspacePipe';
                listname = 'KspaceFuncs';
            end
            
            %If it is a valid function handle, add it to the pipeline:
            if isa(fhandle, 'function_handle')
                theLocation = size(obj.(pipename), 1)+1;
                obj.(pipename)(theLocation, 1) = {fhandle};
            else
                ERRMSG = ['ERROR: processing pipelines must be set with a function handle defining the desired function to modulate the raw kspace or image data'];
                error(ERRMSG);
            end
            
            %Set parameterNames and optionalInputs as appropriate
            if (length(varargin) >= 1)
                if ~isempty(varargin{1})
                    parameterNames = varargin{1};
                else
                    parameterNames = {};
                end
            else
                parameterNames = {};
            end
            if (length(varargin) >= 2)
                if ~isempty(varargin{2})
                    optionalInputs = varargin{2};
                else
                    optionalInputs = {};
                end
            else
                optionalInputs = {};
            end
            
            %Save the name of the function:
            fname = func2str(fhandle);
            if strcmp(fname(1:15), '@(varargin)obj.')
                loc = strfind(fname, '(varargin{:})');
                tempname = '';
                try
                    tempname = fname(16:loc-1);
                catch
                    tempname = fname;
                end
                
                fname = tempname;
            end
            
            %Record the function name in the list of names:
            obj.(listname)(theLocation, 1) = {fname};
            
            %Adjust the length of the obj.userDCEFuncs array:
            if strcmp(pipename, 'SimKspacePipe')
                obj.userDCEFuncs(theLocation) = {[]};
            end
              
            %Determine if the input function is a built in function.  If so, only counter variables are needed (varargin may contain additional options):
            for counter = 1:length(obj.pipelineFuncs)
                if strcmp(func2str(fhandle), ['@(varargin)obj.' obj.pipelineFuncs{counter} '(varargin{:})']) 
                    %Save data for call in the paramCall cell array:
                    paramCall = {'nCoil', 'nTissue', 'nLayer', 'pipelineCounter'};
                    obj.(pipename)(theLocation, 2) = {paramCall};
                    
                    %Save any optional user defined inputs for later use:
                    if exist('optionalInputs', 'var') && ~isempty(optionalInputs)
                        obj.(pipename)(theLocation, 3) = {optionalInputs};
                    end
                    
                    return;
                end
            end
            
            %The function is a user defined function.  Sort out and
            %identify the desired input parameters:
            if ~exist('parameterNames', 'var')
                paramCall = {};
            else
                paramCall = searchForInputs(obj, parameterNames);
            end
            
            %Put the data into the pipeline:
            obj.(pipename)(theLocation, 2) = {paramCall};
            
            %Include any user defined parameters:
            obj.(pipename)(theLocation, 3) = {optionalInputs};
        end
        
        function removePipeFunc(obj, pipeID, funcID)
        % This function removes a function from the processing pipline.  
        % Inputs:
        %   pipeID - string to indicaite which processing pipeline the
        %   function is to be removed from.  valid options are 'image' and
        %   'kspace'
        %
        %   funcID - variable to ID the function to be removed.  This must
        %   either be a number indicaiting the functions location in the
        %   pipline or a string identifying it by name.
        
            %validate the pipeID input:
            valid_pipelines = {'image', 'kspace'};
            pipeID = lower(pipeID);
            
            validatestring(pipeID, valid_pipelines);
            
            if strcmp(pipeID, 'image')
                thenames = obj.ImageFuncs;
                pipename = 'SimImagePipe';
                listname = 'ImageFuncs';
            else
                thenames = obj.KspaceFuncs;
                pipename = 'SimKspacePipe';
                listname = 'KspaceFuncs';
            end
        
            %Validate the funcID input:
            if isnumeric(funcID)
                try
                    validateattributes(funcID, {'numeric'}, {'scalar', 'positive', '<=', length(thenames)});
                catch
                    %Throw a more descriptive error:
                    ERR_MSG = ['Error attempting to remove function from proccessing pipeline: function not found. \n' ...
                                'Provide a valid name or location number of the function to be removed. %s'];
                            
                    error(ERR_MSG, '');
                end
                pipeNum = floor(funcID);
            else
                
                pipeNum = [];
                %Identify the location of the named layer:
                for counter = 1:length(thenames)
                    
                    if strcmp(funcID, thenames{counter})
                        pipeNum = counter;
                        break;
                    end
                end
                
                %If a matching function was not found by examining the list
                %of function names, look at the functions themsleves in
                %case the full name was used:
                if isempty(pipeNum)
                    
                    for counter = 1:length(obj.(pipename))
                        if strcmp(funcID, obj.(pipename){counter})
                            pipeNum = counter;
                            break;
                        end 
                    end
                end
                
                if isempty(pipeNum)
                    ERR_MSG = ['Error attempting to remove function from proccessing pipeline: function not found. \n' ...
                                'Provide a valid name or location number of the function to be removed. %s'];
                            
                    error(ERR_MSG, '');
                end
            end
            
            %Remove the desired function from the desired pipleine:
            obj.(pipename)(pipeNum) = [];
            obj.(listname)(pipeNum) = [];
            
            %Adjust the array containing user DCE functions as appropriate:
            if strcmp(pipename, 'SimKspacePipe')
                obj.userDCEFuncs(pipeNum) = [];
            end
            
        end
        
        function reorderPipeFunc(obj, pipeID, funcID, newLoc)
        % This function moves the function indicaited by funcID to the 
        % location in the pipeline indicaited by newLoc    
        % Inputs:
        %   pipeID - string to indicaite which processing pipeline the
        %   function is to be removed from.  valid options are 'image' and
        %   'kspace'
        %
        %   funcID - variable to ID the function to be removed.  This must
        %   either be a number indicaiting the functions location in the
        %   pipline or a string identifying it by name.
        %
        %   newLoc - a number indicaiting the desired location in the
        %   pipeline for this function.  
        
            %validate the pipeID input:
            valid_piplines = {'image', 'kspace'};
            pipeID = lower(pipeID);
            
            validatestring(pipeID, valid_pipelines);
            
            if strcmp(pipeID, 'image')
                thenames = obj.ImageFuncs;
                pipename = 'SimImagePipe';
                listname = 'ImageFuncs';
            else
                thenames = obj.KspaceFuncs;
                pipename = 'SimKspacePipe';
                listname = 'KspaceFuncs';
            end
            
            %Validate the newLoc input
            validateattributes(newLoc, {'numeric'}, {'scalar', 'positive', '<=', length(thenames)});
        
            %Validate the funcID input:
            if isnumeric(funcID)
                try
                    validateattributes(funcID, {'numeric'}, {'scalar', 'positive', '<=', length(thenames)});
                catch
                    %Throw a more descriptive error:
                    ERR_MSG = ['Error attempting to reorder functions in the processing pipeline: function not found. \n' ...
                                'Provide a valid name or location number of the function to be moved. %s'];
                            
                    error(ERR_MSG, '');
                end
                pipeNum = floor(funcID);
            else
                
                pipeNum = [];
                %Identify the location of the named layer:
                for counter = 1:length(thenames)
                    
                    if strcmp(funcID, thenames{counter})
                        pipeNum = counter;
                        break;
                    end
                end
                
                %If a matching function was not found by examining the list
                %of function names, look at the functions themsleves in
                %case the full name was used:
                if isempty(pipeNum)
                    
                    for counter = 1:length(obj.(pipename))
                        if strcmp(funcID, obj.(pipename){counter})
                            pipeNum = counter;
                            break;
                        end 
                    end
                end
                
                if isempty(pipeNum)
                    ERR_MSG = ['Error attempting to reorder functions in the processing pipeline: function not found. \n' ...
                                'Provide a valid name or location number of the function to be moved. %s'];
                            
                    error(ERR_MSG, '');
                end
            end
        
            %Adjust array tracking user DCE data as appropriate:
            if strcmp(pipename, 'SimKspacePipe')
                obj.reorderUserData(pipeNum, newLoc);
            end
            
            %Remove the function handle from its original location:
            temphandle = obj.(pipename)(pipeNum);
            tempname = obj.(listname)(pipeNum);
            
            obj.(pipename)(pipeNum) =[];
            obj.(listname)(pipeNum) = [];
            
            if newLoc > pipeNum
                newLoc = newLoc - 1;
            end
            
            %put it in its new location:
            if newLoc == 1
                obj.(pipename) = [temphandle, obj.(pipename)];
                obj.(listname) = [tempname, obj.(listname)];
            elseif newLoc > length(obj.(pipename))
                obj.(pipename)(newLoc) = temphandle;
                obj.(listname)(newLoc) = tempname;
            else
                obj.(pipename) = [obj.(pipename)(1:newLoc - 1), temphandle, obj.(pipename)(newLoc:length(obj.(pipename)))];
                obj.(listname) = [obj.(listname)(1:newLoc - 1), tempname, obj.(listname)(newLoc:length(obj.(listname)))];
            end
            
        end
        
        function addUserDCEModel(obj, fhandle, varargin)
            %This function will insert a user specified dynamic enhancement
            %model into the kspace processing pipeline. The function specified by fhandle
            %must take as an input an MRSystem object and a structure
            %containing the specific tissue properties for the tissue in
            %question.  The function must return a time vector and a vector
            %describing the percent enhancement (ranging from 0 to inf,
            %where 1 represents a 100% enhancement over baseline or the
            %signal has doubled in intensity).
            %
            %No error checking is performed on the input model parameters,
            %the function must handle its own error checking to ensure that
            %all necessary variables have been set and that all the results
            %are valid options.
            %
            %Optional: define the names of the variables that are required
            %by the function.  
            
            if ~isa(fhandle, 'function_handle')
                ERRMSG = 'ERROR: applyUserDCEModel must be set with a function handle defining the desired function to modulate the raw kspace data';
                error(ERRMSG);
            end
            
            %Identify the inputs required by the user function:
            if ~exist('varargin', 'var')
                paramCall = {};
            else
                paramCall = searchForInputs(obj, varargin);
            end
            
            %Add call to the pipeline:
            theLocation = length(obj.SimKspacePipe)+1;
            obj.SimKspacePipe(theLocation, 1) = {@obj.applyUserDCEModel};
            obj.SimKspacePipe(theLocation, 2) = {'numTissue', 'nlayer', 'obj.pipelineCounter'};
            
            %Save data for call in the userDCEFuncs cell array:
            obj.userDCEFuncs(length(obj.SimKspacePipe)) = {fhandle, paramCall};
        end
        
        function displayPossibleInputs(obj)
            
            [inputOptions] = searchForInputs(obj, {'-h'});
           
             sectionNames = inputOptions{1};
             parameters = inputOptions(2:size(inputOptions, 1), 1);
             
             if isequal(size(sectionNames, 2), size(parameters, 1))
                 
                 for counter = 1:length(sectionNames) 
                     disp([sectionNames{counter} ':']);
                     
                     %Get the variables for that section:
                     theParams = parameters{counter, 1};
                     for counter2 = 1:length(theParams)
                         display(['     ' theParams{counter2}]);
                     end
                 end
                 
             else
                 %Display the variables for each section:
                 for counter = 1:size(parameters, 1)
                    theParams = parameters{counter, 1};
                     for counter2 = 1:length(theParams)
                         display(['     ' theParams{counter2}]);
                     end
                 end
             end
        end
                
    end
    
    % ----------- Other Methods -------------% 
    methods
       
        function [enhancement_curves] = displayGeneralKineticModel(obj)
            
            %split the vif apart into time and concentration vectors:
            vif = obj.thePhantom.VIF;
            vif_time = vif(:,1);
            vif = vif(:,2);
            
            %Convert vif from M of Gd to mM
            %vif = 1000*vif; %convert from M of Gd to mM
            
            %Set up variable to hold the curves:
            enhancement_curves = cell(size(obj.thePhantom.theTissues));
            
            for nTissue = 1:length(enhancement_curves)
                
                enhancement_curves(nTissue)={zeros([length(obj.thePhantom.theTissues{nTissue}.tissueLayers), length(vif_time)])};
                
                for nLayer = 1:length(obj.thePhantom.theTissues{nTissue}.tissueLayers)

                    if ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 'ktrans') ||...
                       ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 've') ||...
                       ~isfield(obj.thePhantom.theTissues{nTissue}.modelParams, 'vp')
                        continue;
                    end
                    
                    %Get the model parameters:
                    ktrans = obj.thePhantom.theTissues{nTissue}.modelParams.ktrans{nLayer};
                    ve = obj.thePhantom.theTissues{nTissue}.modelParams.ve{nLayer};
                    vp = obj.thePhantom.theTissues{nTissue}.modelParams.vp{nLayer};
                    t10 = obj.thePhantom.theTissues{nTissue}.t10(nLayer);
                    
                    if isempty(ktrans) || isempty(ve) || isempty(vp)
                        continue;
                    end
            
                    %Get the simulated curve:
                    thecurve = gkm_ve([ktrans, ve, vp], vif_time, vif, 0.45);
                        
                    %Convert from Gd concentration to relative signal change:
                    T1 = (1/t10  + obj.theMRSystem.relaxivity*thecurve).^-1;
                    PercentEnhancement = fspgr_full(1, T1, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)./fspgr_full(1, t10, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)-1;

                    enhancement_curves{nTissue}(nLayer, :) = PercentEnhancement;
                end
            end
            
            vif_curve = gkm_ve([1, 0, 0.55], vif_time, vif, 0.45);
            T1 = (1/1440  + obj.theMRSystem.relaxivity*vif_curve).^-1;
            VIFPercentEnhancement = fspgr_full(1, T1, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)./fspgr_full(1, 1440, obj.theMRSystem.flipAngle, obj.theMRSystem.TR*10^-3)-1;

            
            %Plot the Results
            for nTissue = 1:length(enhancement_curves)
                
                figure; 
                hold on;
                plot(vif_time, VIFPercentEnhancement)
                title('General Kinetic Model Percent Enhancement Curves');
                
                legend_title = cell([length(obj.thePhantom.theTissues{nTissue}.tissueLayers) + 1,1]);
                legend_title(1) = {'VIF'};
                color_string = ['g', 'r', 'c', 'm', 'y', 'k', 'b'];
                color_string = repmat(color_string, [1, ceil(length(obj.thePhantom.theTissues{nTissue}.tissueLayers)/length(color_string))]);
                
                for nLayer = 1:length(obj.thePhantom.theTissues{nTissue}.tissueLayers)
                    
                    plot(vif_time, enhancement_curves{nTissue}(nLayer, :), color_string(nLayer));
                    
                    legend_title(nLayer+ 1) = {obj.thePhantom.theTissues{nTissue}.layerNames{nLayer}};
                    
                end
                
                hold off;
                legend(legend_title{1:end});
                
            end
                    

        end
        
        function addWhiteNoise(obj, varargin)
        % function addWhiteNoise(obj, varargin)
        % This function applies adds the white noise function to the
        % processing pipeline.
        %
        % optional inputs allow the user to override the values stored in
        % the phantom object to use the function.  Define these values in
        % the following manner:
        %
        %      addWhiteNoise(..., 'property', value)     
        %
        % Optional inputs:
        %   SNR - the signal to noise ratio.  Default is 20.
        %   signalMean - the mean signal of the image space object.  If not
        %       specified, the algorithm will determine it from the image
        %       space data.
        %
        
            %Check to see if the WhiteNoise function has been added to the
            %processing pipeline.  If not, add it:
            found_whiteNoise = 0;
            for counter = 1:size(obj.KspaceFuncs, 1)
                if strcmp(obj.KspaceFuncs{counter}, 'WhiteNoise')
                    found_whiteNoise = 1;
                end
            end
            
            %If the WhiteNoise function is not there, add it:
            if found_whiteNoise == 0
                obj.addPipeFunc('kspace',  @obj.WhiteNoise);
            end
            
            %Initialize values:
            obj.simSNR = 20;
            obj.signalMean = [];
                
            %Sort out any optional input data:           
            for counter = 1:length(varargin)
                if ~ischar(varargin{counter})
                    continue;
                end
                
                if strcmp(varargin{counter}, 'SNR')
                    set_value = true;
                    try
                        validateattributes(varargin{counter + 1}, {'numeric'}, {'scalar', 'positive', 'nonempty', 'finite'});
                    catch
                        WRN_MSG = 'User Specified SNR is not valid, using default value';
                        warning(WRN_MSG);
                        set_value = false;
                    end
                    
                    if set_value
                        obj.simSNR = varargin{counter+ 1};
                    end
                end
                
                if strcmp(varargin{counter}, 'signalMean')
                    set_value = true;
                    try
                        validateattributes(varargin{counter + 1}, {'numeric'}, {'scalar', 'positive', 'nonempty', 'finite'});
                    catch
                        WRN_MSG = 'User Specified signalMean is not valid, using default value';
                        warning(WRN_MSG);
                        set_value = false;
                    end
                    
                    if set_value
                        obj.signalMean = varargin{counter+ 1};
                    end
                end
            end
            
            %If image space signal mean was not set, determine it:
            if isempty(obj.signalMean)
                %determine signalMean from a subset of points in the center
                %of the phantom:
                numslices = round(obj.thePhantom.xyzDims(3)*1/10);
                if numslices < 1
                    numslices = 1;
                end
                
                theslices = [ceil(obj.thePhantom.xyzDims(3)/2)- floor(numslices/2):ceil(obj.thePhantom.xyzDims(3)/2)+ceil(numslices/2)-1];
                
                thes0 = obj.imageData(:,:,theslices);
                thes0 = thes0(:);
                thes0(abs(thes0) == 0) = [];
                obj.signalMean = abs(mean(thes0));
            end
        
            
%             %This is really a global function, not a tissue or layer
%             %specific function so only apply it once:
%             if ~isequal(nTissue, length(obj.thePhantom.theTissues)) || ~isequal(nLayer, 1)
%                 return;
%             end
%             
%             disp('......Adding White Noise');
%             
%             %Initialize values:
%             SNR = 20;
%             signalMean = [];
%             
%             %Sort out any optional input data:           
%             for counter = 1:length(varargin)
%                 if ~ischar(varargin{counter})
%                     continue;
%                 end
%                 
%                 if strcmp(varargin{counter}, 'SNR')
%                     set_value = true;
%                     try
%                         validateattributes(varargin{counter + 1}, {'numeric'}, {'scalar', 'positive', 'nonempty', 'finite'});
%                     catch
%                         WRN_MSG = 'User Specified SNR is not valid, using default value';
%                         warning(WRN_MSG);
%                         set_value = false;
%                     end
%                     
%                     if set_value
%                         SNR = varargin{counter+ 1};
%                     end
%                 end
%                 
%                 if strcmp(varargin{counter}, 'signalMean')
%                     set_value = true;
%                     try
%                         validateattributes(varargin{counter + 1}, {'numeric'}, {'scalar', 'positive', 'nonempty', 'finite'});
%                     catch
%                         WRN_MSG = 'User Specified signalMean is not valid, using default value';
%                         warning(WRN_MSG);
%                         set_value = false;
%                     end
%                     
%                     if set_value
%                         signalMean = varargin{counter+ 1};
%                     end
%                 end
%             end
%                         
%             %If image space signal mean was not set, determine it:
%             if isempty(signalMean)
%                 %determine signalMean from a subset of points in the center
%                 %of the phantom:
%                 numslices = round(obj.thePhantom.xyzDims(3)*1/10);
%                 if numslices < 1
%                     numslices = 1;
%                 end
%                 
%                 theslices = [ceil(obj.thePhantom.xyzDims(3)/2)- floor(numslices/2):ceil(obj.thePhantom.xyzDims(3)/2)+ceil(numslices/2)-1];
%                 
%                 thes0 = obj.imageData(:,:,theslices);
%                 thes0 = thes0(:);
%                 thes0(abs(thes0) == 0) = [];
%                 signalMean = mean(thes0);
%             end
%           
%             %Add the noise to the kdata
%             variance = abs(signalMean)/SNR;
%             
%             %Scale the variance by the matrix size:
%             variance = variance * obj.theMRSystem.finalMatrixSize(1)*obj.theMRSystem.finalMatrixSize(2)*obj.theMRSystem.finalMatrixSize(3);
%             
%             %Generate normally distributed white gaussian noise:
%             noise_matrix = complex(randn(size(thekdata)), randn(size(thekdata)))*1/sqrt(2)*sqrt(variance);
%             thekdata = thekdata + noise_matrix;
            
        end
       
        
    end
    
    methods (Access = 'protected')
       
        function [thekdata] = obtainKspace(obj, imageData)
            
            %Determine if standard Matlab fftn is being used or
            %if the user has called a custom function:
            if strcmp(func2str(obj.Simfft), 'fftn') || strcmp(func2str(obj.Simfft), 'fft2')
                %Standard fftn function:
                if obj.theMRSystem.is2D
                    imageData = fftshift(fftshift(imageData, 1), 2);
                else
                    imageData = fftshift(fftshift(fftshift(imageData, 1), 2), 3);
                end
                
                data = obj.Simfft(imageData);
                
                %Obtain the desired output points in order:
                thekdata = data(sub2ind(size(imageData), squeeze(obj.simKspacePoints(:,2,:)), squeeze(obj.simKspacePoints(:,1,:)), squeeze(obj.simKspacePoints(:,3,:)))); 
            
            elseif strcmp(func2str(obj.Simfft), 'iFGG_2d_type2')
                %Built in 2D non cartesian function:
                
                %Ensure the mex file has been compiled first:
                if ~exist('FGG_Convolution2D_type2.mexa64', 'file')
                    ERR_MSG = ['c  file for built in 3D noncartesian fourier transform has not been compiled. \n'...
                                'Either set your own function for the noncartesian fourier transform or \n'...
                                'ensure the digital breast phantom files are on your file path and type: \n \n '...
                                '    mex FGG_Convolution2D_type2.c %s  \n'];
                    error(ERR_MSG, '');
                end
                
                display(['......Performing FT of image data using function ' func2str(obj.Simfft)]);
                
                %Locate all the projections from a given slice:
                % ** Assume all slices are in the z direction ** %
                kz_vals = squeeze(obj.theMRSystem.kspacePoints(1,3,:));
                [peaks, slice_locs] = findpeaks(diff(kz_vals));
                slice_locs = [1; slice_locs; size(obj.theMRSystem.kspacePoints, 3)];
                
                thekdata = zeros([size(obj.theMRSystem.kspacePoints, 1), size(obj.theMRSystem.kspacePoints, 3)]);
                
                for counter = 1:length(slice_locs)-1
                    
                    %calculate the desired slice number
                    current_slice = floor((obj.theMRSystem.kspacePoints(1,3, slice_locs(counter)) - obj.theMRSystem.kmin(3))/obj.theMRSystem.dk(3) + 1);
                    
                    % Value for setting up matrix:
                    desired_size = [size(obj.theMRSystem.kspacePoints, 1)*(slice_locs(counter+1) - slice_locs(counter)), 1];
                    
                    %Reshape the kspace points:
                    theimageData = imageData(:,:,current_slice);
                    
                    theSlicekdata = obj.Simfft(permute(theimageData, [2, 1]),[obj.theMRSystem.kmax(1:2); obj.theMRSystem.kmin(1:2); ...
                        [reshape(obj.theMRSystem.kspacePoints(:,1,slice_locs(counter):(slice_locs(counter+1) - 1)), desired_size),...
                         reshape(obj.theMRSystem.kspacePoints(:,2,slice_locs(counter):(slice_locs(counter+1) - 1)), desired_size)]],...
                         6);
                    
                    %remove the first 2 points (kmax and kmin) which were added:
                    theSlicekdata(1:2) = []; 
                    
                    %reshape the data back to the original shape:
                    theSlicekdata = reshape(theSlicekdata, [size(obj.theMRSystem.kspacePoints, 1), (slice_locs(counter+1) - slice_locs(counter))]);
                    
                    thekdata(:,slice_locs(counter):slice_locs(counter+1)-1) = theSlicekdata;
                end
                    
            elseif strcmp(func2str(obj.Simfft), 'iFGG_3d_type2')
                %Built in 3D non cartesian function:
                
                %Ensure the mex file has been compiled first:
                if ~exist('FGG_Convolution3D_type2.mexa64', 'file')
                    ERR_MSG = ['c  file for built in 3D noncartesian fourier transform has not been compiled. \n'...
                                'Either set your own function for the noncartesian fourier transform or \n'...
                                'ensure the digital breast phantom files are on your file path and type: \n \n '...
                                '    mex FGG_Convolution3D_type2.c %s  \n'];
                    error(ERR_MSG, '');
                end
                
                display(['......Performing FT of image data using function ' func2str(obj.Simfft)]);

                desired_size = [size(obj.theMRSystem.kspacePoints, 1)*size(obj.theMRSystem.kspacePoints,3), 1];
                tic
                thekdata = obj.Simfft(permute(imageData, [2, 1, 3]),[obj.theMRSystem.kmax; obj.theMRSystem.kmin; ...
                    [reshape(obj.theMRSystem.kspacePoints(:,1,:), desired_size),...
                     reshape(obj.theMRSystem.kspacePoints(:,2,:), desired_size),...
                     reshape(obj.theMRSystem.kspacePoints(:,3,:), desired_size)]],...
                     6);
                timeiFGG_3d_type2 = toc;
                disp(['Time Forward NUFFT: [' 8 '', num2str(timeiFGG_3d_type2),']' 8 ' sec for a [',num2str(size(imageData)),'] Volume']);
                %remove the first 2 points (kmax and kmin) which were added:
                thekdata(1:2) = [];
                %reshape the data back to the original shape:
                thekdata = reshape(thekdata, [size(obj.theMRSystem.kspacePoints, 1), size(obj.theMRSystem.kspacePoints, 3)]);
                
            elseif strcmp(func2str(obj.Simfft), 'gpuNUFFT')
                %Use GPU NUFFT function:
                
                %Ensure the mex file has been compiled first:
                if ~exist('mex_gpuNUFFT_forw_f.mexw64', 'file')
                    ERR_MSG = ['The executable for the GPU based 3D noncartesian fourier transform has not been compiled. \n'...
                    'Disable the GPU fuctionality (set obj.useGPU = false) or compiled the function for the GPU NUFFT \n'];
                    error(ERR_MSG, '');
                end
                
                
                display(['......Performing FT of image data using function ' func2str(obj.Simfft)]);
                [Spts,nDim,Proj]=size(obj.theMRSystem.kspacePoints);
                k = permute(obj.theMRSystem.kspacePoints, [1, 3, 2]);
                k = reshape(k,[Spts*Proj,nDim]);
                k = ((k - min(k(:)))./ (max(k(:))- min(k(:))))-.5; % range [-.5,.5]
                W = ones(Spts,Proj); % no weighting 
                if(size(imageData,1)>200), osf = 1.25; else, osf = 2; end
                wg = 3; sw = 12;
                FT = obj.Simfft(k',W,osf,wg,sw,size(imageData),[],true);
                FTtime = tic;
                thekdata = FT*imageData;
                timegpuFFT = toc(FTtime);
                disp(['Time Forward NUFFT: [' 8 '', num2str(timegpuFFT),']' 8 ' sec for a [',num2str(size(imageData)),'] Volume']);
                thekdata = reshape(thekdata,[Spts,Proj]);
            else
                %User defined FT function:
                display(['......Performing FT of image data using function ' func2str(obj.Simfft)]);
                
                LoopCounters.TissueCounter = obj.TissueCounter;
                LoopCounters.LayerCounter = obj.LayerCounter;
                LoopCounters.CoilCounter = obj.CoilCounter;
                
                thekdata = obj.Simfft(imageData, obj.theMRSystem, LoopCounters);
                
                %Verify the function has returned a matrix of the correct
                %size:
                if ~isequal([size(obj.theMRSystem.kspacePoints, 1), size(obj.theMRSystem.kspacePoints, 3)], size(thekdata))
                    ERRMSG = ['ERROR: the kspace data obtained via the user defined fft function ' func2str(obj.Simfft) ' must be arranged points per TR by numTRs.'];
                    error(ERRMSG);
                end

            end

        end
        
        function the_points = fftshift_points(obj, the_points, kmax)
        %Adjust coordinates of kspace points to account for fftshift
        %operation:
        
            the_points = the_points + repmat(floor(kmax*0.5), [size(the_points, 1), 1, size(the_points, 3)]);
           
            if ~obj.theMRSystem.is2D
                maxdim = 3;
            else
                maxdim = 2;
            end
            
            for counter = 1:maxdim
                temp_points = the_points(:,counter, :);
                temp_points(temp_points>kmax(counter))=temp_points(temp_points>kmax(counter))-kmax(counter);
                the_points(:,counter, :) = temp_points;
            end

        end
        
        function tissueProps = grabTissueParams(obj, tissueNum, layerNum)
            %This function collects all the parameters of the particular
            %tissue/tissue layer except the images and puts them in a
            %structure that can be passed to other functions:
            
            tissueProps = struct;
            
            %Add relevent properties from the overall Phantom:
            theprops = properties(Phantom);
            numTissues = length(obj.thePhantom.theTissues);
            
            for counter = 1:length(theprops)
                
%                 if strcmp(theprops{counter}, 'patientPosition')
%                     continue;
                if strcmp(theprops{counter}, 'partialVolume_Layer')
                    continue;
                elseif strcmp(theprops{counter}, 'partialVolume_Tissue')
                    continue;
%                 elseif strcmp(theprops{counter}, 'xyzDims')
%                     continue;
                elseif strcmp(theprops{counter}, 'theTissues')
                    continue;
                elseif strcmp(theprops{counter}, 'tissueNames')
                    continue;    
                end
                
                tempprop = obj.thePhantom.(theprops{counter});
                
                if length(tempprop) == numTissues
                    tempprop = tempprop(tissueNum);
                end
                
                tissueProps.(theprops{counter}) = tempprop;
                
            end
            
            %Add properties from specific tissue layer:
            theprops = properties(Tissue);
            numLayers = length(obj.thePhantom.theTissues{tissueNum}.tissueLayers);
            
            for counter = 1:length(theprops)
                
                if strcmp(theprops{counter}, 'tissueLayers')
                    continue;
                elseif strcmp(theprops{counter}, 's0')
                    continue;
                elseif strcmp(theprops{counter}, 'noiseMask')
                    continue;
                elseif strcmp(theprops{counter}, 'modelParams')
                    
                    %Grab the model parameters corresponding to the desired
                    %layer:
                    modelParams = struct;
                    thefields = fieldnames(obj.thePhantom.theTissues{tissueNum}.modelParams);
                    
                    for counter2 = 1:length(thefields)
                        modelParams.(thefields{counter2}) = obj.thePhantom.theTissues{tissueNum}.modelParams.(thefields{counter2}){layerNum};
                    end
                    
                    tissueProps.modelParams = modelParams;
                    
                else

                    tempprop = obj.thePhantom.theTissues{tissueNum}.(theprops{counter});
                
                    if length(tempprop) == numLayers
                        tempprop = tempprop(layerNum);
                    end
                
                    tissueProps.(theprops{counter}) = tempprop;
                end
            end
        end

        function [thekdata]=applyUserDCEModel(obj, thekdata, varargin)
            %This function manages the calling of the user specified DCE
            %function and applies the results to the MR kspace data
            
            %Get the function handle and variable call list out of the
            %userDCEFuncs structure:
            fhandle = obj.userDCEFuncs{obj.pipelineCounter}(1,1);
            var_call = obj.userDCEFuncs{obj.pipelineCounter}(1,2);
            
            %Collect the variables necessary to call the users DCE function:
            input_variables = cell([1, length(var_cal)]);
            
            for counter = 1:length(input_variables)
                input_variables(counter) = {eval(var_cal{counter})};
            end
            
            %Call the users function:
            [time, PercentEnhancement] = fhandle(input_variables{:});
            
            %Verify the returned data:
            validateattributes(time, {'numeric'}, {'vector', 'increasing', 'numel', numel(PercentEnhancement), 'nonnegative'});
            validateattributes(PercentEnhancement, {'numeric'}, {'vector'});
            
            %Apply the time curve to the MRdata:
            
            %interpolate the curve to match the simulation time curve:
            PercentEnhancement = interp1(time, PercentEnhancement, obj.theMRSystem.time);
            
            %Apply the curve to the kdata:
            thekdata = thekdata + thekdata .* repmat(PercentEnhancement, [size(thekdata, 1), 1]);

        end
        
        function reorderUserData(obj, oldLoc, newLoc)
            
            %Remove the data from its original location:
            tempdata = obj.userDCEFuncs(oldLoc);
            obj.userDCEFuncs(oldLoc) = [];

            if newLoc > oldLoc
                newLoc = newLoc - 1;
            end
            
            %put it in its new location:
            if newLoc == 1
                obj.userDCEFuncs = [tempdata, obj.userDCEFuncs];
            elseif newLoc > length(obj.(pipename))
                obj.userDCEFuncs(newLoc) = tempdata;
            else
                obj.userDCEFuncs = [obj.userDCEFuncs(1:newLoc - 1), tempdata, obj.userDCEFuncs(newLoc:length(obj.userDCEFuncs))];
            end
        end
        
        function paramCall = searchForInputs(obj, possibleInputs)
            
            if ~isempty(possibleInputs)

                %verify that all inputs are strings:
                for counter = 1:length(possibleInputs)
                    if ~ischar(possibleInputs{counter})
                        ERR_MSG = ['ERROR: All of the optional, user specified parameters must be strings \n'...
                                    'corresponding to properties the Phantom or MRSystem objects contained in MRSim. %s'];
                        error(ERR_MSG, '');
                    end
                end

                %set up arrays to organize parameters:
                paramCall = {[length(possibleInputs), 1]};

                %Obtain all possible parameters from the MRSim object
                loop_counters = {'nCoil', 'nTissue', 'nLayer'};
                phantom_parameters = properties(obj.thePhantom);
                phantom_user_params = fieldnames(obj.thePhantom.userParams);
                MRsystem_parameters = properties(obj.theMRSystem);
                MRSystem_user_params = fieldnames(obj.theMRSystem.userParams);
                tissue_parameters = properties(Tissue);
                tissue_user_params = {};
                tissue_model_params = {};

                %Obtain a list of tissue parameters 
                for counter = 1:length(obj.thePhantom.theTissues)

                    %If this is the first tissue, add all model parameters
                    %to the list:
                    if counter == 1
                        tissue_model_params = fieldnames(obj.thePhantom.theTissues{counter}.modelParams);
                        tissue_user_params = fieldnames(obj.thePhantom.theTissues{counter}.userParams);
                    else
                        %If this is not the first tissue, only add names of
                        %model parameters not already in the list
                        temp_model_params = fieldnames(obj.thePhantom.theTissues{counter}.modelParams);
                        for counter2 = 1:length(temp_model_params)
                            %check to see if the current parameter is in
                            %the list already:
                            try validatestring(temp_model_params{counter2}, tissue_model_params);
                            catch
                                %The current parameter was not in the list
                                %already so add it:
                                tissue_model_params(length(tissue_model_params) + 1) = temp_model_params(counter2);
                            end
                        end

                        %If this is not the first tissue, only add names of
                        %user parameters not already in the list
                        temp_user_params = fieldnames(obj.thePhantom.theTissues{counter}.userParams);
                        for counter2 = 1:length(temp_user_params)
                            %check to see if the current parameter is in
                            %the list already:
                            try validatestring(temp_user_params{counter2}, tissue_user_params);
                            catch
                                %The current parameter was not in the list
                                %already so add it:
                                tissue_user_params(length(tissue_user_params) + 1) = temp_user_params(counter2);
                            end
                        end
                    end
                end
                
                %Return a list of all possible parameters if the input is
                %'-h':
                if strcmp(possibleInputs{1}, '-h')
                    
                    SectionNames = {'loop counters',...
                            'Phantom Parameters',...
                            'Phantom User Defined Parameters',...
                            'MR System Parameters',...
                            'MR System User Defined Parameters', ...
                            'Tissue Parameters',...
                            'User Defined Tissue Parameters', ...
                            'Tissue Modeling Parameters'};
                    
                    paramCall = {SectionNames;...
                                 loop_counters;...
                                 phantom_parameters; ... 
                                 phantom_user_params; ...
                                 MRsystem_parameters;...
                                 MRSystem_user_params;...
                                 tissue_parameters;...
                                 tissue_user_params;...
                                 tissue_model_params};  
                     return;
                     
                     
                end
                

                %Now search through all the parameters to identify the
                %parameters requested by the users function:
                for counter = 1:length(possibleInputs)

                    %Check for loop counters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, loop_counters);
                    catch
                        valid_string = '';
                    end
                    
                    if ~isempty(valid_string)
                        paramCall(counter) = {[possibleInputs{counter}]};
                        continue;
                    end
                    
                    %Check for phantom parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, phantom_parameters);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.thePhantom.' possibleInputs{counter}]};
                        continue;
                    end

                    %Check for phantom user parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, phantom_user_params);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.thePhantom.userParams.' possibleInputs{counter}]};
                        continue;
                    end 

                    %Check for MR System parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, MRsystem_parameters);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.theMRSystem.' possibleInputs{counter}]};
                        continue;
                    end  

                    %Check for MR System user parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, MRSystem_user_params);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.theMRSystem.userParams' possibleInputs{counter}]};
                        continue;
                    end  

                    %Check for the Tissue parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, tissue_parameters);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        if isequal(length(obj.thePhantom.theTissues{1}.(possibleInputs{counter})), length(obj.thePhantom.theTissues{1}.tissueLayers))
                            if iscell(obj.thePhantom.theTissues{1}.(possibleInputs{counter}))
                                paramCall(counter) = {['obj.thePhantom.theTissues{nTissue}.' possibleInputs{counter} '{nLayer}']};
                            else
                                paramCall(counter) = {['obj.thePhantom.theTissues{nTissue}.' possibleInputs{counter} '(nLayer)']};
                            end
                        else
                            paramCall(counter) = {['obj.thePhantom.theTissues{nTissue}.' possibleInputs{counter}]};
                        end
                        continue;
                    end 

                    %Check for tissue user parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, tissue_user_params);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.thePhantom.theTissue{nTissue}.userParams.' possibleInputs{counter}]};
                        continue;
                    end 


                    %Check for tissue model parameters:
                    try 
                        valid_string = validatestring(possibleInputs{counter}, tissue_model_params);
                    catch
                        valid_string = '';
                    end

                    if ~isempty(valid_string)
                        paramCall(counter) = {['obj.thePhantom.theTissues{nTissue}.modelParams.' possibleInputs{counter} '.{nLayer}']};
                        continue;
                    end 

                    %If you make it to the end of the loop, the desired
                    %parameter was not found, throw error:
                    ERR_MSG = ['ERROR in applyUserDCEModel: specified user parameter ' possibleInputs{counter} ' not found. \n'...
                                'Use the displayPossibleInputs function to view a list of all possible inputs. %s'];
                    error(ERR_MSG, '');

                end
            else
                paramCall = {};
            end
            
        end
                        
    end
    
    
end