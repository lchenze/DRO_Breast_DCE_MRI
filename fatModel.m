

classdef fatModel
    
    %Properties
    properties(Dependent)
        
       %Selected Model Type
        ModelType; 
    end

    properties
        
        %Temperature in degrees C
        Temp = 32.073; %Temperature that makes no adjustment to the model
    end
    
    properties (SetAccess = 'protected')
        %Referenced fat models
        ReferenceModels = struct;
        
        %Absolute ppm of models
        %Model_ppm = [5.29, 5.19, 4.2, 2.75, 2.2, 2.02, 1.6, 1.3, 0.9];  %Absolute ppm of models
    end
    
    properties (SetAccess = 'protected', Dependent)
        
        Model_Relative_ppm = [];
        Model_Relative_Amplitude = [];
        Model_Absolute_ppm = [];
        
    end
    
    properties (SetAccess = 'protected', Hidden)
       
        model_number = 2;
        
    end
        
        % ----------------- Class Constructor ---------------------------%
    methods
        
        function obj = fatModel(varargin)
            
            %Set up the reference fat models:
            obj = obj.setUpReferenceModels;
            
            %Sort through the input as appropriate:
            if ~nargin
                return;
            end
            
            counter = 1;
            while counter <= length(varargin)
                if ~ischar(varargin{counter})
                    counter = counter + 1;
                    continue;
                end
                
                validatestring(varargin{counter}, properties(fatModel));
                
                obj.(varargin{counter}) = varargin{counter+1};
                counter = counter + 2;
                
            end
            
                                    
        end
    end
    
     % ----------------- Set Methods ---------------------------%
    methods
        
        function obj=set.ModelType(obj, str)
            
            %Verify a string has been passed:
            validateattributes(str, {'char'}, {'2d'});
                
            %Get the names of the reference models:
            [model_names, model_names_full] = getReferenceModelNames(obj);
            
            %Match string with Reference Model:
            model_num = 0;
            for counter = 1:length(model_names)
                if strcmp(str, model_names{counter}) || strcmp(str, model_names_full{counter})
                    model_num = counter;
                    break;
                end
            end
                        
            if model_num == 0
                error_msg = ['Invalid model name selected.  Selected model name is: ' str newline ...
                             newline ...
                             'Valid names are: ' newline];
                for counter = 1:length(model_names)
                    error_msg = [error_msg '      ' model_names{counter} newline];
                end
                               
                error(error_msg);
            else
                obj.model_number = model_num;
            end

        end %End set ModelType
        
        
        function obj=set.Temp(obj, val)
            
            %if an empty value is passed, do not reset the temp
            if ~isempty(val)

                %Verify a string has been passed:
                validateattributes(val, {'numeric'}, {'2d', 'nonnegative', 'nonnan', 'numel', 1, '<' 100});

                obj.Temp = val;
            end
            
        end %End set Temp  
        
        
    end
    
    % ----------------- Get Methods ---------------------------%
    methods
        
        function modelname = get.ModelType(obj)
            
            %Get the names of the reference models:
            [~, model_names_full] = getReferenceModelNames(obj);
            
            modelname = model_names_full{obj.model_number};
                        
        end
        
        function absolute_ppm = get.Model_Absolute_ppm(obj)
        
            absolute_ppm = obj.ReferenceModels(obj.model_number).Model_ppm;
            
        end  %End get Model_Relative_ppm
            
        function relative_ppm = get.Model_Relative_ppm(obj)
        
            [deltaPpm] = obj.TemperaturePpmShift(obj.Temp);
            relative_ppm = obj.ReferenceModels(obj.model_number).Model_ppm + deltaPpm-4.7;
            
        end  %End get Model_Relative_ppm
        
        function relative_amp = get.Model_Relative_Amplitude(obj)

            if obj.model_number <= 19
                amplitudes = [2*obj.ReferenceModels(obj.model_number).ndb, ...
                              1, ...
                              4, ...
                              2*obj.ReferenceModels(obj.model_number).nmidb, ...
                              6, ...
                              4*(obj.ReferenceModels(obj.model_number).ndb - obj.ReferenceModels(obj.model_number).nmidb), ...
                              6, ...
                              6*(obj.ReferenceModels(obj.model_number).CL - 4) - 8*obj.ReferenceModels(obj.model_number).ndb + 2*obj.ReferenceModels(obj.model_number).nmidb, ...
                              9];

                relative_amp = amplitudes * 1/sum(amplitudes(:));
            else
                relative_amp = obj.ReferenceModels(obj.model_number).Model_Relative_Amplitude;
            end
            
        end %End get Model_Relative_Amplitude

       % Model_Relative_Amplitue = [];
        
        
        
    end
    
     % ----------------- General Methods ---------------------------%
    methods
        
        function [deltaPpm] = TemperaturePpmShift(obj, Temp)
            
            %Calculated based on the equation given in:
            %Hernando D, Sharma SD, Kramer H, Reeder SB.  On the
            %Confounding Effect of Temperature on Chemical Shift-Encoded
            %Fat Quantification.  MRM 72:464-470(2014)
            
            ppm = -1.085*10^-2 * Temp + 3.748;
            deltaPpm = ppm - 3.4;
            
        end %end TemperaturePpmShift
        
    end
    
    %----------------- Protected Methods ---------------------------%
    methods (Access = 'protected', Hidden = true) 
        
        function obj=setUpReferenceModels(obj)
            
            %Reference fat models:
            %% Models are obtained from: Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046.
                        
            %Peanut Oil
            obj.ReferenceModels(1).name = 'PeanutOil';
            obj.ReferenceModels(1).fullname = 'Peanut Oil';
            obj.ReferenceModels(1).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(1).CL = 17.97;
            obj.ReferenceModels(1).ndb = 3.48;
            obj.ReferenceModels(1).nmidb = 1.01;
                        
            %Human Subcutaneous Fat
            obj.ReferenceModels(2).name = 'SubFat';
            obj.ReferenceModels(2).fullname = 'Human Subcutaneous Fat';
            obj.ReferenceModels(2).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(2).CL = 17.48;
            obj.ReferenceModels(2).ndb = 2.88;
            obj.ReferenceModels(2).nmidb = 0.70;
            
            %Human liver Fat
            obj.ReferenceModels(3).name = 'LiverFat';
            obj.ReferenceModels(3).fullname = 'Human Liver Fat';
            obj.ReferenceModels(3).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(3).CL = 17.45;
            obj.ReferenceModels(3).ndb = 1.92;
            obj.ReferenceModels(3).nmidb = 0.32;
            
            %Human Marrow
            obj.ReferenceModels(4).name = 'Marrow';
            obj.ReferenceModels(4).fullname = 'Human Marrow';
            obj.ReferenceModels(4).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(4).CL = 17.33;
            obj.ReferenceModels(4).ndb = 2.86;
            obj.ReferenceModels(4).nmidb = 0.74;
            
            %Palm Oil
            obj.ReferenceModels(5).name = 'PalmOil';
            obj.ReferenceModels(5).fullname = 'Palm Oil';
            obj.ReferenceModels(5).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(5).CL = 16.98;
            obj.ReferenceModels(5).ndb = 1.68;
            obj.ReferenceModels(5).nmidb = 0.29;

            %Beef Tallow
            obj.ReferenceModels(6).name = 'TallowB';
            obj.ReferenceModels(6).fullname = 'Tallow (beef)';
            obj.ReferenceModels(6).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(6).CL = 17.07;
            obj.ReferenceModels(6).ndb = 1.52;
            obj.ReferenceModels(6).nmidb = 0.13;
            
            %Mutton Tallow
            obj.ReferenceModels(7).name = 'TallowM';
            obj.ReferenceModels(7).fullname = 'Tallow (mutton)';
            obj.ReferenceModels(7).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(7).CL = 17.26;
            obj.ReferenceModels(7).ndb = 1.82;
            obj.ReferenceModels(7).nmidb = 0.32;
            
            %Pork Lard
            obj.ReferenceModels(8).name = 'Lard';
            obj.ReferenceModels(8).fullname = 'Lard (pork)';
            obj.ReferenceModels(8).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(8).CL = 17.33;
            obj.ReferenceModels(8).ndb = 2.10;
            obj.ReferenceModels(8).nmidb = 0.38;
            
            %Chicken Fat
            obj.ReferenceModels(9).name = 'Chicken';
            obj.ReferenceModels(9).fullname = 'Chicken Fat';
            obj.ReferenceModels(9).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(9).CL = 17.34;
            obj.ReferenceModels(9).ndb = 2.72;
            obj.ReferenceModels(9).nmidb = 0.68;
            
            %Cocoa butter
            obj.ReferenceModels(10).name = 'Cocoa';
            obj.ReferenceModels(10).fullname = 'Cocoa Butter';
            obj.ReferenceModels(10).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(10).CL = 17.41;
            obj.ReferenceModels(10).ndb = 1.20;
            obj.ReferenceModels(10).nmidb = 0.09;
            
            %Cottonseed Oil
            obj.ReferenceModels(11).name = 'Cotton';
            obj.ReferenceModels(11).fullname = 'Cottonseed Oil';
            obj.ReferenceModels(11).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(11).CL = 17.43;
            obj.ReferenceModels(11).ndb = 3.74;
            obj.ReferenceModels(11).nmidb = 1.60;
            
            %Olive Oil
            obj.ReferenceModels(12).name = 'Olive';
            obj.ReferenceModels(12).fullname = 'Olive Oil';
            obj.ReferenceModels(12).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(12).CL = 17.74;
            obj.ReferenceModels(12).ndb = 2.89;
            obj.ReferenceModels(12).nmidb = 0.35;
            
            %Corn Oil
            obj.ReferenceModels(13).name = 'Corn';
            obj.ReferenceModels(13).fullname = 'Corn Oil';
            obj.ReferenceModels(13).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(13).CL = 17.76;
            obj.ReferenceModels(13).ndb = 4.31;
            obj.ReferenceModels(13).nmidb = 1.75;
            
            %Soybean Oil
            obj.ReferenceModels(14).name = 'Soy';
            obj.ReferenceModels(14).fullname = 'Soybean Oil';
            obj.ReferenceModels(14).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(14).CL = 17.79;
            obj.ReferenceModels(14).ndb = 4.49;
            obj.ReferenceModels(14).nmidb = 2.00;
            
            %Sesame Oil
            obj.ReferenceModels(15).name = 'Sesame';
            obj.ReferenceModels(15).fullname = 'Sesame Oil';
            obj.ReferenceModels(15).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(15).CL = 17.80;
            obj.ReferenceModels(15).ndb = 3.86;
            obj.ReferenceModels(15).nmidb = 1.32;
            
            %Walnut Oil
            obj.ReferenceModels(16).name = 'Walnut';
            obj.ReferenceModels(16).fullname = 'Walnut Oil';
            obj.ReferenceModels(16).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(16).CL = 17.84;
            obj.ReferenceModels(16).ndb = 5.02;
            obj.ReferenceModels(16).nmidb = 2.32;
            
            %Safflower Oil
            obj.ReferenceModels(17).name = 'Safflower';
            obj.ReferenceModels(17).fullname = 'Safflower Oil';
            obj.ReferenceModels(17).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(17).CL = 17.9;
            obj.ReferenceModels(17).ndb = 5.14;
            obj.ReferenceModels(17).nmidb = 2.34;
            
            %Sunflower Oil
            obj.ReferenceModels(18).name = 'Sunflower';
            obj.ReferenceModels(18).fullname = 'Sunflower Oil';
            obj.ReferenceModels(18).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(18).CL = 17.94;
            obj.ReferenceModels(18).ndb = 3.62;
            obj.ReferenceModels(18).nmidb = 0.92;
            
            %Canola Oil
            obj.ReferenceModels(19).name = 'Canola';
            obj.ReferenceModels(19).fullname = 'Canola Oil';
            obj.ReferenceModels(19).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(19).CL = 17.95;
            obj.ReferenceModels(19).ndb = 3.91;
            obj.ReferenceModels(19).nmidb = 1.14;
            
            %Cod Liver Oil
            obj.ReferenceModels(19).name = 'Cod';
            obj.ReferenceModels(19).fullname = 'Cod Liver Oil';
            obj.ReferenceModels(19).source = 'Bydder M, Girard O, Hamilton G. Mapping the double bonds in triglycerides. Magn Reson Imaging 2011;29:10411046';
            obj.ReferenceModels(19).CL = 18.62;
            obj.ReferenceModels(19).ndb = 5.09;
            obj.ReferenceModels(19).nmidb = 2.86;
            
            %Set absolute model ppm for all models:
            for counter = 1:19
                obj.ReferenceModels(counter).Model_ppm = [5.29, 5.19, 4.2, 2.75, 2.2, 2.02, 1.6, 1.3, 0.9];
            end
            
            %% Simple single peak model:
            obj.ReferenceModels(20).name = 'SinglePeak';
            obj.ReferenceModels(20).fullname = 'Single Peak Fat Model';
            obj.ReferenceModels(20).source = 'PPM difference between fat and water';
            obj.ReferenceModels(20).Model_Relative_Amplitude = 1;
            obj.ReferenceModels(20).Model_ppm = 1.2;
            
            %% Six Peak peanut Oil model:
            obj.ReferenceModels(21).name = 'PeanutOil-6p';
            obj.ReferenceModels(21).fullname = 'Peanut Oil 6 peak';
            obj.ReferenceModels(21).source = 'Yu H, Shimakawa A, McKenzi C, et al. Multiecho Water-Fat Separation and Simultaneous R2* Estimation with Multifrequency Fat Spectrum Modeling Magn Res Med 60:1122-1134 (2008)';
            obj.ReferenceModels(21).Model_Relative_Amplitude = [0.62, 0.15, 0.10, 0.06, 0.03, 0.04];
            obj.ReferenceModels(21).Model_ppm = [-420, -318, 94, -472, -234, -46]*1/(10^-6*42.577*10^6 * 3)+4.7;
            
            %% Six Peak Hamilton Liver Fat Model
            obj.ReferenceModels(22).name = 'LiverFat-6p';
            obj.ReferenceModels(22).fullname = 'Liver Fat 6 peak';
            obj.ReferenceModels(22).source = 'Hamilton G, Yokoo T, Bydder M, et al. In vivo characterization of the liver fat 1H MR Spectrum"';
            obj.ReferenceModels(22).Model_Relative_Amplitude = [0.047, 0.039, 0.006, 0.12, 0.7, 0.088];
            obj.ReferenceModels(22).Model_ppm = [5.3, 4.2, 2.75, 2.1, 1.3, 0.9];
                       
            
            
        end
        
        function [model_names, model_names_full] = getReferenceModelNames(obj)
           
            model_names = cell([length(obj.ReferenceModels), 1]);
            model_names_full = cell([length(obj.ReferenceModels), 1]);
            
            for counter = 1:length(obj.ReferenceModels)
                model_names(counter) = {obj.ReferenceModels(counter).name};
                model_names_full(counter) = {obj.ReferenceModels(counter).fullname};
            end
                
        end
        
        
    end
    
end





