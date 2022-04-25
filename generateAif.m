function [bloodConcentration] = generateAif(time, varargin)
%generateAif generates a simulated arterial input function (AIF) for the 
%digital reference object (DRO)
%
%   The generateAif function produces a simulated AIF based on the 
%   publication by Parker et al Magn Reson Med Nov 2006 [1].In Parker 
%   et al, the population-averaged blood concentration (Cb(t)) is defined as:
%
%       Cb(t) = Sum[A_n/(sigma_n*sqrt(2 *pi))*exp(-(t-T_n)^2/(2*sigm_n^2)]
%                   + alpha * exp(-beta * t) / (1 + exp(-s(t-tau)))
%     
%   Where:
%       n       - Is the number of the Gaussian used to represent the AIF.  
%                   Two are used by default
%       A_n     - Is a scaling constant of the nth Gaussian. 
%                   units = mmol.min
%       T_n     - Is the center of the nth Gaussian.  units = min
%       sigma_n - Is the width of the nth Gaussian.  units = min
%       alpha   - Is the amplitude constant of the exponential.  
%                   units = mmol.
%       beta    - Is the decay constant of the exponential. units = min^-1
%       s       - Is the width of the sigmoid   units = min^-1
%       tau     - Is the center of the sigmoid. units = min
%
%   
%   bloodConcentration = generateAif(time) generates an AIF function 
%   based on the default parameter values. Input and output variables are 
%   defined as follows:
%       bloodConcentration - Is the calculated concentration of the 
%                            contrast agent in the blood in units of mol
%       time - Is the time in seconds after the injection of the contrast 
%              agent
%
%   bloodConcentration = generateAif (time, ‘parameterName1’, 
%   parameterValue1, …) generates the an AIF function based on user 
%   defined values for the parameters in the Cb(t) equation.  One or more 
%   of the parameter values may be set using this format.  Default values 
%   for the parameters are taken from Parker et al:  
%       A_n =   [0.809, 0.330]; %units = mmol.min
%       T_n =   [0.17046, 0.365]; %units = min
%       sigma = [0.0563, 0.132]; %units = min
%       alpha = 1.050;  %units = mmol
%       beta =  0.1685; %units = min^-1
%       s =     38.078; %units = min^-1
%       tau =   0.483;  %units = min
%
%   Examples
%   --------
%   Example 1: Obtain the AIF for a time of 0 to 360 seconds (6 minutes):
%       time = 0:0.05:360 
%       [bloodConcentration] = generateAif(time);
%       plot(time, bloodConcentration*1000);
%       xlabel('Time in seconds'), ylabel('C_b (t) (Mol)');
%
%   Example 2: Shift the peak of the AIF back to 75 seconds (1.25 minutes), 
%   more typical of the timing seen in breast imaging when the contrast 
%   agent must travel from the arm, to the heart and back out to the breast 
%   tissue.
%       time = 0:0.05:360;
%       T_n = [1.25, 1.4445]; % units = min
%       tau = [1.5625]; %units = min
%       [bloodConcentration] = generateAif(time, 'T_n', T_n, 'tau', tau);
%       plot(time, bloodConcentration);
%       xlabel('Time in seconds'), ylabel('C_b (t) (Mol)');
%
%   References:
%   [1] Parker GJ, Roberts C, Macdonald A, et al. "Experimentally-derived 
%       functional form for a population-averaged high-temporal-resolution 
%       arterial input function for dynamic contrast-enhanced MRI." Magn 
%       Reson Med 2006;56(5):993-1000.

%   Author: Leah Henze Bancroft
%        	Department of Radiology
%           University of Wisconsin Madison
%           lhenze@wisc.edu
%
%   

    %Default variable assignment:
    A_n =   [0.809, 0.330]; %units = mmol.min
    T_n =   [0.17046, 0.365]; %units = min
    sigma = [0.0563, 0.132]; %units = min
    alpha = 1.050;  %units = mmol
    beta =  0.1685; %units = min^-1
    s =     38.078; %units = min^-1
    tau =   0.483;  %units = min
        
    %Look for user defined variable assignments: 
    if ~isempty(varargin)
        for i_varargin = 1:2:size(varargin, 2) %increment counter by 2 since
                                              %variables added in pairs
            %Check for a character string containing a parameter name for 
            %the Cb(t) function: 
            if ~ischar(varargin{i_varargin})
                err_msg = ['Improper input type.  Enter the desired parameter name as a character string'];
                error(err_msg)
            end
            
            %Match the character string to a parameter name:
            if strcmp(varargin(i_varargin), 'A_n')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to A_n must be numeric';
                    error(err_msg);
                end
                if size(varargin{i_varargin + 1}, 1) > 1 || ndims(varargin{i_varargin + 1})>2
                    err_msg = 'A_n must be a 1xn array';
                    error(err_msg);
                end
                
                %Assign the new value to A_n:
                A_n = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 'T_n')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to T_n must be numeric';
                    error(err_msg);
                end
                if size(varargin{i_varargin + 1}, 1) > 1 || ndims(varargin{i_varargin + 1})>2
                    err_msg = 'T_n must be a 1xn array';
                    error(err_msg);
                end
                
                %Assign the new value to T_n:
                T_n = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 'sigma')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to signa must be numeric';
                    error(err_msg);
                end
                if size(varargin{i_varargin + 1}, 1) > 1 || ndims(varargin{i_varargin + 1})>2
                    err_msg = 'sigma must be a 1xn array';
                    error(err_msg);
                end
                
                %Assign the new value to sigma:
                sigma = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 'alpha')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to alpha must be numeric';
                    error(err_msg);
                end
                if numelements(varargin{i_varargin + 1}) > 1
                    err_msg = 'alpha must be a 1x1 array';
                    error(err_msg);
                end
                
                %Assign the new value to alpha:
                alpha = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 'beta')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to beta must be numeric';
                    error(err_msg);
                end
                if numelements(varargin{i_varargin + 1}) > 1
                    err_msg = 'beta must be a 1x1 array';
                    error(err_msg);
                end
                
                %Assign the new value to beta:
                beta = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 's')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to s must be numeric';
                    error(err_msg);
                end
                if numelements(varargin{i_varargin + 1}) > 1
                    err_msg = 's must be a 1x1 array';
                    error(err_msg);
                end
                
                %Assign the new value to s:
                s = varargin{i_varargin + 1};
            elseif strcmp(varargin(i_varargin), 'tau')
                %Check for inappropriate input values:
                if ~isnumeric(varargin{i_varargin + 1})
                    err_msg = 'Values assigned to tau must be numeric';
                    error(err_msg);
                end
                if numelements(varargin{i_varargin + 1}) > 1
                    err_msg = 'tau must be a 1x1 array';
                    error(err_msg);
                end
                
                %Assign the new value to alpha:
                tau = varargin{i_varargin + 1};
            else
                err_msg = ['Unrecognized parameter name ' ...
                          varargin{i_varargin}];
                error(err_msg);
            end
            
        end
        
        %After all variables have been entered, verify that A_n, T_n, and
        %sigma all have the same number of elements:
        if ~isequal(numelements(A_n), numelements(T_n)) ||...
           ~isequal(numelements(A_n), numelements(sigma))
       
            err_msg = ['A_n, T_n and sigma all must have the same number' ... 
                       'of elements'];
            error(err_msg);
        end
    end

    
    
    %Convert all varialbes to seconds to match the time vector:
    A_n = A_n * 60; % units = mmol.min * 60 sec/min = mmol.sec
    T_n = T_n * 60; % units = min * 60 sec / min  = sec
    sigma = sigma * 60; %units = min * 60 sec / min = sec
    beta = beta * 1/60; %units = 1/min * 1min / 60 sec = sec ^-1
    s = s * 1/60; %units = 1/min * 1min /60 sec = sec^-1
    tau = tau * 60; %units = min *60 sec / 1 min = sec
    
    %Convert from mmol to mol:
    A_n = A_n * 1/1000; %units = mmol.sec * 1mol/1000 mmol = mol.sec
    alpha = alpha * 1/1000; %units = mmol * 1mol/1000 mmol = mol
    
    %Calculate the contrast agent blood concentration:
    gaussianSum = zeros(size(time)); %initialize to 0
    N = size(A_n, 2); %N is the number of gaussians used to represent the
                      % contrast agent concentration curve
    for iGaussian = 1:N
        gaussianSum = gaussianSum + ...
                      A_n(iGaussian)/(sigma(iGaussian)*sqrt(2*pi)) * ...
                      exp(-(time - T_n(iGaussian)).^2 ./ ...
                          (2*sigma(iGaussian).^2));
    end

    bloodConcentration = gaussianSum + alpha .* exp(-beta * time) ./ ...
                         (1 + exp(-s*(time-tau)));

            
end

% Revision History:
%{
    2022-04-25
        v0.1 released and uploaded to gitHub
%}