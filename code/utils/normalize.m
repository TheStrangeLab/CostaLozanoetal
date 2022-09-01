function [N,C,S] = normalize(A,varargin)
%NORMALIZE   Normalize data.
%   N = NORMALIZE(A) normalizes data in A using the 'zscore' method, which 
%   centers the data to have mean 0 and scales it to have standard 
%   deviation 1. NaN values are ignored. If A is a matrix or a table, 
%   NORMALIZE operates on each column separately. If A is an N-D array,  
%   NORMALIZE operates along the first array dimension whose size does not 
%   equal 1.
%
%   N = NORMALIZE(A,DIM) specifies the dimension to operate along.
%
%   N = NORMALIZE(...,METHOD) normalizes using the normalization method
%   METHOD. NaN values are ignored. METHOD can be one of the following:
%
%     'zscore'    - (default) normalizes by centering the data to have mean 
%                   0 and scaling it to have standard deviation 1.
%
%     'norm'      - normalizes by scaling the data to unit length using the 
%                   vector 2-norm.
%
%     'center'    - normalizes by centering the data to have mean 0.
%
%     'scale'     - normalizes by scaling the data by the standard 
%                   deviation.
%
%     'range'     - normalizes by rescaling the range of the data to the 
%                   interval [0,1].
%
%     'medianiqr' - normalizes by centering the data to have median 0 and 
%                   scaling it to have interquartile range 1.
%
%   N = NORMALIZE(A,'zscore',METHODTYPE) and
%   N = NORMALIZE(A,DIM,'zscore',METHODTYPE) normalizes using the 'zscore'
%   method specified by METHODTYPE. METHODTYPE can be one of the following:
%
%     'std'    - (default) centers the data to have mean 0 and standard 
%                deviation 1.
%     'robust' - centers the data to have median 0 and median absolute 
%                deviation 1.
%
%   N = NORMALIZE(A,'norm',p) and N = NORMALIZE(A,DIM,'norm',p) normalizes
%   the data using the vector p-norm. p can be any positive real value or
%   Inf. p is 2 by default.
%
%   N = NORMALIZE(A,'center',METHODTYPE) and
%   N = NORMALIZE(A,DIM,'center',METHODTYPE) normalizes using the 'center'
%   method specified by METHODTYPE. METHODTYPE can be one of the following:
%
%     'mean'   - (default) centers the data to have mean 0.
%     'median' - centers the data to have median 0.
%        C     - centers the data by the numeric array C, which must have a
%                compatible size with input A.
%        T     - centers the data in the corresponding table variables T.
%                The variables in T are matched to the variables in A by 
%                variable name.
%
%   N = NORMALIZE(A,'scale',METHODTYPE) and 
%   N = NORMALIZE(A,DIM,'scale',METHODTYPE) normalizes using the 'scale'
%   method specified by METHODTYPE. METHODTYPE can be one of the following:
%
%      'std'   - (default) scales the data by the standard deviation.
%      'mad'   - scales the data by the median absolute deviation.
%      'iqr'   - scales the data by the interquartile range.
%      'first' - scales the data by the first element.
%        S     - scales the data by the numeric array S, which must have a
%                compatible size with input A.
%        T     - scales the data in the corresponding table variables T.
%                The variables in T are matched to the variables in A by 
%                variable name.
%
%   N = NORMALIZE(A,'range',[a,b]) and NORMALIZE(A,DIM,'range',[a,b]) 
%   normalizes the output range to [a,b]. The default range is [0,1].
%
%   N = NORMALIZE(A,'center',METHODTYPE,'scale',METHODTYPE2) or
%   N = NORMALIZE(A,'scale',METHODTYPE2,'center',METHODTYPE) normalizes
%   using the 'center' method specified by METHODTYPE and 'scale' method
%   specified by METHODTYPE2. See above for METHODTYPE and METHODTYPE2
%   options. You can use this syntax to specify center and scale values C
%   and S from a previously computed normalization:
%      
%       * When C, S, and A are arrays, NORMALIZE uses the values 
%         to compute N = (A - C) ./ S.
%       * When C, S, and A are tables, NORMALIZE uses the values 
%         to compute N.Var = (A.Var - C.Var) ./ S.Var for each 
%         variable in C and S that is also present in A.
%
%   N = NORMALIZE(...,'DataVariables',DV) normalizes the data only in the
%   table variables specified by DV. Other variables in the table not
%   specified by DV pass through to the output without being operated on.
%   The default is all table variables in A. DV must be a table variable
%   name, a cell array of table variable names, a vector of table variable
%   indices, a logical vector, a function handle that returns a logical
%   scalar (such as @isnumeric), or a table vartype subscript. The output
%   table N has the same size as the input table A.
%
%   [N,C,S] = NORMALIZE(...) also returns the center and scale used to do
%   the normalization with any of the above methods.
%
%       * For array input, C and S are arrays where N = (A - C) ./ S.
%       * For table input, C and S are tables containing the centers and 
%         scales for each table variable that was normalized, 
%         N.Var = (A.Var - C.Var) ./ S.Var. The table variable names of C 
%         and S match the corresponding table variables in the input table.
%
%   EXAMPLE: Compute the z-scores of two data vectors in order to compare 
%   the data on one plot
%       A(:,1) = 3*rand(5,1);
%       A(:,2) = 2000*rand(5,1);
%       N = normalize(A);
%       plot(N);
%
%   EXAMPLE: Scale each column of a matrix to unit length using the vector
%   2-norm
%       A = rand(10);
%       N = normalize(A,'norm');
%
%   See also RESCALE, SMOOTHDATA, FILLOUTLIERS, FILLMISSING, VECNORM.
 
%   Copyright 2017-2020 The MathWorks, Inc.

[dim,method,methodType,dataVars,AisTablular,method2,methodType2] = parseInputs(A,varargin{:});

if ~AisTablular
    % Only check first method since "center" and "scale" have same requirements
    checkSupportedArray(A,method,methodType,false);
    [N,C,S] = normalizeArray(A,method,methodType,dim,method2,methodType2);
else
    N = A;
    if isempty(dataVars)
        C = table.empty(min(size(A,1),1),0);
        S = table.empty(min(size(A,1),1),0);
    else
        C = table();
        S = table();
    end
    varNames = A.Properties.VariableNames;
    for vj = dataVars
        namevj = varNames{vj};
        Avj = A.(namevj);
        % Only check first method since "center" and "scale" have same requirements
        checkSupportedArray(Avj,method,methodType,true);
        if isa(methodType,'tabular')
            mT = methodType.(namevj);
        else
            mT = methodType;
        end
        if isa(methodType2,'tabular')
            mT2 = methodType2.(namevj);
        else
            mT2 = methodType2;
        end
        [N.(namevj),C.(namevj),S.(namevj)] = normalizeArray(A.(namevj),method,mT,dim,method2,mT2);
    end
end
end

%--------------------------------------------------------------------------
function [N,C,S] = normalizeArray(A,method,methodType,dim,method2,methodType2)
% Normalization for arrays - always omit NaNs

if ~isempty(method2)
    if method == "center"
        % Find center
        C = centerArray(A,methodType,dim);
        % Find scale
        S = scaleArray(A,methodType2,dim);
    else
        % Find center
        C = centerArray(A,methodType2,dim);
        % Find scale
        S = scaleArray(A,methodType,dim);
    end
    N = (A - C) ./ S;
elseif isequal("zscore", method)
    if isequal("std",methodType)
        C = mean(A,dim,'omitnan');
        S = std(A,0,dim,'omitnan');
    else % "robust"
        C = median(A,dim,'omitnan');
        S = median(abs(A - C),dim,'omitnan');
    end
    N = (A - C) ./ S;
elseif isequal("norm", method)
    % In order to omit NaNs in this case fill NaNs with 0 to compute norms
    fillA = A;
    fillA(isnan(fillA)) = 0;
    S = vecnorm(fillA,methodType,dim);
    C = zeros(1,'like',S);
    N = A./S;
elseif isequal("center", method)
    C = centerArray(A,methodType,dim);
    S = ones(1,'like',C);
    N = A - C;
elseif isequal("scale", method)
    S = scaleArray(A,methodType,dim);
    C = zeros(1,'like',S);
    N = A ./ S;
elseif isequal("range", method)
    minA = min(A,[],dim);
    maxA = max(A,[],dim);
    if ~isfloat(A)
        minA = double(minA);
        maxA = double(maxA);
    end
    if ~isfloat(methodType)
        methodType = double(methodType);
    end
    a = methodType(1);
    b = methodType(2);
    [N,C,S] = matlab.internal.math.rescaleKernel(A,a,b,minA,maxA);
elseif isequal("medianiqr", method)
    C = median(A,dim,'omitnan');
    S = datafuniqr(A,[25 75],dim,true);
    N = (A - C) ./ S;
end
end

%--------------------------------------------------------------------------
function C = centerArray(A,methodType,dim)
if isequal("mean",methodType)
    C = mean(A,dim,'omitnan');
elseif isequal("median",methodType)
    C = median(A,dim,'omitnan');
else % numeric
    C = methodType;
end
end

%--------------------------------------------------------------------------
function S = scaleArray(A,methodType,dim)
if isequal("std",methodType)
    S = std(A,0,dim,'omitnan');
elseif isequal("mad",methodType)
    S = median(abs(A - median(A,dim,'omitnan')),dim,'omitnan');
elseif isequal("first",methodType)
    if isempty(A)
        S = 1;
    else
        ind = repmat({':'},ndims(A),1);
        ind{dim} = 1;
        S = A(ind{:});
    end
elseif isequal("iqr",methodType)
    S = datafuniqr(A,[25 75],dim,true);
else % numeric
    S = methodType;
end
end

%--------------------------------------------------------------------------
function checkSupportedArray(A,method,methodType,AisTabular)
% Parse input A
if isequal("range",method)
    if (~(isnumeric(A) || islogical(A)) || ~isreal(A))
        if AisTabular
            error(message('MATLAB:normalize:UnsupportedTableVariableRange'));
        else
            error(message('MATLAB:normalize:InvalidFirstInputRange'));
        end
    end
elseif isequal("zscore",method) && isequal("robust",methodType)
    if ~isfloat(A) || ~isreal(A)
        if AisTabular
            error(message('MATLAB:normalize:UnsupportedTableVariableRobust'));
        else
            error(message('MATLAB:normalize:InvalidFirstInputRobust'));
        end
    end
else
    if ~isfloat(A)
        if AisTabular
            error(message('MATLAB:normalize:UnsupportedTableVariable'));
        else
            error(message('MATLAB:normalize:InvalidFirstInput'));
        end
    end
end
end
%--------------------------------------------------------------------------
function [dim,method,methodType,dataVars,AisTabular,method2,methodType2] = parseInputs(A,varargin)
% Parse NORMALIZE inputs
AisTabular = matlab.internal.datatypes.istabular(A);

% Set defaults
method = "zscore";
methodType = "std";
method2 = [];
methodType2 = [];
dataVarsProvided = false;
if ~AisTabular
    dim = find(size(A) ~= 1,1); % default to first non-singleton dimension
    if isempty(dim)
        dim = 2; % dim = 2 for scalar and empty A
    end
    dataVars = []; % not supported for arrays
else
    dim = 1; % Fill each table variable separately
    dataVars = 1:width(A);
end

% NORMALIZE(A,DIM)
% NORMALIZE(A,DIM,METHOD)
% NORMALIZE(A,DIM,METHOD,TYPE)
% NORMALIZE(A,METHOD)
% NORMALIZE(A,METHOD,TYPE)
% NORMALIZE(A,'center','scale')
% NORMALIZE(A,'center',TYPE,'scale')
% NORMALIZE(A,'center','scale',TYPE2)
% NORMALIZE(A,'center',TYPE,'scale',TYPE2)
% NORMALIZE(A,...,'DataVariables',DV) % Not for DIM syntaxes

if nargin > 1
    indStart = 1;
    % Parse dimension - errors for invalid dim
    dimProvided = false;
    methodProvided = false;
    if isnumeric(varargin{indStart}) || islogical(varargin{indStart})
        if AisTabular
            error(message('MATLAB:normalize:TableDIM'));
        end
        dim = varargin{indStart};
        if ~isscalar(dim) || ~isreal(dim) || dim < 1 || ~(fix(dim) == dim) || ~ isfinite(dim)
            error(message('MATLAB:normalize:InvalidDIM'));
        end
        indStart = indStart + 1;
        dimProvided = true;
    end
    if indStart < nargin
        % Parse method - does not error for invalid method
        validMethods = ["zscore","norm","center","scale","range","medianiqr"];
        if checkCharString(varargin{indStart})
            indMethod = startsWith(validMethods, varargin{indStart}, 'IgnoreCase', true);
            if nnz(indMethod) == 1
                method = validMethods(indMethod);
                defaultTypes = {"std",2,"mean","std",[0,1],"none"};
                methodType = defaultTypes{indMethod};
                indStart = indStart + 1;
                methodProvided = true;
                % Parse type - does not error for invalid character/string type
                if indStart < nargin
                    [methodType,indStart] = parseType(varargin{indStart},method,methodType,indStart,AisTabular);
                end
            end
        end
        
        if (indStart < nargin) && (method == "center" || method == "scale")
            % Check if user input both center and scale methods
            if method == "center"
                validSecondMethod = "scale";
                methodType2 = "std";
            else
                validSecondMethod = "center";
                methodType2 = "mean";
            end
            if checkCharString(varargin{indStart})
                indMethod = startsWith(validSecondMethod, varargin{indStart}, 'IgnoreCase', true);
                if nnz(indMethod) == 1
                    method2 = validSecondMethod;
                    indStart = indStart + 1;
                    % Parse type - does not error for invalid character/string type
                    if indStart < nargin
                        [methodType2,indStart] = parseType(varargin{indStart},method2,methodType2,indStart,AisTabular);
                    end
                end
            end
        end

        % Parse name-value pairs
        if rem(nargin-indStart,2) == 0
            for j = indStart:2:length(varargin)
                name = varargin{j};
                if ~AisTabular
                    if any(strncmpi(validMethods,name,1))
                        if methodProvided
                            error(message('MATLAB:normalize:InvalidDoubleMethod'));
                        else
                            error(message('MATLAB:normalize:InvalidMethod'));
                        end
                    else   
                        error(message('MATLAB:normalize:DataVariablesArray'));
                    end
                elseif ~checkCharString(name)
                    error(message('MATLAB:normalize:ParseFlags'));
                elseif startsWith("DataVariables", name, 'IgnoreCase', true)
                    dataVars = matlab.internal.math.checkDataVariables(A, varargin{j+1}, 'normalize');
                    dataVarsProvided = true;
                elseif any(strncmpi(validMethods,name,1))
                    if methodProvided
                        error(message('MATLAB:normalize:InvalidDoubleMethod'));
                    elseif dataVarsProvided
                        error(message('MATLAB:normalize:MethodAfterOptions'));
                    else
                        error(message('MATLAB:normalize:InvalidMethod'));
                    end
                else
                    error(message('MATLAB:normalize:ParseFlags'));
                end
            end
        elseif (nargin < 3) || (dimProvided && nargin < 4)
            error(message('MATLAB:normalize:InvalidMethod'));
        else
            if methodProvided && checkCharString(varargin{indStart}) && ...
                    any(startsWith(validMethods, varargin{indStart}, 'IgnoreCase', true))
                error(message('MATLAB:normalize:InvalidDoubleMethod'));
            elseif ~AisTabular
                error(message('MATLAB:normalize:IncorrectNumInputsArray'));
            else
                error(message('MATLAB:normalize:KeyWithoutValue'));
            end
        end
    end
    methodTypeIsTabular = matlab.internal.datatypes.istabular(methodType);
    methodType2IsTabular = matlab.internal.datatypes.istabular(methodType2);
    % If methodType and/or methodType2 are tables we need to check or
    % update datavars
    if AisTabular && (methodTypeIsTabular || methodType2IsTabular)
        if dataVarsProvided
            % When 'DataVariables' is provided we need to check that those
            % variables are in methodType and/or methodType2
            varNames = A.Properties.VariableNames(dataVars);
            dataVarsInMethodType = true;
            dataVarsInMethodType2 = true;
            if methodTypeIsTabular
                dataVarsInMethodType = all(ismember(varNames,methodType.Properties.VariableNames));
            end
            if methodType2IsTabular
                dataVarsInMethodType2 = all(ismember(varNames,methodType2.Properties.VariableNames));
            end
            if ~(dataVarsInMethodType && dataVarsInMethodType2)
                error(message('MATLAB:normalize:InvalidCenterScaleTypeWithDataVars'));
            end
        else
            % When 'DataVariables' is not provided we need to use the
            % variables in methodType and/or methodType2 as datavars. In
            % this case, methodType and methodType2 must have the same
            % table variable names.
            if methodTypeIsTabular && ~methodType2IsTabular
                vars = methodType.Properties.VariableNames;
            elseif ~methodTypeIsTabular && methodType2IsTabular
                vars = methodType2.Properties.VariableNames;
            elseif methodTypeIsTabular && methodType2IsTabular
                vars = methodType.Properties.VariableNames;
                vars2 = methodType2.Properties.VariableNames;
                if ~isempty(setxor(vars,vars2))
                    error(message('MATLAB:normalize:InvalidCenterScaleType'));
                end
            end
            try
                dataVars = matlab.internal.math.checkDataVariables(A, vars, 'normalize');
            catch
                error(message('MATLAB:normalize:InvalidCenterScaleTypeFirstInput'));
            end
        end
    end
end
end
%--------------------------------------------------------------------------
function [methodType,indStart] = parseType(input,method,methodType,indStart,AisTabular)
% Parse Method Type 
if method == "zscore"
    if checkCharString(input)
        validZscoreType = ["std","robust"];
        indZscoreType = startsWith(validZscoreType, input, 'IgnoreCase', true);
        if nnz(indZscoreType) == 1
            methodType = validZscoreType(indZscoreType);
            indStart = indStart + 1;
        elseif ~AisTabular && ~startsWith("DataVariables", input, 'IgnoreCase', true)
            error(message('MATLAB:normalize:InvalidZscoreType'));
        end
    else
        error(message('MATLAB:normalize:InvalidZscoreType'));
    end
elseif method == "norm"
    if isnumeric(input) || islogical(input)
        if ~isscalar(input) || (input <= 0 || ~isreal(input)) || islogical(input)
            error(message('MATLAB:normalize:InvalidNormType'));
        end
        methodType = input;
        indStart = indStart + 1;
    elseif checkCharString(input)
        if strcmpi("inf",input)
            methodType = "inf";
            indStart = indStart + 1;
        elseif ~AisTabular && ~startsWith("DataVariables", input, 'IgnoreCase', true)
            error(message('MATLAB:normalize:InvalidNormType'));    
        end
    else
        error(message('MATLAB:normalize:InvalidNormType'));
    end
elseif method == "center"
    if isnumeric(input) || islogical(input)
        if islogical(input)
            error(message('MATLAB:normalize:InvalidCenterType'));
        end
        methodType = input;
        indStart = indStart + 1;
        if ~isfloat(methodType)
            methodType = double(methodType);
        end
    elseif checkCharString(input)
        validCenterType = ["mean","median"];
        indCenterType = startsWith(validCenterType, input, 'IgnoreCase', true);
        if nnz(indCenterType) == 1
            methodType = validCenterType(indCenterType);
            indStart = indStart + 1;
        elseif ~AisTabular && ~any(startsWith(["DataVariables","scale"], input, 'IgnoreCase', true))
            error(message('MATLAB:normalize:InvalidCenterType'));
        end
    elseif AisTabular && isa(input,'tabular')
        methodType = input;
        indStart = indStart + 1;
    else
        error(message('MATLAB:normalize:InvalidCenterType'));
    end
elseif method == "scale"
    if isnumeric(input) || islogical(input)
        if islogical(input)
            error(message('MATLAB:normalize:InvalidScaleType'));
        end
        methodType = input;
        indStart = indStart + 1;
        if ~isfloat(methodType)
            methodType = double(methodType);
        end 
    elseif checkCharString(input)
        validScaleType = ["std","mad","first","iqr"];
        indScaleType = startsWith(validScaleType, input, 'IgnoreCase', true);
        if nnz(indScaleType) == 1
            methodType = validScaleType(indScaleType);
            indStart = indStart + 1;
        elseif ~AisTabular && ~any(startsWith(["DataVariables","center"], input, 'IgnoreCase', true))
            error(message('MATLAB:normalize:InvalidScaleType'));
        end
    elseif AisTabular && isa(input,'tabular')
        methodType = input;
        indStart = indStart + 1;
    else
        error(message('MATLAB:normalize:InvalidScaleType'));
    end
elseif method == "range"
    if isnumeric(input) || islogical(input)
        if ~isvector(input) || length(input) ~= 2 || ~isreal(input) ||...
                islogical(input) || diff(input) < 0
            error(message('MATLAB:normalize:InvalidRangeType'));
        end
        methodType = input;
        indStart = indStart + 1;
    elseif checkCharString(input)
        if ~AisTabular && ~startsWith("DataVariables", input, 'IgnoreCase', true)
            error(message('MATLAB:normalize:InvalidRangeType'));
        end
    else
        error(message('MATLAB:normalize:InvalidRangeType'));
    end
end
end

%--------------------------------------------------------------------------
function flag = checkCharString(inputName)
flag = (ischar(inputName) && isrow(inputName)) || (isstring(inputName) && isscalar(inputName) ...
    && strlength(inputName) ~= 0);
end

