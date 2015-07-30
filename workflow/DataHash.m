function Hash = DataHash(Data, Opt)
% DATAHASH - Checksum for Matlab array of any type
% This function creates a hash value for an input of any type. The type and
% dimensions of the input are considered as default, such that UINT8([0,0]) and
% UINT16(0) have different hash values. Nested STRUCTs and CELLs are parsed
% recursively.
%
% Hash = DataHash(Data, Opt)
% INPUT:
%   Data: Array of these built-in types:
%           (U)INT8/16/32/64, SINGLE, DOUBLE, (real or complex)
%           CHAR, LOGICAL, CELL (nested), STRUCT (scalar or array, nested),
%           function_handle.
%   Opt:  Struct to specify the hashing algorithm and the output format.
%         Opt and all its fields are optional.
%         Opt.Method: String, known methods for Java 1.6 (Matlab 2009a):
%              'SHA-1', 'SHA-256', 'SHA-384', 'SHA-512', 'MD2', 'MD5'.
%            Known methods for Java 1.3 (Matlab 6.5):
%              'MD5', 'SHA-1'.
%            Default: 'MD5'.
%         Opt.Format: String specifying the output format:
%            'hex', 'HEX':      Lower/uppercase hexadecimal string.
%            'double', 'uint8': Numerical vector.
%            'base64':          Base64 encoded string, only printable ASCII
%                               characters, shorter than 'hex', no padding.
%            Default: 'hex'.
%         Opt.Input: Type of the input as string, not case-sensitive:
%            'array': The contents, type and size of the input [Data] are
%                     considered  for the creation of the hash. Nested CELLs
%                     and STRUCT arrays are parsed recursively. Empty arrays of
%                     different type reply different hashs.
%            'file':  [Data] is treated as file name and the hash is calculated
%                     for the files contents.
%            'bin':   [Data] is a numerical, LOGICAL or CHAR array. Only the
%                     binary contents of the array is considered, such that
%                     e.g. empty arrays of different type reply the same hash.
%            Default: 'array'.
%
% OUTPUT:
%   Hash: String, DOUBLE or UINT8 vector. The length depends on the hashing
%         method.
%
% EXAMPLES:
% % Default: MD5, hex:
%   DataHash([])                % 7de5637fd217d0e44e0082f4d79b3e73
% % MD5, Base64:
%   Opt.Format = 'base64';
%   Opt.Method = 'MD5';
%   DataHash(int32(1:10), Opt)  % bKdecqzUpOrL4oxzk+cfyg
% % SHA-1, Base64:
%   S.a = uint8([]);
%   S.b = {{1:10}, struct('q', uint64(415))};
%   Opt.Method = 'SHA-1';
%   DataHash(S, Opt)            % ZMe4eUAp0G9TDrvSW0/Qc0gQ9/A
% % SHA-1 of binary values:
%   Opt.Method = 'SHA-1';
%   Opt.Input  = 'bin';
%   DataHash(1:8, Opt)          % 826cf9d3a5d74bbe415e97d4cecf03f445f69225
%
% NOTES:
%   Function handles and user-defined objects cannot be converted uniquely:
%   - The subfunction ConvertFuncHandle uses the built-in function FUNCTIONS,
%     but the replied struct can depend on the Matlab version.
%   - It is tried to convert objects to UINT8 streams in the subfunction
%     ConvertObject. A conversion by STRUCT() might be more appropriate.
%   Adjust these subfunctions on demand.
%
%   MATLAB CHARs have 16 bits! In consequence the string 'hello' is treated as
%   UINT16('hello') for the binary input method. Use this to get the hash of an
%   ASCII-string (Result as defined in RFC 1321!):
%     Opt.Method = 'MD5'; Opt.Input = 'bin';
%     DataHash(uint8('abc'), Opt);    % '900150983cd24fb0d6963f7d28e17f72'
%
%   DataHash uses James Tursa's smart and fast TYPECASTX, if it is installed:
%     http://www.mathworks.com/matlabcentral/fileexchange/17476
%   As fallback the built-in TYPECAST is used automatically, but for large
%   inputs this can be more than 100 times slower.
%
% Tested: Matlab 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2011-2015 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also: TYPECAST, CAST.
% FEX:
% Michael Kleder, "Compute Hash", no structs and cells:
%   http://www.mathworks.com/matlabcentral/fileexchange/8944
% Tim, "Serialize/Deserialize", converts structs and cells to a byte stream:
%   http://www.mathworks.com/matlabcentral/fileexchange/29457
% Jan Simon, "CalcMD5", MD5 only, much faster C-mex:
%   http://www.mathworks.com/matlabcentral/fileexchange/25921

% $JRev: R-v V:022 Sum:68JCxGh2/Q0N Date:30-Mar-2015 01:35:37 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLFile\DataHash.m $
% History:
% 001: 01-May-2011 21:52, First version.
% 007: 10-Jun-2011 10:38, [Opt.Input], binary data, complex values considered.
% 011: 26-May-2012 15:57, Fixed: Failed for binary input and empty data.
% 014: 04-Nov-2012 11:37, Consider Mex-, MDL- and P-files also.
%      Thanks to David (author 243360), who found this bug.
%      Jan Achterhold (author 267816) suggested to consider Java objects.
% 016: 01-Feb-2015 20:53, Java heap space exhausted for large files.
%      Now files are process in chunks to save memory.
% 017: 15-Feb-2015 19:40, Collsions: Same hash for different data.
%      Examples: zeros(1,1) and zeros(1,1,0)
%                complex(0) and zeros(1,1,0,0)
%      Now the number of dimensions is included, to avoid this.
% 022: 30-Mar-2015 00:04, Bugfix: Failed for strings and [] without TYPECASTX.
%      Ross found these 2 bugs, which occur when TYPECASTX is not installed.
%      If you need the base64 format padded with '=' characters, adjust
%      fBase64_enc as you like.

% OPEN BUGS:
% Nath wrote:
% function handle refering to struct containing the function will create
% infinite loop. Is there any workaround ?
% Example:
%   d= dynamicprops();
%   addprop(d,'f');
%   d.f= @(varargin) struct2cell(d);
%   DataHash(d.f) % infinite loop

% Main function: ===============================================================
% typecastx creates a shared data copy instead of the deep copy as Matlab's
% TYPECAST - for a [1000x1000] DOUBLE array this is 100 times faster!
persistent usetypecastx
if isempty(usetypecastx)
   % Java is needed:
   if ~usejava('jvm')
      Error_L('NoJava', 'Java is required.');
   end

   usetypecastx = ~isempty(which('typecastx'));  % Run the slow WHICH once only
end

% Default options: -------------------------------------------------------------
Method    = 'MD5';
OutFormat = 'hex';
isFile    = false;
isBin     = false;

% Check number and type of inputs: ---------------------------------------------
nArg = nargin;
if nArg == 2
   if isa(Opt, 'struct') == 0   % Bad type of 2nd input:
      Error_L('BadInput2', '2nd input [Opt] must be a struct.');
   end
   
   % Specify hash algorithm:
   if isfield(Opt, 'Method')
      Method = upper(Opt.Method);
   end
   
   % Specify output format:
   if isfield(Opt, 'Format')
      OutFormat = Opt.Format;
   end
   
   % Check if the Input type is specified - default: 'array':
   if isfield(Opt, 'Input')
      if strcmpi(Opt.Input, 'File')
         isFile = true;
         if ischar(Data) == 0
            Error_L('CannotOpen', '1st input FileName must be a string');
         end
         
      elseif strncmpi(Opt.Input, 'bin', 3)  % Accept 'binary' also
         isBin = true;
         if (isnumeric(Data) || ischar(Data) || islogical(Data)) == 0 || ...
               issparse(Data)
            Error_L('BadDataType', ...
               '1st input must be numeric, CHAR or LOGICAL for binary input.');
         end
      end
   end
   
elseif nArg ~= 1  % Bad number of arguments:
   Error_L('BadNInput', '1 or 2 inputs required.');
end

% Create the engine: -----------------------------------------------------------
try
   Engine = java.security.MessageDigest.getInstance(Method);
catch
   Error_L('BadInput2', 'Invalid algorithm: [%s].', Method);
end

% Create the hash value: -------------------------------------------------------
if isFile
   % Check existence of file:
   Found = FileExist(Data);
   if ~Found
      Error_L('FileNotFound', 'File not found: %s.', Data);
   end
   
   % Open the file:
   FID = fopen(Data, 'r');
   if FID < 0
      Error_L('CannotOpen', 'Cannot open file: %s.', Data);
   end
   
   % Read file in chunks to save memory and Java heap space:
   Chunk         = 1e6;
   [Data, Count] = fread(FID, Chunk, '*uint8');
   Engine.update(Data);
   while Count == Chunk
      [Data, Count] = fread(FID, Chunk, '*uint8');
      Engine.update(Data);
   end
   fclose(FID);
   
   % Calculate the hash:
   if usetypecastx
      Hash = typecastx(Engine.digest, 'uint8');
   else
      Hash = typecast(Engine.digest, 'uint8');
   end

elseif isBin             % Contents of an elementary array, type tested already:
   if isempty(Data)      % Nothing to do, Engine.update fails for empty input!
      if usetypecastx    % Bugfix: Consider missing typecastx
         Hash = typecastx(Engine.digest, 'uint8');
      else
         Hash = typecast(Engine.digest, 'uint8');
      end
      
   elseif usetypecastx   % Faster typecastx:
      if isreal(Data)
         Engine.update(typecastx(Data(:), 'uint8'));
      else
         Engine.update(typecastx(real(Data(:)), 'uint8'));
         Engine.update(typecastx(imag(Data(:)), 'uint8'));
      end
      Hash = typecastx(Engine.digest, 'uint8');
      
   else                  % Matlab's TYPECAST is less elegant:
      if isnumeric(Data)
         if isreal(Data)
            Engine.update(typecast(Data(:), 'uint8'));
         else
            Engine.update(typecast(real(Data(:)), 'uint8'));
            Engine.update(typecast(imag(Data(:)), 'uint8'));
         end
      elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
         Engine.update(typecast(uint8(Data(:)), 'uint8'));
      elseif ischar(Data)                  % TYPECAST cannot handle CHAR
         Engine.update(typecast(uint16(Data(:)), 'uint8'));
         % Bugfix: Line removed
      end
      Hash = typecast(Engine.digest, 'uint8');
   end
   
elseif usetypecastx  % Faster typecastx:
   Engine = CoreHash_(Data, Engine);
   Hash   = typecastx(Engine.digest, 'uint8');
   
else                 % Slower built-in TYPECAST:
   Engine = CoreHash(Data, Engine);
   Hash   = typecast(Engine.digest, 'uint8');
end

% Convert hash specific output format: -----------------------------------------
switch OutFormat
   case 'hex'
      Hash = sprintf('%.2x', double(Hash));
   case 'HEX'
      Hash = sprintf('%.2X', double(Hash));
   case 'double'
      Hash = double(reshape(Hash, 1, []));
   case 'uint8'
      Hash = reshape(Hash, 1, []);
   case 'base64'
      Hash = fBase64_enc(double(Hash));
   otherwise
      Error_L('BadOutFormat', ...
         '[Opt.Format] must be: HEX, hex, uint8, double, base64.');
end

% return;

% ******************************************************************************
function Engine = CoreHash_(Data, Engine)
% This method uses the faster typecastx version.

% Consider the type and dimensions of the array to distinguish arrays with the
% same data, but different shape: [0 x 0] and [0 x 1], [1,2] and [1;2],
% DOUBLE(0) and SINGLE([0,0]):
% <  v016: [class, size, data]. BUG! 0 and zeros(1,1,0) had the same hash!
% >= v016: [class, ndim, size, data]
Engine.update([uint8(class(Data)), ...
              typecastx([ndims(Data), size(Data)], 'uint8')]);

% Special treatment for sparse arrays - store the indices at first and the
% values afterwards:
if issparse(Data)
   % Replace Data by vector of non-zero elements:
   [i1, i2, Data] = find(Data);
   Engine.update(typecastx(i1, 'uint8'));
   Engine.update(typecastx(i2, 'uint8'));
end

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash_(F, Engine);    % Consider the fieldnames
   
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash_(Data(iS).(F{iField}), Engine);
      end
   end
   
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash_(Data{iS}, Engine);
   end
      
elseif isnumeric(Data) || islogical(Data) || ischar(Data)
   if isempty(Data) == 0
      if isreal(Data)                % TRUE for LOGICAL and CHAR also:
         Engine.update(typecastx(Data(:), 'uint8'));
      else                           % Complex input:
         Engine.update(typecastx(real(Data(:)), 'uint8'));
         Engine.update(typecastx(imag(Data(:)), 'uint8'));
      end
   end
   
elseif isa(Data, 'function_handle')
   Engine = CoreHash_(ConvertFuncHandle(Data), Engine);
   
elseif (isobject(Data) || isjava(Data)) && ismethod(Data, 'hashCode')
   Engine = CoreHash_(Data.hashCode, Engine);
   
else  % Most likely this is a user-defined object:
   try
      Engine = CoreHash_(ConvertObject(Data), Engine);
   catch
      warning(['JSimon:', mfilename, ':BadDataType'], ...
         ['Type of variable not considered: ', class(Data)]);
   end
end

% return;

% ******************************************************************************
function Engine = CoreHash(Data, Engine)
% This methods uses the slower TYPECAST of Matlab
% See CoreHash_ for comments.

Engine.update([uint8(class(Data)), ...
              typecast([ndims(Data), size(Data)], 'uint8')]);

if isstruct(Data)                    % Hash for all array elements and fields:
   F      = sort(fieldnames(Data));  % Ignore order of fields
   Engine = CoreHash(F, Engine);     % Catch the fieldnames
   for iS = 1:numel(Data)            % Loop over elements of struct array
      for iField = 1:length(F)       % Loop over fields
         Engine = CoreHash(Data(iS).(F{iField}), Engine);
      end
   end
elseif iscell(Data)                  % Get hash for all cell elements:
   for iS = 1:numel(Data)
      Engine = CoreHash(Data{iS}, Engine);
   end
elseif isempty(Data)
elseif isnumeric(Data)
   if isreal(Data)
      Engine.update(typecast(Data(:), 'uint8'));
   else
      Engine.update(typecast(real(Data(:)), 'uint8'));
      Engine.update(typecast(imag(Data(:)), 'uint8'));
   end
elseif islogical(Data)               % TYPECAST cannot handle LOGICAL
   Engine.update(typecast(uint8(Data(:)), 'uint8'));
elseif ischar(Data)                  % TYPECAST cannot handle CHAR
   Engine.update(typecast(uint16(Data(:)), 'uint8'));
elseif isa(Data, 'function_handle')
   Engine = CoreHash(ConvertFuncHandle(Data), Engine);
elseif (isobject(Data) || isjava(Data)) && ismethod(Data, 'hashCode')
   Engine = CoreHash(Data.hashCode, Engine);
else  % Most likely a user-defined object:
   try
      Engine = CoreHash(ConvertObject(Data), Engine);
   catch
      warning(['JSimon:', mfilename, ':BadDataType'], ...
         ['Type of variable not considered: ', class(Data)]);
   end
end

% return;

% ******************************************************************************
function FuncKey = ConvertFuncHandle(FuncH)
%   The subfunction ConvertFuncHandle converts function_handles to a struct
%   using the Matlab function FUNCTIONS. The output of this function changes
%   with the Matlab version, such that DataHash(@sin) replies different hashes
%   under Matlab 6.5 and 2009a.
%   An alternative is using the function name and name of the file for
%   function_handles, but this is not unique for nested or anonymous functions.
%   If the MATLABROOT is removed from the file's path, at least the hash of
%   Matlab's toolbox functions is (usually!) not influenced by the version.
%   Finally I'm in doubt if there is a unique method to hash function handles.
%   Please adjust the subfunction ConvertFuncHandles to your needs.

% The Matlab version influences the conversion by FUNCTIONS:
% 1. The format of the struct replied FUNCTIONS is not fixed,
% 2. The full paths of toolbox function e.g. for @mean differ.
FuncKey = functions(FuncH);

% ALTERNATIVE: Use name and path. The <matlabroot> part of the toolbox functions
% is replaced such that the hash for @mean does not depend on the Matlab
% version.
% Drawbacks: Anonymous functions, nested functions...
% funcStruct = functions(FuncH);
% funcfile   = strrep(funcStruct.file, matlabroot, '<MATLAB>');
% FuncKey    = uint8([funcStruct.function, ' ', funcfile]);

% Finally I'm afraid there is no unique method to get a hash for a function
% handle. Please adjust this conversion to your needs.

% return;

% ******************************************************************************
function DataBin = ConvertObject(DataObj)
% Convert a user-defined object to a binary stream. There cannot be a unique
% solution, so this part is left for the user...

% Perhaps a direct conversion is implemented:
DataBin = uint8(DataObj);

% Or perhaps this is better:
% DataBin = struct(DataObj);

% return;

% ******************************************************************************
function Out = fBase64_enc(In)
% Encode numeric vector of UINT8 values to base64 string.
% The intention of this is to create a shorter hash than the HEX format.
% Therefore a padding with '=' characters is omitted on purpose.

Pool = [65:90, 97:122, 48:57, 43, 47];  % [0:9, a:z, A:Z, +, /]
v8   = [128; 64; 32; 16; 8; 4; 2; 1];
v6   = [32, 16, 8, 4, 2, 1];

In  = reshape(In, 1, []);
X   = rem(floor(In(ones(8, 1), :) ./ v8(:, ones(length(In), 1))), 2);
Y   = reshape([X(:); zeros(6 - rem(numel(X), 6), 1)], 6, []);
Out = char(Pool(1 + v6 * Y));

% return;

% ******************************************************************************
function Ex = FileExist(FileName)
% A more reliable version of EXIST(FileName, 'file'):
dirFile = dir(FileName);
if length(dirFile) == 1
   Ex = ~(dirFile.isdir);
else
   Ex = false;
end

% return;

% ******************************************************************************
function Error_L(ID, varargin)

error(['JSimon:', mfilename, ':', ID], ['*** %s: ', varargin{1}], ...
   mfilename, varargin{2:nargin - 1});

% return;
