function state=bsxops(enable)
% function state=bsxops(enable)
%
% PURPOSE: Enable/disable bsxfun mode for MATLAB standard operators
%
% bsxops(1): enable flexible auto-expansion mode
% bsxops(0): restore "compatible" mode
%
% INQUIRE STATE: state=bsxops() 
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Last update: 21/April/2009

global BSXENABLE

if nargin<1
    if isempty(BSXENABLE)
        enable=1;
    else
        enable=BSXENABLE; % Use the previous state
    end
end

bsxopsroot=fileparts(mfilename('fullpath'));
bsxclasspath=[bsxopsroot filesep 'classops'];
if enable
    addpath(bsxclasspath);
elseif BSXENABLE
    rmpath(bsxclasspath);
end

% Keep track of state
BSXENABLE = enable;

% Return the state (inquiring call)
if nargout>0 || nargin==0
    state=enable;
end

end