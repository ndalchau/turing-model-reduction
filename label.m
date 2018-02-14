function ha = label(str,location)
%LABEL Simple construction of text annotations
%   [HA] = LABEL(STRING,POSITION) where POSITION = [left bottom width height]
%   places the text in STRING at the specified POSITION.

% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

ha = annotation('textbox',location,'String',str,'EdgeColor','none',...
	'FontSize',12,'VerticalAlignment','middle','FontName','Arial');

return