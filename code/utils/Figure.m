classdef Figure < handle
    % Publication ready figure.
    %   This class helps automating the creation of figures that can be
    %   easily exported with minimal effort.
    %
    % AE 2013-03-05
    
    properties
        handle
        fontSizeScreen = 12
        fontSizePrint = 7;
        tickLength = 3; % pt
        fileName = '';
    end
    
    properties (Dependent)
        size
    end
    
    properties (Access = private)
        sizepx
    end
    
    properties (Constant)
        mmPerInch = 25.4
        mmPerPt = 25.4 / 72
        pxPerInch = get(0, 'ScreenPixelsPerInch')
        isHg2 = datenum(version('-date')) > 735856;
    end
    
    
    methods
        
        function self = Figure(varargin)
            
            % open figure/bring it to focus
            first = double(nargin && isnumeric(varargin{1}));
            self.handle = figure(varargin{1 : first});
            clf
%             set(0,'defaulttextinterpreter','none');
            % set properties passed by name/value
            k = first + 1;
            while k < numel(varargin)
                current = varargin{k};
                if any(strcmp(properties(self), current))
                    self.(current) = varargin{k + 1};
                    varargin(k : k + 1) = [];
                else
                    k = k + 2;
                end
            end
            pos = get(self.handle, 'position');
            self.sizepx = pos(3 : 4);
            
            % set figure defaults/forward overrides passed by name/value
            builtin('set', self.handle, ...
                'Color', 'w', ...
                'DefaultAxesTickDir', 'out', ...
                'DefaultAxesTickLength', [0.025 0.025], ...
                'DefaultAxesBox', 'off', ...
                'DefaultAxesFontName', 'Helvetica', ...
                'DefaultAxesFontSize', self.fontSizeScreen, ...
                'ResizeFcn', @(varargin) self.resizeFcn(), ...
                'DeleteFcn', @(varargin) self.deleteFcn(), ...
                'UserData', self, ...
                varargin{first + 1 : end});
            
            % for HG2 set default font size multipliers to 1
            if self.isHg2
                builtin('set', self.handle, ...
                    'DefaultAxesTitleFontSizeMultiplier', 1, ...
                    'DefaultAxesTitleFontWeight', 'normal', ...
                    'DefaultAxesLabelFontSizeMultiplier', 1, ...
                    'DefaultAxesXColor', 'k', ...
                    'DefaultAxesYColor', 'k', ...
                    'DefaultAxesZColor', 'k', ...
                    'DefaultAxesGridColor', 'k', ...
                    'DefaultAxesMinorGridColor', 'k');
            end
        end
        
        
        function set.size(self, sz)
            scale = self.fontSizeScreen / self.fontSizePrint;
            sz = sz * self.pxPerInch / self.mmPerInch * scale;
            self.sizepx = sz;
            pos = get(self.handle, 'position');
            pos = [pos(1), pos(2) + pos(4) - sz(2), sz];
            set(self.handle, 'position', pos)
        end
        
        
        function sz = get.size(self)
            scale = self.fontSizeScreen / self.fontSizePrint;
            sz = self.sizepx;
            sz = sz / self.pxPerInch * self.mmPerInch / scale;
        end
        
        
        function set(self, varargin)
            set(self.handle, varargin{:})
        end
        
        
        function val = get(self, name)
            val = get(self.handle, name);
        end
        
        
        function save(self, file, varargin)
            % Save figure to file using given style.
            %   fig.save(filename) saves the figure using the given file
            %   name.
           
%             % check for Axis Exponents
%             if self.handle.CurrentAxes.XAxis.Exponent
% %                self.handle.CurrentAxes.XAxis.TickLabelFormat = '%.3f';
% %               text(self.handle.CurrentAxes.XLim(2), self.handle.CurrentAxes.XAxis.Label.Position(2),...
% %                 ['x $10^{\:{' num2str(self.handle.CurrentAxes.XAxis.Exponent) '}}$'], ...
% %                 'VerticalAlignment', 'top','HorizontalAlignment', 'right', ...
% %                 'Interpreter', 'latex')
%             end
%             if self.handle.CurrentAxes.YAxis.Exponent
%               self.handle.CurrentAxes.YAxis.Exponent = double(self.handle.CurrentAxes.YAxis.Exponent);
% %               self.handle.CurrentAxes.YAxis.TickLabelFormat = '%.3f';
% %             text(self.handle.CurrentAxes.YLim(1), self.handle.CurrentAxes.YAxis.Limits(2),...
% %                 ['x$10^{\:{' num2str(self.handle.CurrentAxes.YAxis.Exponent) '}}$'], ...
% %                 'VerticalAlignment', 'top','HorizontalAlignment', 'right', ...
% %                 'Interpreter', 'latex')
%             end
            
            
            if nargin < 2
                file = self.fileName;
            else
                self.fileName = file;
            end
            [~, ~, ext] = fileparts(file);
            if isempty(ext)
                ext = 'eps';
            else
                ext = ext(2 : end);
            end
            
            % define style
            s.Version = 1;
            s.Format = ext;
            s.Preview = 'none';
            s.Units = 'centimeters';
            s.Color = 'rgb';
            s.Background = get(self.handle, 'color');
            s.ScaledFontSize = 'auto';
            s.FontMode = 'fixed';
            s.FontSizeMin = 8;
            s.ScaledLineWidth = 'auto';
            s.LineMode = 'none';
            s.LineWidthMin = 0.5;
            s.FontName = 'Helvetica';
            s.FontWeight = 'auto';
            s.FontAngle = 'auto';
            s.FontEncoding = 'latin1';
            s.PSLevel = 2;
            s.LineStyleMap = 'none';
            s.ApplyStyle = false;
            s.Bounds = 'loose';
            s.LockAxes = 'on';
            s.LockAxesTicks = 'on';
            s.ShowUI = 'off';
            s.SeparateText = 'off';
            for i = 1 : 2 : length(varargin)
                s.(varargin{i}) = varargin{i + 1};
            end
            
            % set renderer and resolution
            switch ext
                case {'eps', 'pdf'}
                    s.Renderer = 'opengl';
                    s.Resolution = 'auto';
                    s.FixedFontSize = self.fontSizePrint;
                    if self.isHg2
                        % Workaround for Matlab R2014b and following:
                        % From R2014b on Matlab cannot save line widths
                        % <1pt in eps files. Therefore, we scale everything
                        % up up a factor of 2 in order to compensate
                        s.FixedFontSize = s.FixedFontSize * 2;
                    end
                case {'png', 'jpg', 'tif'}
                    s.Renderer = 'zbuffer';
                    s.Resolution = '72';
                    s.FixedFontSize = self.fontSizeScreen;
                otherwise
                    error('Unknown format: %s', ext)
            end
            scale = s.FixedFontSize / self.fontSizePrint;
            s.Width = self.size(1) / 10 * scale;
            s.Height = self.size(2) / 10 * scale;
            
            if self.isHg2
                pos = get(self.handle, 'Position');
                scr = get(0, 'ScreenSize');
                if any(pos(3 : 4) * s.FixedFontSize / self.fontSizeScreen > scr(3 : 4))
                    warning('Figure too large. Font sizes will be incorrect due to Matlab bug. Split the figure.')
                end
            end
            
            % save figure
            self.fixFontSize(s.FixedFontSize);
            hgexport(self.handle, file, s);
            self.fixFontSize(self.fontSizeScreen);
        end
        
        
        function setsub(self, varargin)
            % Set property for all subplots
            %   fig.setsub('name1', value1, 'name2', value2, ...) sets the
            %   given properties for all subplots.
            
            set(self.subplots(), varargin{:})
        end
        
        
        function fixticks(self)
            % Fix length of tick marks
            %   fig.fixticks() fixes all tick marks to the same standard
            %   length.
            
            for hdl = self.subplots()
                set(hdl, 'TickDir', 'out');
                pos = get(hdl, 'Position');
                
                sz = pos(3 : 4) .* self.size;
                if strcmp(get(hdl, 'DataAspectRatioMode'), 'manual')
                    ar = get(hdl, 'DataAspectRatio');
                    ar = [diff(get(hdl, 'xlim')) diff(get(hdl, 'ylim'))] ./ ar(1 : 2);
                elseif strcmp(get(hdl, 'PlotBoxAspectRatioMode'), 'manual')
                    ar = get(hdl, 'PlotBoxAspectRatio');
                else
                    ar = sz;
                end
                if ar(1) / ar(2) > sz(1) / sz(2)
                    sz(2) = sz(1) * ar(2) / ar(1);
                else
                    sz(1) = sz(2) * ar(1) / ar(2);
                end
                [~, index] = max(sz);
                ticklen = self.tickLength / sz(index) * self.mmPerPt;
                set(hdl, 'TickLength', ticklen * [1 1]);
            end
        end
        
        
        function cleanup(self)
            % Clean up figure.
            %   fig.cleanup() removes boxes on subplots and adjusts tick
            %   length to proper values.
            
            self.setsub('Box', 'off');
            self.fixticks();
            self.fixFontSize(self.fontSizeScreen);
        end
        
    
        function fixFontSize(self, sz)
            % Set font sizes of legends and colorbars in R2014b
            %   fig.fixFontSize(sz) sets the font size of legends and
            %   colorbars to sz. Applies only to >R2014b (HG2).
            
            if self.isHg2
                hdl = [findobj(self.handle, 'Type', 'Legend'); ...
                       findobj(self.handle, 'Type', 'Colorbar')];
                set(hdl, 'FontSize', sz);
            end
        end
        
    end
    
    
    methods (Access = private)
        
        function resizeFcn(self)
            pos = get(self.handle, 'Position');
            self.sizepx = pos(3 : 4);
        end
        
        function deleteFcn(self)
            delete(self)
        end
        
        function handles = subplots(self)
            handles = findobj(self.handle, 'type', 'axes', '-and', '-not', 'tag', 'legend')';
        end
    end
end
