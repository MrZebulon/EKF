classdef PoseViewer < HelperPoseViewer
    properties (Hidden)
        thefig
    end
    methods (Access = protected)
        function setupImpl(obj, varargin)
            obj.thefig.Visible = true;
            setupImpl@HelperPoseViewer(obj, varargin{:});
        end
        
        function createAppWindow(obj)
            fig = figure('Name', 'Pose Viewer', ...
                'NumberTitle', 'off', ...
                'DockControls','off', ...
                'Units', 'normalized', ...
                'OuterPosition', [0 0.25 0.5 0.5], ...
                'Visible', 'off', ...
                'HandleVisibility', 'on', ...
                'NextPlot', 'new', ...
                'IntegerHandle', 'off', ...
                'CloseRequestFcn', @(x,~)set(x,'Visible', 'off'));
            u1 = uipanel('Title', 'Pose', 'Parent', fig);
            obj.AppWindow = u1;
            obj.thefig = fig;
        end
    end
    
end

function checkboxcb(src, ~)
[~,f] = gcbo;
f.UserData.(src.String) = src.Value;
end

function pulldowncb(src, ~)
[~,f] = gcbo;
str = strsplit(src.String{src.Value});
f.UserData.(src.Tag) = str2double(str{1});
end
