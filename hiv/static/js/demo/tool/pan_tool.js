(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "common/logging", "./tool", "./event_generators"], function(_, Backbone, Logging, Tool, EventGenerators) {
    var PanTool, PanToolView, PanTools, TwoPointEventGenerator, logger, _ref, _ref1, _ref2;
    TwoPointEventGenerator = EventGenerators.TwoPointEventGenerator;
    logger = Logging.logger;
    window.render_count = 0;
    PanToolView = (function(_super) {
      __extends(PanToolView, _super);

      function PanToolView() {
        _ref = PanToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      PanToolView.prototype.initialize = function(options) {
        var dims;
        PanToolView.__super__.initialize.call(this, options);
        dims = this.mget('dimensions');
        if (dims.length === 0) {
          return logger.warn("pan tool given empty dimensions");
        } else if (dims.length === 1) {
          if (dims[0] === 'width') {
            return this.evgen_options.buttonText = "Pan (x-axis)";
          } else if (dims[0] === 'height') {
            return this.evgen_options.buttonText = "Pan (y-axis)";
          } else {
            return logger.warn("pan tool given unrecognized dimensions: " + dims);
          }
        } else if (dims.length === 2) {
          if (dims.indexOf('width') < 0 || dims.indexOf('height') < 0) {
            return logger.warn("pan tool given unrecognized dimensions: " + dims);
          }
        } else {
          return logger.warn("pan tool given more than two dimensions: " + dims);
        }
      };

      PanToolView.prototype.bind_bokeh_events = function() {
        return PanToolView.__super__.bind_bokeh_events.call(this);
      };

      PanToolView.prototype.eventGeneratorClass = TwoPointEventGenerator;

      PanToolView.prototype.toolType = "PanTool";

      PanToolView.prototype.evgen_options = {
        keyName: null,
        buttonText: "Pan",
        buttonHook: "pan",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyRpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNS4xIE1hY2ludG9zaCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDpCRTI5MDhEODIwQjUxMUU0ODREQUYzNzM5QTM2MjBCRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDpCRTI5MDhEOTIwQjUxMUU0ODREQUYzNzM5QTM2MjBCRSI+IDx4bXBNTTpEZXJpdmVkRnJvbSBzdFJlZjppbnN0YW5jZUlEPSJ4bXAuaWlkOkJFMjkwOEQ2MjBCNTExRTQ4NERBRjM3MzlBMzYyMEJFIiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOkJFMjkwOEQ3MjBCNTExRTQ4NERBRjM3MzlBMzYyMEJFIi8+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+OXzPwwAAAKNJREFUeNrsVsEKgCAM3cyj0f8fuwT9XdEHrLyVIOKYY4kPPDim0+fenF+3HZi4nhFec+Rs4oCPAALwjDVUsKMWA6DNAFX6YXcMYIERdRWIYBzAZbKYGsSKex6mVUAK8Za0TphgoFTbpSvlx3/I0EQOILO2i/ibegLk/mgVONM4JvuBVizgkGH3XTGrR/xlV0ycbO8qCeMN54wdtVQwSTFwCzAATqEZUn8W8W4AAAAASUVORK5CYII=",
        cursor: "move",
        auto_deactivate: true,
        restrict_to_innercanvas: true
      };

      PanToolView.prototype.tool_events = {
        UpdatingMouseMove: "_drag",
        SetBasepoint: "_set_base_point"
      };

      PanToolView.prototype.mouse_coords = function(e, x, y) {
        var x_, y_, _ref1;
        _ref1 = [this.plot_view.canvas.sx_to_vx(x), this.plot_view.canvas.sy_to_vy(y)], x_ = _ref1[0], y_ = _ref1[1];
        return [x_, y_];
      };

      PanToolView.prototype._set_base_point = function(e) {
        var _ref1;
        _ref1 = this.mouse_coords(e, e.bokehX, e.bokehY), this.x = _ref1[0], this.y = _ref1[1];
        return null;
      };

      PanToolView.prototype._drag = function(e) {
        var dims, end, mapper, name, pan_info, sdx, sdy, start, sx0, sx1, sx_high, sx_low, sy0, sy1, sy_high, sy_low, x, xdiff, xr, xrs, y, ydiff, yr, yrs, _ref1, _ref2, _ref3, _ref4, _ref5, _ref6;
        _ref1 = this.mouse_coords(e, e.bokehX, e.bokehY), x = _ref1[0], y = _ref1[1];
        xdiff = x - this.x;
        ydiff = y - this.y;
        _ref2 = [x, y], this.x = _ref2[0], this.y = _ref2[1];
        xr = this.plot_view.frame.get('h_range');
        sx_low = xr.get('start') - xdiff;
        sx_high = xr.get('end') - xdiff;
        yr = this.plot_view.frame.get('v_range');
        sy_low = yr.get('start') - ydiff;
        sy_high = yr.get('end') - ydiff;
        dims = this.mget('dimensions');
        if (dims.indexOf('width') > -1) {
          sx0 = sx_low;
          sx1 = sx_high;
          sdx = -xdiff;
        } else {
          sx0 = xr.get('start');
          sx1 = xr.get('end');
          sdx = 0;
        }
        if (dims.indexOf('height') > -1) {
          sy0 = sy_low;
          sy1 = sy_high;
          sdy = ydiff;
        } else {
          sy0 = yr.get('start');
          sy1 = yr.get('end');
          sdy = 0;
        }
        xrs = {};
        _ref3 = this.plot_view.frame.get('x_mappers');
        for (name in _ref3) {
          mapper = _ref3[name];
          _ref4 = mapper.v_map_from_target([sx0, sx1]), start = _ref4[0], end = _ref4[1];
          xrs[name] = {
            start: start,
            end: end
          };
        }
        yrs = {};
        _ref5 = this.plot_view.frame.get('y_mappers');
        for (name in _ref5) {
          mapper = _ref5[name];
          _ref6 = mapper.v_map_from_target([sy0, sy1]), start = _ref6[0], end = _ref6[1];
          yrs[name] = {
            start: start,
            end: end
          };
        }
        pan_info = {
          xrs: xrs,
          yrs: yrs,
          sdx: sdx,
          sdy: sdy
        };
        this.plot_view.update_range(pan_info);
        return null;
      };

      return PanToolView;

    })(Tool.View);
    PanTool = (function(_super) {
      __extends(PanTool, _super);

      function PanTool() {
        _ref1 = PanTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      PanTool.prototype.default_view = PanToolView;

      PanTool.prototype.type = "PanTool";

      PanTool.prototype.defaults = function() {
        return {
          dimensions: ["width", "height"]
        };
      };

      PanTool.prototype.display_defaults = function() {
        return PanTool.__super__.display_defaults.call(this);
      };

      return PanTool;

    })(Tool.Model);
    PanTools = (function(_super) {
      __extends(PanTools, _super);

      function PanTools() {
        _ref2 = PanTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      PanTools.prototype.model = PanTool;

      return PanTools;

    })(Backbone.Collection);
    return {
      "Model": PanTool,
      "Collection": new PanTools(),
      "View": PanToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=pan_tool.js.map
*/