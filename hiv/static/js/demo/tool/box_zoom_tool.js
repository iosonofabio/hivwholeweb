(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "./tool", "./event_generators"], function(_, Backbone, Tool, EventGenerators) {
    var BoxZoomTool, BoxZoomToolView, BoxZoomTools, TwoPointEventGenerator, _ref, _ref1, _ref2;
    TwoPointEventGenerator = EventGenerators.TwoPointEventGenerator;
    BoxZoomToolView = (function(_super) {
      __extends(BoxZoomToolView, _super);

      function BoxZoomToolView() {
        _ref = BoxZoomToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      BoxZoomToolView.prototype.initialize = function(options) {
        return BoxZoomToolView.__super__.initialize.call(this, options);
      };

      BoxZoomToolView.prototype.bind_bokeh_events = function() {
        return BoxZoomToolView.__super__.bind_bokeh_events.call(this);
      };

      BoxZoomToolView.prototype.eventGeneratorClass = TwoPointEventGenerator;

      BoxZoomToolView.prototype.toolType = "BoxZoomTool";

      BoxZoomToolView.prototype.evgen_options = {
        keyName: "ctrlKey",
        buttonText: "Box Zoom",
        buttonHook: "box-zoom",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACIAAAAgCAYAAAB3j6rJAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyRpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNS4xIE1hY2ludG9zaCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDozMjFERDhEMjIwQjIxMUU0ODREQUYzNzM5QTM2MjBCRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDozMjFERDhEMzIwQjIxMUU0ODREQUYzNzM5QTM2MjBCRSI+IDx4bXBNTTpEZXJpdmVkRnJvbSBzdFJlZjppbnN0YW5jZUlEPSJ4bXAuaWlkOjMyMUREOEQwMjBCMjExRTQ4NERBRjM3MzlBMzYyMEJFIiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOjMyMUREOEQxMjBCMjExRTQ4NERBRjM3MzlBMzYyMEJFIi8+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+a2Q0KAAAAmVJREFUeNq8V19EpFEUvzOtmKfpJSJKDL2WiLJExKaUEq0eeikiaolZLT2lVUpPydqHqIlIo1ilFOmphxj1miKWWHppnobIt7+zeyZ3jjvz/bnf9OPHd8/9d77z3XN+94ts7ew6SqksWKX+w1GFiLjYdVSAfeAQ2Ag2sf0GvAXT4C/wle1x3lt9UOGBNk6BrYa+FuYIeAWOsmNviGqe6W+q081OmAGvizgh0cpjZ3RjGBFZBpMG+xn4wM8NYJfWFwNXwXrwS96RiIUTwwYn6AxMgb+FvQ5c4zOUxzR4Ce5GLZyo5LfSsQP2G5xQbKO+bWFfoLWinA1OAEcoM2rFRpMe5sloJWgtm4j0iPZcPhVdkOWxBWvZONIi2uc+5sqxbTaO1Ij2o4+5T6JdGy1SF4Kg2mLsi01E/oh2l4+5HTKaNlmTEe0ka40XyNqTsYnIkWiTwC16rMRNci0bR0hJ7w1veizqy9uB5D4ZDZKBtI3WvLCCJoT9E3jHny4j1DdmWOcbrWWjNYuGoqaL2kdmKayTztio7yzTJprz4A/9PuI3a8YMh5IKVC9fetxAY5rB79pNzXdESMJ/GrSjm8/DCTjAgpjQZCDDh5I+w4HuQBBHOsE9USty4KB2KF85m9J+v5XX9KXr3T7fQZS26WefYlcU+ayJlxhDIT40jBnn21hQOPrfgFtEqAhdGETqK7gZ4h/Av4g4Jf5TUoYquQSuqJDhFpEJca3b4EoYOtyyhrSkHTzlcj4R4t4FZ9NL+j6yMzlT/ocZES9aky3D3r6y5t2gaw3xWXgs7XFhdyzsgSpr2fFXgAEAmp2J9DuX/WgAAAAASUVORK5CYII=",
        cursor: "crosshair",
        auto_deactivate: true,
        restrict_to_innercanvas: true
      };

      BoxZoomToolView.prototype.tool_events = {
        SetBasepoint: "_start_selecting",
        UpdatingMouseMove: "_selecting",
        DragEnd: "_dragend"
      };

      BoxZoomToolView.prototype.pause = function() {
        return null;
      };

      BoxZoomToolView.prototype.view_coords = function(sx, sy) {
        var vx, vy, _ref1;
        _ref1 = [this.plot_view.canvas.sx_to_vx(sx), this.plot_view.canvas.sy_to_vy(sy)], vx = _ref1[0], vy = _ref1[1];
        return [vx, vy];
      };

      BoxZoomToolView.prototype._start_selecting = function(e) {
        var vx, vy, _ref1;
        this.plot_view.pause();
        this.trigger('startselect');
        _ref1 = this.view_coords(e.bokehX, e.bokehY), vx = _ref1[0], vy = _ref1[1];
        this.mset({
          'start_vx': vx,
          'start_vy': vy,
          'current_vx': null,
          'current_vy': null
        });
        return this.basepoint_set = true;
      };

      BoxZoomToolView.prototype._get_selection_range = function() {
        var xrange, yrange;
        if (this.mget('select_x')) {
          xrange = [this.mget('start_vx'), this.mget('current_vx')];
          xrange = [_.min(xrange), _.max(xrange)];
        } else {
          xrange = null;
        }
        if (this.mget('select_y')) {
          yrange = [this.mget('start_vy'), this.mget('current_vy')];
          yrange = [_.min(yrange), _.max(yrange)];
        } else {
          yrange = null;
        }
        return [xrange, yrange];
      };

      BoxZoomToolView.prototype._selecting = function(e, x_, y_) {
        var vx, vy, _ref1, _ref2;
        _ref1 = this.view_coords(e.bokehX, e.bokehY), vx = _ref1[0], vy = _ref1[1];
        this.mset({
          'current_vx': vx,
          'current_vy': vy
        });
        _ref2 = this._get_selection_range(), this.xrange = _ref2[0], this.yrange = _ref2[1];
        this.trigger('boxselect', this.xrange, this.yrange);
        this.plot_view._render_levels(this.plot_view.canvas_view.ctx, ['overlay']);
        return null;
      };

      BoxZoomToolView.prototype._dragend = function() {
        this._select_data();
        this.basepoint_set = false;
        this.plot_view.unpause();
        return this.trigger('stopselect');
      };

      BoxZoomToolView.prototype._select_data = function() {
        var end, mapper, name, start, xrs, yrs, zoom_info, _ref1, _ref2, _ref3, _ref4;
        if (!this.basepoint_set) {
          return;
        }
        xrs = {};
        _ref1 = this.plot_view.frame.get('x_mappers');
        for (name in _ref1) {
          mapper = _ref1[name];
          _ref2 = mapper.v_map_from_target([this.xrange[0], this.xrange[1]]), start = _ref2[0], end = _ref2[1];
          xrs[name] = {
            start: start,
            end: end
          };
        }
        yrs = {};
        _ref3 = this.plot_view.frame.get('y_mappers');
        for (name in _ref3) {
          mapper = _ref3[name];
          _ref4 = mapper.v_map_from_target([this.yrange[0], this.yrange[1]]), start = _ref4[0], end = _ref4[1];
          yrs[name] = {
            start: start,
            end: end
          };
        }
        zoom_info = {
          xrs: xrs,
          yrs: yrs
        };
        return this.plot_view.update_range(zoom_info);
      };

      return BoxZoomToolView;

    })(Tool.View);
    BoxZoomTool = (function(_super) {
      __extends(BoxZoomTool, _super);

      function BoxZoomTool() {
        _ref1 = BoxZoomTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      BoxZoomTool.prototype.default_view = BoxZoomToolView;

      BoxZoomTool.prototype.type = "BoxZoomTool";

      BoxZoomTool.prototype.defaults = function() {
        return _.extend(BoxZoomTool.__super__.defaults.call(this), {
          renderers: [],
          select_x: true,
          select_y: true,
          select_every_mousemove: false,
          data_source_options: {}
        });
      };

      BoxZoomTool.prototype.display_defaults = function() {
        return BoxZoomTool.__super__.display_defaults.call(this);
      };

      return BoxZoomTool;

    })(Tool.Model);
    BoxZoomTools = (function(_super) {
      __extends(BoxZoomTools, _super);

      function BoxZoomTools() {
        _ref2 = BoxZoomTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      BoxZoomTools.prototype.model = BoxZoomTool;

      return BoxZoomTools;

    })(Backbone.Collection);
    return {
      "Model": BoxZoomTool,
      "Collection": new BoxZoomTools(),
      "View": BoxZoomToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=box_zoom_tool.js.map
*/