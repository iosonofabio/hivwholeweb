(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "./tool", "./event_generators"], function(_, Backbone, Tool, EventGenerators) {
    var ResizeTool, ResizeToolView, ResizeTools, TwoPointEventGenerator, _ref, _ref1, _ref2;
    TwoPointEventGenerator = EventGenerators.TwoPointEventGenerator;
    ResizeToolView = (function(_super) {
      __extends(ResizeToolView, _super);

      function ResizeToolView() {
        _ref = ResizeToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      ResizeToolView.prototype.initialize = function(options) {
        ResizeToolView.__super__.initialize.call(this, options);
        return this.active = false;
      };

      ResizeToolView.prototype.bind_events = function(plotview) {
        return ResizeToolView.__super__.bind_events.call(this, plotview);
      };

      ResizeToolView.prototype.eventGeneratorClass = TwoPointEventGenerator;

      ResizeToolView.prototype.toolType = "ResizeTool";

      ResizeToolView.prototype.evgen_options = {
        keyName: "",
        buttonText: "Resize",
        buttonHook: "resize",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACIAAAAgCAYAAAB3j6rJAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyRpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNS4xIE1hY2ludG9zaCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDpBODVDNDBCQjIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDpBODVDNDBCQzIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRSI+IDx4bXBNTTpEZXJpdmVkRnJvbSBzdFJlZjppbnN0YW5jZUlEPSJ4bXAuaWlkOjMyMUREOEQ4MjBCMjExRTQ4NERBRjM3MzlBMzYyMEJFIiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOkE4NUM0MEJBMjBCMzExRTQ4NERBRjM3MzlBMzYyMEJFIi8+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+nIbQ0AAAAIJJREFUeNpiXLhs5X8G7ICRgTYAq31MDIMEwBzyERoCyJhWAN2ej4MqRFiIjUMahczgSyMsNE4PxACBQZlrcAFsuYkcLECpQwZNiIw6ZNQhow4ZdcioQ0YdMuoQerRZkQE/vdqwgypqQD7+MIBuANn9f1CnEcbRXIMjd4zM0QCAAAMAbdAPQaze1JcAAAAASUVORK5CYII=",
        cursor: "move"
      };

      ResizeToolView.prototype.tool_events = {
        activated: "_activate",
        deactivated: "_deactivate",
        UpdatingMouseMove: "_drag",
        SetBasepoint: "_set_base_point"
      };

      ResizeToolView.prototype.render = function() {
        var ch, ctx, cw, line_width;
        if (!this.active) {
          return;
        }
        ctx = this.plot_view.canvas_view.ctx;
        cw = this.plot_view.canvas.get('width');
        ch = this.plot_view.canvas.get('height');
        line_width = 8;
        ctx.save();
        ctx.strokeStyle = 'transparent';
        ctx.globalAlpha = 0.7;
        ctx.lineWidth = line_width;
        ctx.setLineDash([]);
        ctx.beginPath();
        ctx.rect(line_width, line_width, cw - line_width * 2, ch - line_width * 2);
        ctx.moveTo(line_width, line_width);
        ctx.lineTo(cw - line_width, ch - line_width);
        ctx.moveTo(line_width, ch - line_width);
        ctx.lineTo(cw - line_width, line_width);
        ctx.stroke();
        return ctx.restore();
      };

      ResizeToolView.prototype.mouse_coords = function(e, x, y) {
        return [x, y];
      };

      ResizeToolView.prototype._activate = function(e) {
        var bbar, ch, cw, plotarea;
        if (this.active) {
          return;
        }
        this.active = true;
        this.popup = $('<div class="resize_bokeh_plot pull-right hide"/>');
        bbar = this.plot_view.$el.find('.bokeh_canvas_wrapper');
        plotarea = this.plot_view.$el.find('.plotarea');
        this.popup.appendTo(bbar);
        ch = this.plot_view.canvas.get('height');
        cw = this.plot_view.canvas.get('width');
        this.plot_view.request_render(true);
        return null;
      };

      ResizeToolView.prototype._deactivate = function(e) {
        this.active = false;
        this.popup.remove();
        this.request_render();
        this.plot_view.request_render();
        return null;
      };

      ResizeToolView.prototype._set_base_point = function(e) {
        var _ref1;
        _ref1 = this.mouse_coords(e, e.bokehX, e.bokehY), this.x = _ref1[0], this.y = _ref1[1];
        return null;
      };

      ResizeToolView.prototype._drag = function(e) {
        var ch, cw, x, xdiff, y, ydiff, _ref1, _ref2;
        this.plot_view.pause();
        _ref1 = this.mouse_coords(e, e.bokehX, e.bokehY), x = _ref1[0], y = _ref1[1];
        xdiff = x - this.x;
        ydiff = y - this.y;
        _ref2 = [x, y], this.x = _ref2[0], this.y = _ref2[1];
        ch = this.plot_view.canvas.get('height');
        cw = this.plot_view.canvas.get('width');
        this.plot_view.canvas._set_dims([cw + xdiff, ch + ydiff]);
        this.plot_view.request_render();
        this.plot_view.unpause(true);
        return null;
      };

      return ResizeToolView;

    })(Tool.View);
    ResizeTool = (function(_super) {
      __extends(ResizeTool, _super);

      function ResizeTool() {
        _ref1 = ResizeTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      ResizeTool.prototype.default_view = ResizeToolView;

      ResizeTool.prototype.type = "ResizeTool";

      ResizeTool.prototype.display_defaults = function() {
        return ResizeTool.__super__.display_defaults.call(this);
      };

      return ResizeTool;

    })(Tool.Model);
    ResizeTools = (function(_super) {
      __extends(ResizeTools, _super);

      function ResizeTools() {
        _ref2 = ResizeTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      ResizeTools.prototype.model = ResizeTool;

      return ResizeTools;

    })(Backbone.Collection);
    return {
      "Model": ResizeTool,
      "Collection": new ResizeTools(),
      "View": ResizeToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=resize_tool.js.map
*/