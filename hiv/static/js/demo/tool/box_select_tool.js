(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "./tool", "./event_generators"], function(_, Backbone, Tool, EventGenerators) {
    var BoxSelectTool, BoxSelectToolView, BoxSelectTools, TwoPointEventGenerator, _ref, _ref1, _ref2;
    TwoPointEventGenerator = EventGenerators.TwoPointEventGenerator;
    BoxSelectToolView = (function(_super) {
      __extends(BoxSelectToolView, _super);

      function BoxSelectToolView() {
        _ref = BoxSelectToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      BoxSelectToolView.prototype.initialize = function(options) {
        BoxSelectToolView.__super__.initialize.call(this, options);
        return this.select_every_mousemove = this.mget('select_every_mousemove');
      };

      BoxSelectToolView.prototype.bind_bokeh_events = function() {
        var renderer, rendererview, _i, _len, _ref1, _results;
        BoxSelectToolView.__super__.bind_bokeh_events.call(this);
        _ref1 = this.mget('renderers');
        _results = [];
        for (_i = 0, _len = _ref1.length; _i < _len; _i++) {
          renderer = _ref1[_i];
          rendererview = this.plot_view.renderers[renderer.id];
          this.listenTo(rendererview.xrange(), 'change', this.select_callback);
          this.listenTo(rendererview.yrange(), 'change', this.select_callback);
          _results.push(this.listenTo(renderer, 'change', this.select_callback));
        }
        return _results;
      };

      BoxSelectToolView.prototype.eventGeneratorClass = TwoPointEventGenerator;

      BoxSelectToolView.prototype.toolType = "BoxSelectTool";

      BoxSelectToolView.prototype.evgen_options = {
        keyName: "shiftKey",
        buttonText: "Select",
        buttonHook: "select",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACQAAAAgCAYAAAB6kdqOAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyRpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNS4xIE1hY2ludG9zaCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDpBODVDNDBCRjIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDpBODVDNDBDMDIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRSI+IDx4bXBNTTpEZXJpdmVkRnJvbSBzdFJlZjppbnN0YW5jZUlEPSJ4bXAuaWlkOkE4NUM0MEJEMjBCMzExRTQ4NERBRjM3MzlBMzYyMEJFIiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOkE4NUM0MEJFMjBCMzExRTQ4NERBRjM3MzlBMzYyMEJFIi8+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+hdQ7dQAAAJdJREFUeNpiXLhs5X8GBPgIxAJQNjZxfiD+wIAKGCkUZ0SWZGIYZIAF3YVoPkEHH6kojhUMyhD6jydEaAlgaWnwh9BAgf9DKpfxDxYHjeay0Vw2bHMZw2guG81lwyXKRnMZWlt98JdDTFAX/x9NQwPkIH6kGMAVEyjyo7lstC4jouc69Moh9L42rlyBTZyYXDS00xBAgAEAqsguPe03+cYAAAAASUVORK5CYII=",
        cursor: "crosshair",
        restrict_to_innercanvas: true
      };

      BoxSelectToolView.prototype.tool_events = {
        SetBasepoint: "_start_selecting",
        UpdatingMouseMove: "_selecting",
        deactivated: "_stop_selecting",
        DragEnd: "_dragend"
      };

      BoxSelectToolView.prototype.pause = function() {
        return null;
      };

      BoxSelectToolView.prototype.view_coords = function(sx, sy) {
        var vx, vy, _ref1;
        _ref1 = [this.plot_view.canvas.sx_to_vx(sx), this.plot_view.canvas.sy_to_vy(sy)], vx = _ref1[0], vy = _ref1[1];
        return [vx, vy];
      };

      BoxSelectToolView.prototype._stop_selecting = function() {
        this.trigger('stopselect');
        this.basepoint_set = false;
        return this.plot_view.unpause();
      };

      BoxSelectToolView.prototype._start_selecting = function(e) {
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

      BoxSelectToolView.prototype._get_selection_range = function() {
        var range, xrange, yrange;
        if (this.mget('select_x')) {
          xrange = [this.mget('start_vx'), this.mget('current_vx')];
          xrange = [_.min(xrange), _.max(xrange)];
        } else {
          range = this.plot_view.frame.get('h_range');
          xrange = [range.get('start'), range.get('end')];
        }
        if (this.mget('select_y')) {
          yrange = [this.mget('start_vy'), this.mget('current_vy')];
          yrange = [_.min(yrange), _.max(yrange)];
        } else {
          range = this.plot_view.frame.get('v_range');
          yrange = [range.get('start'), range.get('end')];
        }
        return [xrange, yrange];
      };

      BoxSelectToolView.prototype._selecting = function(e, x_, y_) {
        var vx, vy, _ref1, _ref2;
        _ref1 = this.view_coords(e.bokehX, e.bokehY), vx = _ref1[0], vy = _ref1[1];
        this.mset({
          'current_vx': vx,
          'current_vy': vy
        });
        _ref2 = this._get_selection_range(), this.xrange = _ref2[0], this.yrange = _ref2[1];
        this.trigger('boxselect', this.xrange, this.yrange);
        if (this.select_every_mousemove) {
          this._select_data();
        }
        this.plot_view._render_levels(this.plot_view.canvas_view.ctx, ['overlay']);
        return null;
      };

      BoxSelectToolView.prototype._dragend = function() {
        return this._select_data();
      };

      BoxSelectToolView.prototype._select_data = function() {
        var datasource, datasource_id, datasource_selections, datasources, ds, geometry, k, renderer, selected, v, _i, _j, _len, _len1, _ref1, _ref2;
        if (!this.basepoint_set) {
          return;
        }
        geometry = {
          type: 'rect',
          vx0: this.xrange[0],
          vx1: this.xrange[1],
          vy0: this.yrange[0],
          vy1: this.yrange[1]
        };
        datasources = {};
        datasource_selections = {};
        _ref1 = this.mget('renderers');
        for (_i = 0, _len = _ref1.length; _i < _len; _i++) {
          renderer = _ref1[_i];
          datasource = renderer.get('data_source');
          datasources[datasource.id] = datasource;
        }
        _ref2 = this.mget('renderers');
        for (_j = 0, _len1 = _ref2.length; _j < _len1; _j++) {
          renderer = _ref2[_j];
          datasource_id = renderer.get('data_source').id;
          _.setdefault(datasource_selections, datasource_id, []);
          selected = this.plot_view.renderers[renderer.id].hit_test(geometry);
          datasource_selections[datasource_id].push(selected);
        }
        for (k in datasource_selections) {
          if (!__hasProp.call(datasource_selections, k)) continue;
          v = datasource_selections[k];
          selected = _.intersection.apply(_, v);
          ds = datasources[k];
          ds.save({
            selected: selected
          }, {
            patch: true
          });
          this.plot_view.unpause();
        }
        return null;
      };

      return BoxSelectToolView;

    })(Tool.View);
    BoxSelectTool = (function(_super) {
      __extends(BoxSelectTool, _super);

      function BoxSelectTool() {
        _ref1 = BoxSelectTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      BoxSelectTool.prototype.default_view = BoxSelectToolView;

      BoxSelectTool.prototype.type = "BoxSelectTool";

      BoxSelectTool.prototype.defaults = function() {
        return _.extend(BoxSelectTool.__super__.defaults.call(this), {
          renderers: [],
          select_x: true,
          select_y: true,
          select_every_mousemove: false,
          data_source_options: {}
        });
      };

      BoxSelectTool.prototype.display_defaults = function() {
        return BoxSelectTool.__super__.display_defaults.call(this);
      };

      return BoxSelectTool;

    })(Tool.Model);
    BoxSelectTools = (function(_super) {
      __extends(BoxSelectTools, _super);

      function BoxSelectTools() {
        _ref2 = BoxSelectTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      BoxSelectTools.prototype.model = BoxSelectTool;

      return BoxSelectTools;

    })(Backbone.Collection);
    return {
      "Model": BoxSelectTool,
      "Collection": new BoxSelectTools(),
      "View": BoxSelectToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=box_select_tool.js.map
*/