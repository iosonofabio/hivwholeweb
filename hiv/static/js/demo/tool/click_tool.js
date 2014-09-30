(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "./tool"], function(_, Backbone, Tool) {
    var ClickTool, ClickToolView, ClickTools, _ref, _ref1, _ref2;
    ClickToolView = (function(_super) {
      __extends(ClickToolView, _super);

      function ClickToolView() {
        _ref = ClickToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      ClickToolView.prototype.initialize = function(options) {
        ClickToolView.__super__.initialize.call(this, options);
        this.listenTo(this, 'clicked', (function(selected, ds) {
          return console.log(selected, ds);
        }));
        return this.active = false;
      };

      ClickToolView.prototype.view_coords = function(sx, sy) {
        var vx, vy, _ref1;
        _ref1 = [this.plot_view.canvas.sx_to_vx(sx), this.plot_view.canvas.sy_to_vy(sy)], vx = _ref1[0], vy = _ref1[1];
        return [vx, vy];
      };

      ClickToolView.prototype.bind_bokeh_events = function() {
        var tool_name,
          _this = this;
        tool_name = "click_tool";
        if (!this.mget('always_active')) {
          this.tool_button = $("<button class='bk-toolbar-button'> Click </button>");
          this.plot_view.$el.find('.bk-button-bar').append(this.tool_button);
          this.tool_button.click(function() {
            if (_this.active) {
              return _this.plot_view.eventSink.trigger("clear_active_tool");
            } else {
              return _this.plot_view.eventSink.trigger("active_tool", tool_name);
            }
          });
          this.plot_view.eventSink.on("" + tool_name + ":deactivated", function() {
            _this.active = false;
            return _this.tool_button.removeClass('active');
          });
          this.plot_view.eventSink.on("" + tool_name + ":activated", function() {
            _this.active = true;
            return _this.tool_button.addClass('active');
          });
        }
        this.plot_view.canvas_view.canvas_wrapper.bind("mousedown", function(e) {
          if (!_this.active && !_this.mget('always_active')) {
            return;
          }
          _this.start_posx = e.pageX;
          return _this.start_posy = e.pageY;
        });
        return this.plot_view.canvas_view.canvas_wrapper.bind("mouseup", function(e) {
          var left, offset, top, vx, vy, _ref1;
          if (!_this.active && !_this.mget('always_active')) {
            return;
          }
          if (_this.start_posx !== e.pageX || _this.start_posy !== e.pageY) {
            return;
          }
          offset = $(e.currentTarget).offset();
          left = offset != null ? offset.left : 0;
          top = offset != null ? offset.top : 0;
          e.bokehX = e.pageX - left;
          e.bokehY = e.pageY - top;
          _ref1 = _this.view_coords(e.bokehX, e.bokehY), vx = _ref1[0], vy = _ref1[1];
          return _this._select(vx, vy, e);
        });
      };

      ClickToolView.prototype._select = function(vx, vy, e) {
        var datasource, datasource_id, datasource_selections, datasources, ds, geometry, renderer, renderers, selected, x, xmapper, y, ymapper, _i, _j, _len, _len1;
        geometry = {
          type: 'point',
          vx: vx,
          vy: vy
        };
        datasources = {};
        datasource_selections = {};
        renderers = this.mget('renderers');
        for (_i = 0, _len = renderers.length; _i < _len; _i++) {
          renderer = renderers[_i];
          datasource = renderer.get('data_source');
          datasources[datasource.id] = datasource;
        }
        for (_j = 0, _len1 = renderers.length; _j < _len1; _j++) {
          renderer = renderers[_j];
          datasource_id = renderer.get('data_source').id;
          _.setdefault(datasource_selections, datasource_id, []);
          selected = this.plot_view.renderers[renderer.id].hit_test(geometry);
          ds = datasources[datasource_id];
          xmapper = this.plot_view.frame.get('x_mappers')[renderer.get('x_range_name')];
          ymapper = this.plot_view.frame.get('y_mappers')[renderer.get('y_range_name')];
          x = xmapper.map_from_target(vx);
          y = ymapper.map_from_target(vy);
          if (selected === null) {
            continue;
          }
          if (selected.length === 0) {
            continue;
          }
          this.trigger('clicked', selected, ds);
        }
        return null;
      };

      return ClickToolView;

    })(Tool.View);
    ClickTool = (function(_super) {
      __extends(ClickTool, _super);

      function ClickTool() {
        _ref1 = ClickTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      ClickTool.prototype.default_view = ClickToolView;

      ClickTool.prototype.type = "ClickTool";

      ClickTool.prototype.initialize = function(attrs, options) {
        var all_renderers, names, r, renderers;
        ClickTool.__super__.initialize.call(this, attrs, options);
        names = this.get('names');
        all_renderers = this.get('plot').get('renderers');
        renderers = (function() {
          var _i, _len, _results;
          _results = [];
          for (_i = 0, _len = all_renderers.length; _i < _len; _i++) {
            r = all_renderers[_i];
            if (r.type === "Glyph") {
              _results.push(r);
            }
          }
          return _results;
        })();
        if (names.length > 0) {
          renderers = (function() {
            var _i, _len, _results;
            _results = [];
            for (_i = 0, _len = renderers.length; _i < _len; _i++) {
              r = renderers[_i];
              if (names.indexOf(r.get('name')) >= 0) {
                _results.push(r);
              }
            }
            return _results;
          })();
        }
        return this.set('renderers', renderers);
      };

      ClickTool.prototype.defaults = function() {
        return _.extend(ClickTool.__super__.defaults.call(this), {
          renderers: [],
          names: [],
          always_active: []
        });
      };

      ClickTool.prototype.display_defaults = function() {
        return ClickTool.__super__.display_defaults.call(this);
      };

      return ClickTool;

    })(Tool.Model);
    ClickTools = (function(_super) {
      __extends(ClickTools, _super);

      function ClickTools() {
        _ref2 = ClickTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      ClickTools.prototype.model = ClickTool;

      return ClickTools;

    })(Backbone.Collection);
    return {
      "Model": ClickTool,
      "Collection": new ClickTools(),
      "View": ClickToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=click_tool.js.map
*/