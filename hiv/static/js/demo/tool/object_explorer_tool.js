(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "jquery", "bootstrap/modal", "backbone", "./tool", "./event_generators", "./object_explorer_tool_template", "widget/object_explorer"], function(_, $, $$1, Backbone, Tool, EventGenerators, object_explorer_tool_template, ObjectExplorer) {
    var ButtonEventGenerator, ObjectExplorerTool, ObjectExplorerToolView, ObjectExplorerTools, _ref, _ref1, _ref2;
    ButtonEventGenerator = EventGenerators.ButtonEventGenerator;
    ObjectExplorerToolView = (function(_super) {
      __extends(ObjectExplorerToolView, _super);

      function ObjectExplorerToolView() {
        _ref = ObjectExplorerToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      ObjectExplorerToolView.prototype.initialize = function(options) {
        return ObjectExplorerToolView.__super__.initialize.call(this, options);
      };

      ObjectExplorerToolView.prototype.eventGeneratorClass = ButtonEventGenerator;

      ObjectExplorerToolView.prototype.evgen_options = {
        buttonText: "Object Explorer",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABwAAAAcCAYAAAByDd+UAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAABx0RVh0U29mdHdhcmUAQWRvYmUgRmlyZXdvcmtzIENTNui8sowAAAJHSURBVEiJvZbPTlNREManQmtSvTvqI7gw4R1cQpP6ABoahKWURxBfgo0ulRfgGYA9lI1dIdEFMSo0sZLY3p+L+514eu85t7dKOskk58838+XMmTNnaoAtUu4tlM3MDLCqpwQ2KMrGPDzREwLLQJKmaU3zB2a2pu1fUjOzNaBpZpamaQ1IgOWZzLm1OrADfAQ6wCOgDfwAJsCudKK1dWGeyeYVUA/y5AmBhsh8+QqMNb4EHksvtTYRxpcdoFFKCCx5ZCOgD9xofgOcAm8UgbrGpzlMX7aOdKmMMAEGAve1tgUcA92Sa+kKs6V5Xz4GQBIl1KV3gN/AN2Azevlx8k3ZjoGOl3ThO9SmC9HJPxCeuPDm1ovPQiFtm9l9MxuZ2X7AYYssY9tAK8C5L9uGMA/zDpx2gQ/e6c4CZKvAoZeJh8BqAHfmJdF7soIxHVKm5ROwl3OyAhxp/4sUra3ksHvy4UuQcET2oJ9QfLhtYa6Ap9Irra3nsHX52OXvE4mWNqRVpAwX9hMJ6QXwOodreSH9LC0L6cWskLqkGQpwF0kzlM9ugVDARHd1C/wEngeclT4L4IVsb4WJVxrPyD2N4/D1xIV5Hr4rbWbWNLPvZvZuXkIzeyvbJl5pKzBrHCveR8B2yam2hXmpebXiTfn3NATOKX5P516iXTPP96SNBtBjWlz1h6yCuA/YVZOxML70mPUBe5uuxRiQtQ2uxbgm+9170lCLMZBttRbDAxSaKODAC7cL2wEVmqiZhJHk+O82sVaV7K5k4Z33H/QTdNyD5wyAAAAAAElFTkSuQmCC"
      };

      ObjectExplorerToolView.prototype.toolType = "ObjectExplorerTool";

      ObjectExplorerToolView.prototype.tool_events = {
        activated: "_activated",
        deactivated: "_close_modal"
      };

      ObjectExplorerToolView.prototype._activated = function(e) {
        var _this = this;
        this.$modal = $(object_explorer_tool_template({}));
        this.$object_explorer_view = new ObjectExplorer.View({
          el: this.$modal.find(".bk-bs-modal-body")
        });
        $('body').append(this.$modal);
        this.$modal.on('hidden', function() {
          return _this.plot_view.eventSink.trigger("clear_active_tool");
        });
        return this.$modal.modal({
          show: true
        });
      };

      ObjectExplorerToolView.prototype._close_modal = function() {
        return this.$modal.remove();
      };

      return ObjectExplorerToolView;

    })(Tool.View);
    ObjectExplorerTool = (function(_super) {
      __extends(ObjectExplorerTool, _super);

      function ObjectExplorerTool() {
        _ref1 = ObjectExplorerTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      ObjectExplorerTool.prototype.default_view = ObjectExplorerToolView;

      ObjectExplorerTool.prototype.type = "ObjectExplorerTool";

      ObjectExplorerTool.prototype.display_defaults = function() {
        return ObjectExplorerTool.__super__.display_defaults.call(this);
      };

      return ObjectExplorerTool;

    })(Tool.Model);
    ObjectExplorerTools = (function(_super) {
      __extends(ObjectExplorerTools, _super);

      function ObjectExplorerTools() {
        _ref2 = ObjectExplorerTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      ObjectExplorerTools.prototype.model = ObjectExplorerTool;

      return ObjectExplorerTools;

    })(Backbone.Collection);
    return {
      Model: ObjectExplorerTool,
      Collection: new ObjectExplorerTools(),
      View: ObjectExplorerToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=object_explorer_tool.js.map
*/