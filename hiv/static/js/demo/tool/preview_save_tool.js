(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "jquery", "bootstrap/modal", "backbone", "common/bulk_save", "./tool", "./event_generators", "./preview_save_tool_template"], function(_, $, $$1, Backbone, bulk_save, Tool, EventGenerators, preview_save_tool_template) {
    var ButtonEventGenerator, PreviewSaveTool, PreviewSaveToolView, PreviewSaveTools, _ref, _ref1, _ref2;
    ButtonEventGenerator = EventGenerators.ButtonEventGenerator;
    PreviewSaveToolView = (function(_super) {
      __extends(PreviewSaveToolView, _super);

      function PreviewSaveToolView() {
        _ref = PreviewSaveToolView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      PreviewSaveToolView.prototype.initialize = function(options) {
        return PreviewSaveToolView.__super__.initialize.call(this, options);
      };

      PreviewSaveToolView.prototype.eventGeneratorClass = ButtonEventGenerator;

      PreviewSaveToolView.prototype.evgen_options = {
        buttonText: "Preview/Save",
        buttonHook: "preview-save",
        buttonIcon: "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyRpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENTNS4xIE1hY2ludG9zaCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDozMjFERDhENjIwQjIxMUU0ODREQUYzNzM5QTM2MjBCRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDozMjFERDhENzIwQjIxMUU0ODREQUYzNzM5QTM2MjBCRSI+IDx4bXBNTTpEZXJpdmVkRnJvbSBzdFJlZjppbnN0YW5jZUlEPSJ4bXAuaWlkOjMyMUREOEQ0MjBCMjExRTQ4NERBRjM3MzlBMzYyMEJFIiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOjMyMUREOEQ1MjBCMjExRTQ4NERBRjM3MzlBMzYyMEJFIi8+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+h5hT8AAAAKBJREFUeNpiWbhs5QcGBgZ+hgECTAwDDGAO+AjEjGj4Lw5xUrAAkl3ocr8IhQAzjT3PRu0o+I+EHw65NDDqgJHrABYC8t9JMIuRmiHACS2IKC0LOKH0X1JDAOTzs0BsBs3XlIKz5KSBRCA+RQXLjwNxNDlp4BoQm9Mo7fGPZsNRB4w6YNQBI94BfwfaAV9G08CoA9DbA/xUavkMvRAACDAAaPgYViexODkAAAAASUVORK5CYII="
      };

      PreviewSaveToolView.prototype.toolType = "PreviewSaveTool";

      PreviewSaveToolView.prototype.tool_events = {
        activated: "_activated",
        deactivated: "_close_modal"
      };

      PreviewSaveToolView.prototype._activated = function(e) {
        var data_uri,
          _this = this;
        data_uri = this.plot_view.canvas_view.canvas[0].toDataURL();
        this.plot_model.set('png', this.plot_view.canvas_view.canvas[0].toDataURL());
        this.$modal = $(preview_save_tool_template({
          data_uri: data_uri
        }));
        $('body').append(this.$modal);
        this.$modal.on('hidden', function() {
          return _this.plot_view.eventSink.trigger("clear_active_tool");
        });
        return this.$modal.modal({
          show: true
        });
      };

      PreviewSaveToolView.prototype._close_modal = function() {
        return this.$modal.remove();
      };

      return PreviewSaveToolView;

    })(Tool.View);
    PreviewSaveTool = (function(_super) {
      __extends(PreviewSaveTool, _super);

      function PreviewSaveTool() {
        _ref1 = PreviewSaveTool.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      PreviewSaveTool.prototype.default_view = PreviewSaveToolView;

      PreviewSaveTool.prototype.type = "PreviewSaveTool";

      PreviewSaveTool.prototype.display_defaults = function() {
        return PreviewSaveTool.__super__.display_defaults.call(this);
      };

      return PreviewSaveTool;

    })(Tool.Model);
    PreviewSaveTools = (function(_super) {
      __extends(PreviewSaveTools, _super);

      function PreviewSaveTools() {
        _ref2 = PreviewSaveTools.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      PreviewSaveTools.prototype.model = PreviewSaveTool;

      return PreviewSaveTools;

    })(Backbone.Collection);
    return {
      Model: PreviewSaveTool,
      Collection: new PreviewSaveTools(),
      View: PreviewSaveToolView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=preview_save_tool.js.map
*/