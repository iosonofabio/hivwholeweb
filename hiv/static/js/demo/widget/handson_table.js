(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "jquery", "handsontable", "backbone", "common/has_properties", "common/continuum_view"], function(_, $, $$1, Backbone, HasProperties, ContinuumView) {
    var HandsonTable, HandsonTableView, HandsonTables, _ref, _ref1, _ref2;
    HandsonTableView = (function(_super) {
      __extends(HandsonTableView, _super);

      function HandsonTableView() {
        _ref = HandsonTableView.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      HandsonTableView.prototype.initialize = function(options) {
        var source,
          _this = this;
        HandsonTableView.__super__.initialize.call(this, options);
        this.render();
        this.listenTo(this.model, 'change', function() {
          return _this.renderFn();
        });
        source = this.mget("source");
        this.listenTo(source, 'change:data', function() {
          return _this.renderFn();
        });
        return this.listenTo(source, 'change:selected', function() {
          return _this.changeSelection();
        });
      };

      HandsonTableView.prototype.changeSelection = function() {
        var i, j, n, selected;
        this.ht.deselectCell();
        selected = this.mget("source").get("selected");
        i = _.min(selected);
        j = _.max(selected);
        n = this.ht.countCols();
        return this.ht.selectCell(i, 0, j, n - 1, true);
      };

      HandsonTableView.prototype.renderFn = function() {
        var column, columns, headers, source, _i, _len, _ref1,
          _this = this;
        source = this.mget("source");
        if (source != null) {
          headers = [];
          columns = [];
          _ref1 = this.mget("columns");
          for (_i = 0, _len = _ref1.length; _i < _len; _i++) {
            column = _ref1[_i];
            if (column != null) {
              headers.push(column.get("header"));
              columns.push({
                data: column.get("field"),
                type: column.get("type"),
                format: column.get("format"),
                source: column.get("source"),
                strict: column.get("strict"),
                checkedTemplate: column.get("checked"),
                uncheckedTemplate: column.get("unchecked")
              });
            }
          }
          this.$el.handsontable({
            data: source.datapoints(),
            colHeaders: headers,
            columns: columns,
            columnSorting: this.mget("sorting"),
            rowHeaders: true,
            width: this.mget("width"),
            height: this.mget("height"),
            afterChange: function(changes, source) {
              if (source === "edit") {
                return _this.editData(changes);
              }
            }
          });
        } else {
          this.$el.handsontable();
        }
        return this.ht = this.$el.handsontable("getInstance");
      };

      HandsonTableView.prototype.render = function() {
        var display, interval,
          _this = this;
        display = this.$el.css("display");
        return interval = setInterval(function() {
          if (_this.$el.css("display") !== display) {
            clearInterval(interval);
            return _this.renderFn();
          }
        }, 100);
      };

      HandsonTableView.prototype.editData = function(changes) {
        var array, change, column, data, i, index, new_val, old_val, source, _i, _j, _len, _ref1;
        source = this.mget("source");
        data = source.get("data");
        for (_i = 0, _len = changes.length; _i < _len; _i++) {
          change = changes[_i];
          index = change[0], column = change[1], old_val = change[2], new_val = change[3];
          array = _.clone(data[column]);
          if (index < array.length) {
            array[index] = new_val;
          } else {
            for (i = _j = 0, _ref1 = array.length - index; 0 <= _ref1 ? _j < _ref1 : _j > _ref1; i = 0 <= _ref1 ? ++_j : --_j) {
              array.push(NaN);
            }
            array.push(new_val);
          }
          data[column] = array;
        }
        return source.set(data);
      };

      return HandsonTableView;

    })(ContinuumView.View);
    HandsonTable = (function(_super) {
      __extends(HandsonTable, _super);

      function HandsonTable() {
        _ref1 = HandsonTable.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      HandsonTable.prototype.type = 'HandsonTable';

      HandsonTable.prototype.default_view = HandsonTableView;

      HandsonTable.prototype.defaults = function() {
        return {
          source: null,
          columns: [],
          width: null,
          height: null
        };
      };

      return HandsonTable;

    })(HasProperties);
    HandsonTables = (function(_super) {
      __extends(HandsonTables, _super);

      function HandsonTables() {
        _ref2 = HandsonTables.__super__.constructor.apply(this, arguments);
        return _ref2;
      }

      HandsonTables.prototype.model = HandsonTable;

      return HandsonTables;

    })(Backbone.Collection);
    return {
      Model: HandsonTable,
      Collection: new HandsonTables(),
      View: HandsonTableView
    };
  });

}).call(this);

/*
//@ sourceMappingURL=handson_table.js.map
*/