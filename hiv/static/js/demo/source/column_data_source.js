(function() {
  var __hasProp = {}.hasOwnProperty,
    __extends = function(child, parent) { for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; };

  define(["underscore", "backbone", "common/has_properties"], function(_, Backbone, HasProperties) {
    var ColumnDataSource, ColumnDataSources, _ref, _ref1;
    ColumnDataSource = (function(_super) {
      __extends(ColumnDataSource, _super);

      function ColumnDataSource() {
        _ref = ColumnDataSource.__super__.constructor.apply(this, arguments);
        return _ref;
      }

      ColumnDataSource.prototype.type = 'ColumnDataSource';

      ColumnDataSource.prototype.get_column = function(colname) {
        var _ref1;
        return (_ref1 = this.get('data')[colname]) != null ? _ref1 : null;
      };

      ColumnDataSource.prototype.get_length = function() {
        var data;
        data = this.get('data');
        if (_.keys(data).length === 0) {
          return 0;
        }
        return data[_.keys(data)[0]].length;
      };

      ColumnDataSource.prototype.columns = function() {
        return _.keys(this.get('data'));
      };

      ColumnDataSource.prototype.datapoints = function() {
        var data, field, fields, i, point, points, _i, _j, _len, _ref1;
        data = this.get('data');
        fields = _.keys(data);
        if (fields.length === 0) {
          return [];
        }
        points = [];
        for (i = _i = 0, _ref1 = data[fields[0]].length; 0 <= _ref1 ? _i < _ref1 : _i > _ref1; i = 0 <= _ref1 ? ++_i : --_i) {
          point = {};
          for (_j = 0, _len = fields.length; _j < _len; _j++) {
            field = fields[_j];
            point[field] = data[field][i];
          }
          points.push(point);
        }
        return points;
      };

      ColumnDataSource.prototype.defaults = function() {
        return {
          data: {}
        };
      };

      return ColumnDataSource;

    })(HasProperties);
    ColumnDataSources = (function(_super) {
      __extends(ColumnDataSources, _super);

      function ColumnDataSources() {
        _ref1 = ColumnDataSources.__super__.constructor.apply(this, arguments);
        return _ref1;
      }

      ColumnDataSources.prototype.model = ColumnDataSource;

      return ColumnDataSources;

    })(Backbone.Collection);
    return {
      "Model": ColumnDataSource,
      "Collection": new ColumnDataSources()
    };
  });

}).call(this);

/*
//@ sourceMappingURL=column_data_source.js.map
*/