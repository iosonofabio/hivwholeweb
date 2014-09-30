(function() {
  require(['main', 'underscore'], function(Bokeh, _) {
    var N, data, hover, i, options, plot, r, scatter, t, text, val, x, y;
    r = new Bokeh.Random(123456789);
    N = 30;
    x = (function() {
      var _i, _len, _ref, _results;
      _ref = _.range(N);
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        i = _ref[_i];
        _results.push(r.randf() * 100);
      }
      return _results;
    })();
    y = (function() {
      var _i, _len, _ref, _results;
      _ref = _.range(N);
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        i = _ref[_i];
        _results.push(r.randf() * 100);
      }
      return _results;
    })();
    t = (function() {
      var _i, _ref, _results;
      _results = [];
      for (i = _i = 0, _ref = x.length; 0 <= _ref ? _i <= _ref : _i >= _ref; i = 0 <= _ref ? ++_i : --_i) {
        _results.push("" + i);
      }
      return _results;
    })();
    data = {
      x: x,
      y: y,
      radius: (function() {
        var _i, _len, _ref, _results;
        _ref = _.range(N);
        _results = [];
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          i = _ref[_i];
          _results.push(r.randf() * 2 + 2);
        }
        return _results;
      })(),
      color: (function() {
        var _i, _len, _ref, _results;
        _ref = _.zip(x, y);
        _results = [];
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          val = _ref[_i];
          _results.push("rgb(" + (Math.floor(50 + 2 * val[0])) + ", " + (Math.floor(30 + 2 * val[1])) + ", 150)");
        }
        return _results;
      })(),
      text: t
    };
    scatter = {
      type: 'circle',
      x: 'x',
      y: 'y',
      radius: 'radius',
      radius_units: 'data',
      fill_color: 'color',
      fill_alpha: 0.6,
      line_color: null,
      name: "mydata"
    };
    text = {
      type: 'text',
      x: 'x',
      y: 'y',
      angle: 0,
      text: {
        field: 'text'
      },
      text_font_size: '8pt',
      text_color: 'black',
      text_alpha: 0.4,
      text_align: 'center'
    };
    options = {
      title: "Hover Demo",
      dims: [600, 600],
      xrange: [0, 100],
      yrange: [0, 100],
      xaxes: "below",
      yaxes: "left",
      tools: "pan,wheel_zoom,select,resize,preview,reset,box_zoom,hover",
      legend: false
    };
    plot = Bokeh.Plotting.make_plot([scatter, text], data, options);
    hover = _.find(plot.get('tools'), function(t) {
      return t.type === "HoverTool";
    });
    hover.set('tooltips', {
      "index": "$index",
      "color": "$color[hex,swatch]:color",
      "radius": "@radius",
      "data (x, y)": "(@x, @y)",
      "cursor (x, y)": "($x, $y)",
      "(vx, vy)": "($vx, $vy)",
      "data (x, y)": "($x, $y)",
      "canvas (x, y)": "($sx, $sy)"
    });
    return Bokeh.Plotting.show(plot);
  });

}).call(this);

/*
//@ sourceMappingURL=hover.js.map
*/
