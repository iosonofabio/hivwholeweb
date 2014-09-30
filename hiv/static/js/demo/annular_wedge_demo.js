(function() {
  require(['main'], function(Bokeh) {
    var N, color, data, i, options, plot, r, rects, val, x, xs, ys;
    N = 630;
    N = 10;
    r = new Bokeh.Random(123456789);
    xs = (function() {
      var _i, _len, _ref, _results;
      _ref = _.range(N);
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        x = _ref[_i];
        _results.push(x / 50);
      }
      return _results;
    })();
    ys = (function() {
      var _i, _len, _results;
      _results = [];
      for (_i = 0, _len = xs.length; _i < _len; _i++) {
        x = xs[_i];
        _results.push(Math.sin(x));
      }
      return _results;
    })();
    color = (function() {
      var _i, _len, _results;
      _results = [];
      for (_i = 0, _len = ys.length; _i < _len; _i++) {
        val = ys[_i];
        _results.push("rgb(" + (Math.floor(155 + 100 * val)) + ", " + (Math.floor(100 + 50 * val)) + ", " + (Math.floor(150 - 50 * val)) + ")");
      }
      return _results;
    })();
    xs = (function() {
      var _i, _len, _ref, _results;
      _ref = _.range(N);
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        i = _ref[_i];
        _results.push(r.randf() * 2);
      }
      return _results;
    })();
    ys = (function() {
      var _i, _len, _ref, _results;
      _ref = _.range(N);
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        i = _ref[_i];
        _results.push(r.randf() * 2);
      }
      return _results;
    })();
    data = {
      x: xs,
      y: ys,
      color: color
    };
    rects = {
      type: 'annular_wedge',
      x: 'x',
      y: 'y',
      size: .2,
      inner_radius: .1,
      outer_radius: .20,
      start_angle: .5,
      end_angle: 1.25,
      line_alpha: 1,
      fill_color: '#0F0FF0',
      do_fill: true
    };
    options = {
      title: "Annular Wedge Demo",
      dims: [800, 500],
      xaxes: "below",
      yaxes: "left",
      tools: "pan,zoom,resize,preview,select",
      legend: false
    };
    plot = Bokeh.Plotting.make_plot(rects, data, options);
    return Bokeh.Plotting.show(plot);
  });

}).call(this);

/*
//@ sourceMappingURL=annular_wedge_demo.js.map
*/