define(function(){
  var template = function(__obj) {
  var _safe = function(value) {
    if (typeof value === 'undefined' && value == null)
      value = '';
    var result = new String(value);
    result.ecoSafe = true;
    return result;
  };
  return (function() {
    var __out = [], __self = this, _print = function(value) {
      if (typeof value !== 'undefined' && value != null)
        __out.push(value.ecoSafe ? value : __self.escape(value));
    }, _capture = function(callback) {
      var out = __out, result;
      __out = [];
      callback.call(this);
      result = __out.join('');
      __out = out;
      return _safe(result);
    };
    (function() {
      _print(_safe('<div class=\'bokeh_canvas_wrapper_outer\'>\n\t<table>\n\t\t<tr>\n\t\t\t<td></td>\n\t\t\t<td class=\'bk-plot-above\'></td>\n\t\t\t<td></td>\n\t\t</tr>\n\t\t<tr>\n\t\t\t<td class="bk-plot-left"></td>\n\t\t\t<td class=\'bk-plot-canvas-wrapper\'></td>\n\t\t\t<td class="bk-plot-right"></td>\n\t\t</tr>\n\t\t<tr>\n\t\t\t<td></td>\n\t\t\t<td class=\'bk-plot-below\'></td>\n\t\t\t<td></td>\n\t\t</tr>\n\t</table>\n \n\n</div>\n'));
    
    }).call(this);
    
    return __out.join('');
  }).call((function() {
    var obj = {
      escape: function(value) {
        return ('' + value)
          .replace(/&/g, '&amp;')
          .replace(/</g, '&lt;')
          .replace(/>/g, '&gt;')
          .replace(/"/g, '&quot;');
      },
      safe: _safe
    }, key;
    for (key in __obj) obj[key] = __obj[key];
    return obj;
  })());
};
  return template;
});
