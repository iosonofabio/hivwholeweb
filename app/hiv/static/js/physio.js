function emptyPhysio(id, keepData) {

 var svg = d3.select('#'+id);
 svg.selectAll("*").remove();

 if ((typeof(keepData) == "undefined") | (!keepData)) svg.datum(null);

}

var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹",
    formatPower = function(d) { 
        var tmp='';
        if (d<0){tmp+='⁻';}
        return tmp+(Math.abs(d) + "").split("").map(function(c) { return superscript[c]; }).join(""); };

function updatePhysio(id, data) {

 var svg = d3.select('#'+id),
     divWidth = $('#'+id).parent().width();

 var margin = {top: 10, right: 80, bottom: 60, left: 80},
     width = 0.9 * divWidth - margin.left - margin.right,
     height = 330 - margin.top - margin.bottom;

 if (typeof(data) != "undefined") svg.datum(data);

 data = svg.datum();

 var colors = {"vl": "black", "cc": "steelblue"};

 var tmax = get_tmax(data);
 var vlmax = get_ymax(data.vl);
 var ccmax = get_ymax(data.cc);

 var y_vl = d3.scale.log()
      .domain([40, 1.1 * vlmax])
      .range([height, 0]);

 var y_cc = d3.scale.linear()
      .domain([0, 1.1 * ccmax])
      .range([height, 0]);

 var y = {"vl": y_vl, "cc": y_cc};
 
 var x = d3.scale.linear()
      .domain([0, 1.1 * tmax])
      .range([0, width]);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");
 
 var xAxisTop = d3.svg.axis()
     .scale(x)
     .orient("bottom")
     .tickFormat("");
 
 var yAxis_vl = d3.svg.axis()
     .scale(y_vl)
     .orient("left")
     .ticks(Math.round(Math.log(vlmax) / Math.LN10)-1, function(d) { return 10 + formatPower(Math.round(Math.log(d) / Math.LN10)); });

 var yAxis_cc = d3.svg.axis()
     .scale(y_cc)
     .ticks(5)
     .orient("right");

 svg.attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.bottom + margin.top);

 var chart = svg.append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 chart.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 50)
      .style("text-anchor", "middle")
      .text("Time from infection [days]");

 // Draw the x Grid lines
 chart.append("g")
     .attr("class", "grid")
     .attr("transform", "translate(0," + height + ")")
     .call(make_x_axis()
         .tickSize(-height, 0, 0)
         .tickFormat("")
     );

 chart.append("g")
      .attr("class", "d3-axis")
      .call(xAxisTop);

  
 var vlTextBox = chart.append("g")
      .attr("class", "d3-axis")
      .call(yAxis_vl)
      .append("g")
      .attr("class", "d3-axis-textbox")
      .attr("transform", "rotate(-90)");

 vlTextBox.append("text")
      .attr("dy", "-4.8em")
      .attr("x", -height / 2)
      .style("text-anchor", "middle")
      .text("Viral load [copies/ml]");
 
  vlTextBox.append("circle")
      .attr("class", "circles vlAxis")
      .attr("cy", "-5.1em")
      .attr("cx", -0.85 * height)
      .attr("r", 7)
      .style("fill", colors.vl);

 var ccTextBox = chart.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(" + width + " ,0)")
      .call(yAxis_cc)
      .append("g")
      .attr("class", "d3-axis-textbox")
      .attr("transform", "rotate(+90)");

 ccTextBox.append("text")
      .attr("dy", "-4em")
      .attr("x", height / 2)
      .style("text-anchor", "middle")
      .text("CD4+ counts [cells/ul]");

 ccTextBox.append("rect")
      .attr("class", "circles ccAxis")
      .attr("y", "-4.8em")
      .attr("x", 0.1 * height)
      .attr("width", 14)
      .attr("height", 14)
      .style("fill", colors.cc);

 chart.append("g")
      .attr("class", "circles VL")
      .selectAll()
      .data(data.vl)
      .enter()
      .append("circle")
      .attr("class", "circle")
      .attr("cx", function(d) { return x(d[0]); })
      .attr("cy", function(d) { return y_vl(d[1]); })
      .attr("r", 6)
      .style("fill", colors.vl);

 chart.append("g")
      .attr("class", "circles CC")
      .selectAll()
      .data(data.cc)
      .enter()
      .append("rect")
      .attr("class", "rect")
      .attr("x", function(d) { return x(d[0]) - 3; })
      .attr("y", function(d) { return y_cc(d[1]) - 3; })
      .attr("width", 12)
      .attr("height", 12)
      .style("fill", colors.cc);

 plotLine("vl");
 plotLine("cc");

 function plotLine(dtype, d) {
     chart.append("path")
         .attr("class", "data-line-" + dtype)
         .attr("d", lineFunction(data[dtype], y[dtype]))
	 .style("stroke-width", 3)
	 .style("stroke", colors[dtype])
	 .style("opacity", 0.6)
	 .style("fill", "none");
 }

 function lineFunction(data, yscale) {
  return d3.svg.line()
            .x(function(d) { return x(d[0]); })
            .y(function(d) { return yscale(d[1]); })
            .interpolate("linear")(data);
 }

 
 // function for the x grid lines
 function make_x_axis() {
     return d3.svg.axis()
         .scale(x)
         .orient("bottom");
 }

 function get_tmax(data) {
   var tvl = d3.max(data.vl, function(d) { return d[0]; });
   var tcc = d3.max(data.cc, function(d) { return d[0]; });
   return d3.max([tvl, tcc]);
 }
 
 function get_ymax(data) {
  return d3.max(data, function(d) { return d[1]; });
 }

}

