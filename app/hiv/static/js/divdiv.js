 function get_tmax(data) {
  var tmax = 0;
  for (i = 0; i < data.dg.length; i++) {
    if (data.dg[i][0] > tmax) {
     tmax = data.dg[i][0];
    }
  }
  for (i = 0; i < data.ds.length; i++) {
    if (data.ds[i][0] > tmax) {
     tmax = data.ds[i][0];
    }
  }
  return tmax;
 }

 function get_ymax(data) {
  var ymax = data[0][1];
  for (i = 1; i < data.length; i++) {
    if (data[i][1] > ymax) {
     ymax = data[i][1];
    }
  }
  return ymax;
 }

 function get_ymin(data) {
  var ymin = data[0][1];
  for (i = 1; i < data.length; i++) {
    if (data[i][1] < ymin) {
     ymin = data[i][1];
    }
  }
  return ymin;
 }

 function update(data, id) {

 var colors = {"dg": "black", "ds": "steelblue"};

 var div_width = $('.svg-container').width();


  var margin = {top: 10, right: 80, bottom: 60, left: 80},
      width = 0.9 * div_width - margin.left - margin.right,
      height = 500 - margin.top - margin.bottom;

  var tmax = get_tmax(data);
  var dgmax = get_ymax(data.dg);
  var dgmin = get_ymin(data.dg);
  var dsmax = get_ymax(data.ds);
  var dsmin = get_ymin(data.ds);

  var y = d3.scale.log()
       .domain([1e-3, 0.5])
       .range([height, 0]);

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

 var xAxisGrid = d3.svg.axis()
      .scale(x)
      .orient("bottom")
      .tickSize(-height, 0, 0)
      .tickFormat("");
 
  var yAxis_dg = d3.svg.axis()
      .scale(y)
      .orient("left");

  var yAxis_ds = d3.svg.axis()
      .scale(y)
      .orient("right");

 var yAxisGrid = d3.svg.axis()
      .scale(y)
      .orient("right")
      .tickSize(width, 0, 0)
      .tickFormat("");

  var chart_ext = d3.select("."+id)
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.bottom + margin.top);

  var chart = chart_ext.append("g")
       .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 // Draw the grid lines (they should stay behind)
 chart.append("g")
     .attr("class", "grid")
     .attr("transform", "translate(0," + height + ")")
     .call(xAxisGrid);

 chart.append("g")
     .attr("class", "grid")
     .call(yAxisGrid);

 chart.append("g")
       .attr("class", "d3-axis")
       .attr("transform", "translate(0," + height + ")")
       .call(xAxis)
       .append("text")
       .attr("x", width / 2)
       .attr("y", 40)
       .style("text-anchor", "middle")
       .text("Time from infection [days]");
   
 chart.append("g")
      .attr("class", "d3-axis")
      .call(xAxisTop);
 
 var dgTextBox = chart.append("g")
       .attr("class", "d3-axis")
       .call(yAxis_dg)
       .append("g")
       .attr("class", "d3-axis-textbox")
       .attr("transform", "rotate(-90)");

 dgTextBox.append("text")
       .attr("dy", "-4.7em")
       .attr("x", -height / 2)
       .style("text-anchor", "middle")
       .text("Divergence");

 dgTextBox.append("circle")
      .attr("class", "circles dgAxis")
      .attr("cy", "-5.1em")
      .attr("cx", -0.65 * height)
      .attr("r", 7)
      .style("fill", colors.dg);

 var dsTextBox = chart.append("g")
       .attr("class", "d3-axis")
       .attr("transform", "translate(" + width + " ,0)")
       .call(yAxis_ds)
       .append("g")
       .attr("class", "d3-axis-textbox")
       .attr("transform", "rotate(+90)");

 dsTextBox.append("text")
       .attr("dy", "-5em")
       .attr("x", height / 2)
       .style("text-anchor", "middle")
       .text("Diversity");

 dsTextBox.append("rect")
      .attr("class", "circles dsAxis")
      .attr("y", "-5.8em")
      .attr("x", 0.35 * height)
      .attr("width", 14)
      .attr("height", 14)
      .style("fill", colors.ds);

  chart.append("g")
       .attr("class", "circles DG")
       .selectAll()
       .data(data.dg)
       .enter()
       .append("circle")
       .attr("class", "circle")
       .attr("cx", function(d) { return x(d[0]); })
       .attr("cy", function(d) { return y(d[1]); })
       .attr("r", 6);

  chart.append("g")
       .attr("class", "circles DS")
       .selectAll()
       .data(data.ds)
       .enter()
       .append("rect")
       .attr("class", "rect")
       .attr("x", function(d) { return x(d[0]) - 3; })
       .attr("y", function(d) { return y(d[1]) - 3; })
       .attr("width", 12)
       .attr("height", 12)
       .style("fill", "steelblue");

 }

