function emptyDivDivLocal(id, keepData) {

 var svg = d3.select('#'+id);
 svg.selectAll("*").remove();

 if ((typeof(keepData) == "undefined") | (!keepData)) svg.datum(null);

}

function updateDivDivLocal(id, data) {

    var svg = d3.select('#'+id),
        divWidth = $('#'+id).parent().width();

    var margin = {top: 10, right: 40, bottom: 150, left: 80},
        width = divWidth - margin.left - margin.right,
        height = 700 - margin.top - margin.bottom,
        height_genome = 200,
        vpad = 25,
        height_single = (height - height_genome - 2 * vpad) / 2;

    svg.attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.bottom + margin.top);

    var datal = data.divdiv.dg.length;
  
    var y = d3.scale.log()
         .domain([0.0001, 1.1])
         .range([height_single, 0]);
    
    var x = d3.scale.linear()
         .domain([-100, data.divdiv.dg[0][1].length * data.divdiv.dx + data.divdiv.block_length + 100])
         .range([0, width]);
    
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");
    
    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left")
        .ticks(4, function(d) { return 10 + formatPower(Math.round(Math.log(d) / Math.LN10)); });

    var yAxisRight = d3.svg.axis()
        .scale(y)
        .orient("right")
        .ticks(4, function(d) { return 10 + formatPower(Math.round(Math.log(d) / Math.LN10)); });
  
    var yGrid = d3.svg.axis()
        .scale(y)
        .orient("left")
        .tickSize(-width, 0, 0)
        .tickFormat("");

    var vis = svg.append("g")
         .attr("class", "d3-chart")
         .attr("id", "divdiv-chart"+id)
         .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var charts = {"dg": vis.append("g").attr("id", "divDivLocalDg-chart"+id),
                  "ds": vis.append("g").attr("id", "divDivLocalDs-chart"+id)
                        .attr("transform", "translate(0," + (height_single + vpad) + ")"),
                 };
       
    // Divergence Y axis
    charts.dg.append("g")
         .attr("class", "d3-axis")
         .call(yAxis)
         .append("text")
         .attr("transform", "rotate(-90)")
         .attr("dy", "-4.5em")
         .attr("x", -height_single / 2)
         .style("text-anchor", "middle")
         .text("Divergence");
  
    charts.dg.append("g")
         .attr("class", "d3-axis")
         .attr("transform", "translate(" + width + " ,0)")
         .call(yAxisRight);
  
    charts.dg.append("g")
         .attr("class", "grid")
         .call(yGrid);
  
    // Diversity Y axis
    charts.ds.append("g")
         .attr("class", "d3-axis")
         .call(yAxis)
         .append("text")
         .attr("transform", "rotate(-90)")
         .attr("dy", "-4.5em")
         .attr("x", -height_single / 2)
         .style("text-anchor", "middle")
         .text("Diversity");
    
    charts.ds.append("g")
         .attr("class", "d3-axis")
         .attr("transform", "translate(" + width + " ,0)")
         .call(yAxisRight);
  
    charts.ds.append("g")
         .attr("class", "grid")
         .call(yGrid);
  
  
    var lineFunc = d3.svg.line()
                         .x(function(d, i) { return x(i * data.divdiv.dx + 0.5 * data.divdiv.block_length); })
  		       .y(function(d) { return y(d); })
  		       .interpolate('monotone');
  
    var colors = d3.scale.linear()
    .domain([0, data.divdiv.len / 4, data.divdiv.len / 3, data.divdiv.len / 2,
             2 * data.divdiv.len / 3, 3 * data.divdiv.len / 4, data.divdiv.len])
    .interpolate(d3.interpolateRgb)
    .range(["darkblue", "blue", "cyan", "green", "yellow", "orange", "red"]);

    // Initial data
    charts.dg.append("svg:path")
          .attr("class", "lineDivergence")
          .attr("d", lineFunc(data.divdiv.dg[0][1]))
          .attr("stroke", "black")
          .attr("stroke-width", 2)
          .attr("fill", "none")
          .attr("opacity", 0.8);
   
    charts.ds.append("svg:path")
          .attr("class", "lineDiversity")
          .attr("d", lineFunc(data.divdiv.ds[0][1]))
          .attr("stroke", "black")
          .attr("stroke-width", 2)
          .attr("fill", "none")
          .attr("opacity", 0.8);
   
    // Slider
    var xsl = d3.scale.linear()
      .domain([1, datal])
      .range([0, width])
      .clamp(true);
   
    vis.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + (height + 70) + ")")
        .call(d3.svg.axis()
          .scale(xsl)
          .orient("bottom")
          .ticks(datal)
          .tickFormat(function(d, i){ return data.divdiv.dg[i][0]; })
          .tickPadding(12))
      .select(".domain")
      .select(function() { return this.parentNode.appendChild(this.cloneNode(true)); })
        .attr("class", "halo");
   
    vis.append("text")
        .attr("x", width / 2)
        .attr("y", height + 130)
        .style("text-anchor", "middle")
        .text("Time from infection [days]");
   
    var brush = d3.svg.brush()
        .x(xsl)
        .extent([1, 1])
        .on("brush", brushed);
   
    var slider = vis.append("g")
        .attr("class", "slider")
        .call(brush);
    
    slider.selectAll(".extent,.resize")
        .remove();
    
    var handle = slider.append("circle")
        .attr("class", "handle")
        .attr("transform", "translate(0," + (height + 70) + ")")
        .attr("r", 9);
    
    slider.call(brush.event)
        .transition() // gratuitous intro!
        .duration(750)
        .call(brush.extent([1, 1]))
        .call(brush.event);
   
    function get_itime(value) {
     return Math.min(datal - 1, Math.max(0, Math.round(value) - 1)); 
    }
    
    function brushed() {
      var value = brush.extent()[0];
    
      if (d3.event.sourceEvent) { // not a programmatic event
        value = xsl.invert(d3.mouse(this)[0]);
        brush.extent([value, value]);
      }
   
      var i_time = get_itime(value);
      handle.attr("cx", xsl(i_time + 1));
      charts.dg.selectAll(".lineDivergence")
          .transition()
          .duration(1000)
          .attr("d", lineFunc(data.divdiv.dg[i_time][1]));
   
      charts.ds.selectAll(".lineDiversity")
          .transition()
          .duration(1000)
          .attr("d", lineFunc(data.divdiv.ds[i_time][1]));
   
    }
 
    // update the genome
    var gChart = genomeChart().resizeSvg(false)
        .drawBorderTop(false)
        .margin({left: margin.left, right: margin.right, bottom: margin.bottom,
                 top: margin.top + height - height_genome})
        .vpadBlockTop(vpad)
        .width(width)
        .height(height_genome)
        .x(x);
        //.zoomCallbacks({'zoomin': {'middle': zoomIn},
        //                'zoomout': {'post': zoomOut}});

    svg.datum(data)
        .call(gChart);

}
