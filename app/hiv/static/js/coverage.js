function updateCoverage(data) {
    var id = data.id;
    var div_width = $('.svg-container').width();
   
    var margin = {top: 10, right: 60, bottom: 150, left: 80},
        width = 0.9 * div_width - margin.left - margin.right,
        vpad = 20, height_genome = 200,
        height = $('.'+id).height() - margin.top - margin.bottom,
        height_cov = height - vpad - height_genome,
        currentITime = 0;
   
    var datal = data.cov.length;
    var covmax = get_ymax(data.cov);
   
    var y = d3.scale.log()
         .domain([0.8, covmax])
         .range([height_cov, 0]);
    
    var x = d3.scale.linear()
         .domain([-0.05 * data.cov[0][1].length, 1.05 * data.cov[0][1].length])
         .range([0, width]);
    
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");
    
    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");
   
    var yAxisRight = d3.svg.axis()
        .scale(y)
        .orient("right");
   
    var yAxisGrid = d3.svg.axis()
         .scale(y)
         .orient("right")
         .tickSize(width, 0, 0)
         .tickFormat("");
   
    var svg = d3.select("."+id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.bottom + margin.top);
   
    var chart = svg.append("g")
         .attr("id", "main-chart"+id)
         .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
   
    var charts = {"cov": chart.append("g")
                      .attr("id", "cov-chart"+id),
                  };
   
    // Draw the grid lines (they should stay behind)
    charts.cov.append("g")
        .attr("class", "grid")
        .call(yAxisGrid);
        
    charts.cov.append("g")
         .attr("class", "d3-axis")
         .call(yAxis)
         .append("text")
         .attr("transform", "rotate(-90)")
         .attr("dy", "-4.5em")
         .attr("x", -height_cov / 2)
         .style("text-anchor", "middle")
         .text("Coverage");
   
    charts.cov.append("g")
         .attr("class", "d3-axis")
         .attr("transform", "translate(" + width + " ,0)")
         .call(yAxisRight);
   
    // Initial data
    charts.cov.append("g")
         .attr("class", "circles COV")
         .selectAll("cov_point")
         .data(data.cov[currentITime][1])
         .enter()
         .append("circle")
         .attr("class", "circle cov_point")
         .attr("cx", function(d, i) { return x(i); })
         .attr("cy", function(d) { return y(d); })
         .attr("r", 3);
   
    // Add number of templates line if present
    charts.cov.selectAll(".ntemplates")
      .data(data.ntemplates)
      .enter()
      .append("line")
      .filter(function(d) { return d[0] == data.cov[0][0]; })
      .attr("class", "ntemplates")
      .attr("x1", x(0))
      .attr("x2", function(d) { return x(data.cov[0][1].length); })
      .attr("y1", function(d) { return y(d[1]); })
      .attr("y2", function(d) { return y(d[1]); })
      .attr("stroke", "darkred")
      .attr("stroke-width", 15)
      .attr("opacity", 0.3);
   
    // Slider
    var xsl = d3.scale.linear()
      .domain([1, datal])
      .range([0, width])
      .clamp(true);
   
    chart.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + (height + 70) + ")")
        .call(d3.svg.axis()
          .scale(xsl)
          .orient("bottom")
          .ticks(datal)
          .tickFormat(function(d, i){ return data.cov[i][0]; })
          .tickPadding(12))
      .select(".domain")
      .select(function() { return this.parentNode.appendChild(this.cloneNode(true)); })
        .attr("class", "halo");
   
    chart.append("text")
        .attr("x", width / 2)
        .attr("y", height + 130)
        .style("text-anchor", "middle")
        .text("Time from infection [days]");
   
    var brush = d3.svg.brush()
        .x(xsl)
        .extent([1, 1])
        .on("brush", brushed);
   
    var slider = chart.append("g")
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
   
      currentITime = get_itime(value);
      handle.attr("cx", xsl(currentITime + 1));
      // we could put a transition here, but it the code must be optimized first, or else it's
      // going to stutter a lot
      charts.cov.selectAll(".cov_point")
         .data(data.cov[currentITime][1])
         .attr("cy", function(d) { return y(d); });
   
      // add template number line
      charts.cov.selectAll(".ntemplates").remove();
      charts.cov.selectAll(".ntemplates")
        .data(data.ntemplates)
        .enter()
        .append("line")
        .filter(function(d) { return d[0] == data.cov[currentITime][0]; })
        .attr("class", "ntemplates")
        .attr("x1", x(0))
        .attr("x2", function(d) { return x(data.cov[currentITime][1].length); })
        .attr("y1", function(d) { return y(d[1]); })
        .attr("y2", function(d) { return y(d[1]); })
        .attr("stroke", "darkred")
        .attr("stroke-width", 15)
        .attr("opacity", 0.3);
   
    }
   
    function zoomIn(zoomData) {

        // delete missing points
        charts.cov.selectAll(".cov_point")
            .remove();

        // rescale the other ones
        charts.cov.select(".COV")
            .selectAll(".cov_point")
            .data(data.cov[currentITime][1].slice(zoomData.start, zoomData.end))
            .enter()
            .append("circle")
            .attr("class", "circle cov_point")
            .attr("cy", function(d) { return y(d); })
            .attr("cx", function(d, i) { return x(i + zoomData.start); })
            .attr("r", 3);
    }

    function zoomOut() {
        // delete all points and restart (could be done better)
        charts.cov.selectAll(".cov_point")
            .remove();

        charts.cov.select(".COV")
             .selectAll(".cov_point")
             .data(data.cov[currentITime][1])
             .enter()
             .append("circle")
             .attr("class", "circle cov_point")
             .attr("cx", function(d, i) { return x(i); })
             .attr("cy", function(d) { return y(d); })
             .attr("r", 3);
    
    }

    // add reusable, dynamic genome chart
    var gChart = genomeChart().resizeSvg(false)
        .drawBorderTop(false)
        .margin({left: margin.left, right: margin.right, bottom: margin.bottom,
            top: margin.top + height_cov})
        .vpadBlockTop(vpad)
        .width(width)
        .height(height_genome)
        .x(x)
        .zoomCallbacks({'zoomin': {'middle': zoomIn},
                        'zoomout': {'post': zoomOut}});

    svg.datum(data)
        .call(gChart);


    function get_ymax(data) {
        var ymax = d3.max(data[0][1]);
        var ytmp;
        for (i = 1; i < data.length; i++) {
            ytmp = d3.max(data[i][1]);
            if (ytmp > ymax)
                ymax = ytmp;
        }
        return ymax;
    }


} 
