function updateAlleleFrequencies(data) {
    var datal = data.len,
        id = data.id;
   
    var div_width = $('.svg-container').width();
   
    var margin = {top: 45, right: 50, bottom: 50, left: 80},
        width = div_width - margin.left - margin.right,
        height = $('.'+id).height() - margin.top - margin.bottom,
        vpad = 50, height_genome=200,
        heightAft = (height - height_genome - vpad - vpad) / 2,
        currentTime = data.times[0],
        zoomArea = [-1, 10000];
   
    var svg = d3.select("."+id)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.bottom + margin.top);
   
    var vis = svg.append("g")
         .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
   
    var visAft = vis.append("g")
         .attr("class", "d3-chart aft-chart"+id);
   
    var visAfg = vis.append("g")
         .attr("class", "d3-chart afg-chart"+id)
         .attr("transform", "translate(0," + (heightAft + 1.5 * vpad) + ")");
   
    var visGenome = vis.append("g")
             .attr("class", "d3-chart genome-chart"+id)
             .attr("transform", "translate(0," + (2 * heightAft + 2 * vpad) + ")");
   
    var y = d3.scale.log()
         .domain([0.0009, 0.9991])
         .range([heightAft, 0])
         .clamp(true);
   
    var x = d3.scale.linear()
         .domain([-0.05 * data.genome.len, 1.05 * data.genome.len])
         .range([0, width]);
    
    var xT = d3.scale.linear()
         .domain([-10, 1.05 * data.tmax])
         .range([0, width])
         .clamp(true);
   
    var xTSlider = d3.scale.linear()
        .domain([data.times[0], data.times[data.times.length - 1]])
        .range([xT(data.times[0]), xT(data.times[data.times.length - 1])])
        .clamp(true);
   
    var xAxis = d3.svg.axis()
        .scale(xT)
        .orient("bottom");
    
    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");
   
    var yAxisRight = d3.svg.axis()
        .scale(y)
        .orient("right");
   
    visAft.append("g")
         .attr("class", "d3-axis")
         .attr("transform", "translate(0," + heightAft + ")")
         .call(xAxis)
         .append("text")
         .attr("x", width / 2)
         .attr("y", 40)
         .style("text-anchor", "middle")
         .text("Time from infection [days]");
     
    visAft.append("g")
         .attr("class", "d3-axis")
         .call(yAxis)
         .append("text")
         .attr("transform", "rotate(-90)")
         .attr("dy", "-4.5em")
         .attr("x", -(heightAft) / 2)
         .style("text-anchor", "middle")
         .text("Frequency");
   
    visAft.append("g")
         .attr("class", "d3-axis")
         .attr("transform", "translate(" + width + ",0)")
         .call(yAxisRight);
   
    var aftLine = d3.svg.line()
        .x(function(d) { return xT(d[0]); })
        .y(function(d) { return y(d[1]); })
        .interpolate('monotone');
   
    var colors = d3.scale.linear()
        .domain([0, data.len / 4, data.len / 3,
                data.len / 2, 2 * data.len / 3,
                3 * data.len / 4, data.len])
        .interpolate(d3.interpolateRgb)
        .range(["darkblue", "blue", "cyan",
               "green", "yellow",
               "orange", "red"]);
   
    visAft.selectAll(".aft")
         .data(data.aft)
         .enter()
         .append("svg:path")
         .attr("class", "aft")
         .attr("id", function(d, i) { return "aft-" + (i+1); })
         .attr("d", function(d) { return aftLine(d[2]); })
         .attr("stroke", function(d) { return colors(d[0]); })
         .attr("stroke-width", 2)
         .attr("fill", "none")
         .attr("opacity", 0.5);
   
    // Add number of templates line if present
    visAft.append("svg:path")
      .attr("class", "ntemplates")
      .attr("d", d3.svg.line()
                   .x(function(d) { return xT(d[0]); })
   		.y(function(d) { return y(1.0 / d[1]); })
   		.interpolate("monotone")(data.ntemplates)
   		)
      .attr("stroke", "darkred")
      .attr("stroke-width", 15)
      .attr("fill", "none")
      .attr("opacity", 0.3);
   
    // Add max depth for sequencing errors
    visAft.append("line")
        .attr("class", "maxDepth")
        .attr("x1", xT(0))
        .attr("x2", x.range()[1])
        .attr("y1", y(2e-3))
        .attr("y2", y(2e-3))
        .attr("stroke", "steelblue")
        .attr("stroke-width", 15)
        .attr("opacity", 0.4);


    // update the allele frequencies along the genome
    visAfg.datum(data)
        .call(chartAlleleFrequencyGenome);

    // add reusable, dynamic genome chart
    var gChart = genomeChart().resizeSvg(false)
        .drawBorderTop(false)
        .margin({left: margin.left, right: margin.right, bottom: margin.bottom,
            top: margin.top + 2 * heightAft + 1.5 * vpad})
        .vpadBlockTop(vpad / 2)
        .width(width)
        .height(height_genome)
        .x(x)
        .zoomCallbacks({'zoomin': {'middle': zoomIn},
                        'zoomout': {'post': zoomOut}});

    svg.datum(data)
        .call(gChart);

    function zoomIn(zoomData) {
        zoomArea = [zoomData.start, zoomData.end];
        updatePlots();
    }

    function zoomOut() {
        zoomArea = [-1, 10000];
        updatePlots();    
    }

    function updatePlots() {
        // remove and recreate allele frequencies along the genome
        vis.select(".afg")
            .selectAll(".afg_point")
            .remove();

        vis.select(".afg")
            .selectAll(".afg_point")
            .data(data.aft)
            .enter()
            .append("circle")
            .filter(function(d) {
                return (hasTime(d[2], currentTime)) &
                    (d[0] >= zoomArea[0]) & (d[0] < zoomArea[1]);
            })
            .attr("class", "circle afg_point")
            .attr("cx", function(d) { return x(d[0]); })
            .attr("cy", function(d) { return y(getFreq(d[2], currentTime)); })
            .attr("r", 3)
            .attr("fill", function(d) { return colors(d[0]); })
            .attr("stroke", "grey")
            .attr("stroke-width", 0.3);

        // make only lines in the zoom area visible
        vis.selectAll(".aft")
            .attr("opacity", function(daft) {
                if ((daft[0] >= zoomArea[0]) & (daft[0] < zoomArea[1]))
                    return 0.5;
                else
                    return 0;
            });
    }

    function chartAlleleFrequencyGenome(selection) {
        selection.each(function (data) {
       
            var y = d3.scale.log()
                 .domain([0.0009, 0.9991])
                 .range([heightAft, 0])
                 .clamp(true);
       
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
       
            var vis = d3.select(this);
       
            // add slider
            vis.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + (-heightAft - 1.5 * vpad - 20) + ")")
                .call(d3.svg.axis()
                .scale(xTSlider)
                .orient("top"))
                .select(".domain")
                .select(function() { return this.parentNode.appendChild(this.cloneNode(true)); })
                .attr("class", "halo");
       
            var brush = d3.svg.brush()
                .x(xTSlider)
                .extent([1, 1])
                .on("brush", brushed);
       
            var slider = vis.append("g")
               .attr("class", "slider")
               .call(brush);
            
            slider.selectAll(".extent,.resize")
                .remove();
       
            var handle = slider.append("circle")
                .attr("class", "handle")
                .attr("transform", "translate(0," + (- 1.5 * vpad - heightAft - 20) + ")")
                .attr("r", 5);
       
            slider.call(brush.event)
                .call(brush.extent([1, 1]))
                .call(brush.event);
       
            // NOTE: we should make transitions here, but this involves enumerating the cases and is therefore very tedious
            function brushed() {
                var value = brush.extent()[0];
            
                if (d3.event.sourceEvent) { // not a programmatic event
                    value = xTSlider.invert(d3.mouse(this)[0]);
                    brush.extent([value, value]);
                }
       
                currentTime = getTime(data.times, value);
                handle.attr("cx", xTSlider(currentTime));
       
                // update time tracker
                vis.select("#aftTimeTracker")
                    .attr("x1", xTSlider(currentTime))
                    .attr("x2", xTSlider(currentTime));
       
                updatePlots();
       
                // update ntemplates
                vis.select("#afgTemplates")
                    .remove();
             
                if (hasTime(data.ntemplates, currentTime)) {
                    vis.append("svg:line")
                        .attr("id", "afgTemplates")
                        .datum(1.0 / getFreq(data.ntemplates, currentTime))
                        .attr("x1", xT.range()[0])
                        .attr("x2", xT.range()[1])
                        .attr("y1", function(d) { return y(d); })
                        .attr("y2", function(d) { return y(d); })
                        .attr("stroke", "darkred")
                        .attr("stroke-width", 15)
                        .attr("fill", "none")
                        .attr("opacity", 0.3);
                }
       
            }
            // end of slider
       
            vis.append("g")
                .attr("class", "grid")
                .call(yAxisGrid);
       
            vis.append("g")
                .attr("class", "d3-axis")
                .call(yAxis)
                .append("text")
                .attr("transform", "rotate(-90)")
                .attr("dy", "-4.5em")
                .attr("x", -(heightAft) / 2)
                .style("text-anchor", "middle")
                .text("Frequency");
       
            vis.append("g")
                .attr("class", "d3-axis")
                .attr("transform", "translate(" + width + ",0)")
                .call(yAxisRight);
       
            // Initial data
            vis.append("g")
                 .attr("class", "circles afg")
                 .selectAll(".afg_point")
                 .data(data.aft)
                 .enter()
                 .append("circle")
                 .filter(function(d) { return hasTime(d[2], currentTime); })
                 .attr("class", "circle afg_point")
                 .attr("cx", function(d) { return x(d[0]); })
                 .attr("cy", function(d) { return y(d[2][0][1]); })
                 .attr("r", 3)
                 .attr("fill", function(d) { return colors(d[0]); })
                 .attr("stroke", "grey")
                 .attr("stroke-width", 0.3);
       
            // number of templates if appropriate
            if (hasTime(data.ntemplates, currentTime)) {
             vis.append("svg:line")
                .attr("id", "afgTemplates")
                .datum(1.0 / getFreq(data.ntemplates, currentTime))
                .attr("x1", xT.range()[0])
                .attr("x2", xT.range()[1])
                .attr("y1", function(d) { return y(d); })
                .attr("y2", function(d) { return y(d); })
                .attr("stroke", "darkred")
                .attr("stroke-width", 15)
                .attr("fill", "none")
                .attr("opacity", 0.3);
            }
       
             // add sequencing depth
             vis.append("svg:line")
                .attr("id", "afgSeqDepth")
                .datum(0.002)
                .attr("x1", xT.range()[0])
                .attr("x2", xT.range()[1])
                .attr("y1", function(d) { return y(d); })
                .attr("y2", function(d) { return y(d); })
                .attr("stroke", "steelblue")
                .attr("stroke-width", 15)
                .attr("fill", "none")
                .attr("opacity", 0.4);
       
       
            // vertical line for time tracking
            vis.append("svg:line")
                .attr("id", "aftTimeTracker")
                .attr("x1", xTSlider(0))
                .attr("x2", xTSlider(0))
                .attr("y1", -heightAft - 1.5 * vpad - 15)
                .attr("y2", - 1.5 * vpad)
                .attr("stroke", "grey")
                .attr("stroke-width", 5)
                .attr("opacity", 0.5);
    
        });
    }

    function getTime(times, value) {
        var t = times[0],
            d = Math.abs(value - times[0]),
            dnew;
        for (i = 1; i < times.length; i++) {
            dnew = Math.abs(value - times[i]);
            if (dnew < d) {
                t = times[i];
                d = dnew;
            }
        }
        return t; 
    }
     
    function hasTime(arr, time) {
        for (i = 0; i < arr.length; i++)
            if (arr[i][0] == time) return true;
        return false;
    }
    
    function getFreq(arr, time) {
        for (i = 0; i < arr.length; i++)
            if (arr[i][0] == time) return arr[i][1];
        return undefined;
    }

}
