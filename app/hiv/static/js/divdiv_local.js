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


//    // FIXME: this is done in a reusable chart now!
//    update_genome(data_gn, id, x);
//
//    function update_genome(data, id, x) {
//    
//     var tip = d3.tip()
//          .attr('class', 'd3-tip')
//          .html(function(d) { return d.name + ": " + (+d.location[0][0] + 1) + ", " + d.location[0][1]; });
//    
//      var chart = d3.select(".genome-chart"+id);
//      var charty = +(chart.attr("transform").split(",")[1].split(")")[0]);
//    
//      var chartdd = d3.select("#divdiv-chart"+id);
//      
//      // add tooltip
//      chart.call(tip);
//     
//      var div_width = $('.svg-container').width();
//    
//      var height_block = 20, vpad_block = 5;
//    
//      function get_n_feature_types(data) {
//       features = [];
//       var found = 0;
//       for(i = 0; i < data.features.length; i++) {
//        found = 0;
//        for(j = 0; j < features.length; j++) {
//         if (features[j] == data.features[i].type) {
//          found = 1;
//         }
//        }
//        if (found == 0) {
//         features.push(data.features[i].type);
//        }
//       }
//       return features.length; 
//      }
//    
//      function get_feature_group(groupname) {
//       var group = [];
//       var gname = "";
//       var tmpfea;
//       for(i = 0; i < data.features.length; i++) {
//        gname = data.features[i].type
//        if (gname == groupname) {
//         group.push(data.features[i]);
//        } else if (gname + '-single' == groupname) {
//         if (['vpu', 'vpr', 'nef', 'vif'].indexOf(data.features[i].name) != -1) {
//          group.push(data.features[i]);
//         }
//        } else if (gname + '-exon1' == groupname) {
//         if (['rev', 'tat'].indexOf(data.features[i].name) != -1) {
//          tmpfea = {"name": data.features[i].name + ' exon 1',
//                    "type": "gene",
//                    "location": [data.features[i].location[0]]};
//          group.push(tmpfea);
//         }
//        } else if (gname + '-exon2' == groupname) {
//         if (['rev', 'tat'].indexOf(data.features[i].name) != -1) {
//          tmpfea = {"name": data.features[i].name + ' exon 2',
//                    "type": "gene",
//                    "location": [data.features[i].location[1]]};
//          group.push(tmpfea);
//         }
//        }
//       }
//       return group;
//      }
//    
//      function height_feature(feature) {
//       if (feature.type == "fragment") {
//        return height_fragment(feature);
//       } else if ((feature.type == "gene") | (feature.type == "protein")) {
//        return height_inframe(feature); 
//       } else {
//        return 0;
//       }
//      }
//    
//      function height_fragment(fragment) {
//       var fn = +(fragment.name[1]);
//       var height = vpad_block;
//       if ((fn % 2) == 0) {
//        height += height_block + vpad_block;
//       }
//       return height;
//      }
//    
//      function height_inframe(feature) {
//       var frame;
//       if (['tat exon 2', 'rev exon 2'].indexOf(feature.name) != -1) {
//         frame = (+(feature.location[0][1]) - (+data.framestart)) % 3;
//       } else {
//         frame = (+(feature.location[0][0]) - (+data.framestart)) % 3;
//       } 
//       return frame * (height_block + vpad_block);
//      }
//    
//      function moverFeature(d) {
//        // find position of this box (FIXME: do it better)
//        var boxy = +(d3.select(this).attr("transform").split(",")[1].split(")")[0]);
//    
//        // change color
//        var rect = d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect");
//    
//        rect.style("fill", "darkred");  
//    
//        // show tooltip
//        tip.show(d);
//    
//        // show line connecting the exons
//        var exons = [];
//        if (d.name.indexOf('exon 1') != -1) {
//         exons = ["exon-1", "exon-2"];
//        } else if (d.name.indexOf('exon 2') != -1) {
//         exons = ["exon-2", "exon-1"];
//        }
//        if (exons.length > 0) {
//          var rect2 = d3.select("#" + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "rect");
//          var d2 = rect2.data()[0];
//          var box2 = d3.select('#' + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "-box");
//          // FIXME: rewrite (?): we extract the y coordinate from the affine transformation STRING!
//          var box2y = +(box2.attr("transform").split(",")[1].split(")")[0]);
//    
//          chart.append('line')
//           .attr("class", "exon-line")
//           .attr("x1", 0.5 * (x(d.location[0][1]) + x(d.location[0][0])))
//           .attr("x2", 0.5 * (x(d2.location[0][1]) + x(d2.location[0][0])))
//           .attr("y1", boxy + 0.5 * height_block)
//           .attr("y2", box2y + 0.5 * height_block)
//           .style("opacity", 0.7)
//           .style("stroke", "darkred")
//           .style("stroke-width", 2);
//        }
//    
//       // show vertical "grid" lines at the edges of the region
//       chartdd.selectAll()
//            .data([d.location[0][0], d.location[0][1]])
//    	.enter()
//    	.append("line")
//            .attr("class", "grid feature-edge-line")
//            .attr("x1", function(dd) { return x(dd); })
//            .attr("x2", function(dd) { return x(dd); })
//            .attr("y1", 0)
//            .attr("y2", charty + boxy)
//    	.style("fill", "none")
//    	.style("stroke", "grey")
//    	.style("stroke-width", 1);
//    
//      }
//    
//      function moutFeature(d) {
//        // change color back
//        d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect")
//         .style("fill", "steelblue");
//    
//        // hide tooltip
//        tip.hide(d);   
//    
//        // hide exon line
//        if ((d.name.indexOf('exon 1') != -1) | (d.name.indexOf('exon 2') != -1)) {
//          chart.selectAll(".exon-line").remove();
//        }
//    
//        // hide vertical "grid" lines
//       chartdd.selectAll(".feature-edge-line")
//              .remove();
//      }
//    
//      function plot_group(groupname, dy) {
//       var group = get_feature_group(groupname);
//       var fea = chart.selectAll()
//                      .data(group)
//            	   .enter()
//            	   .append("g")
//                      .attr("class", "featurebox " + groupname)
//           	   .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "-box"; })
//            	  .attr("transform", function(d) { return "translate(" + x(d.location[0][0]) + "," + (dy + height_feature(d)) + ")";})
//            .on('mouseover', moverFeature)
//            .on('mouseout', moutFeature);
//    
//       var fearect = fea.append("rect")
//            .attr("class", "featurerect")
//            .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "rect"; })
//            .attr("x", 0)
//            .attr("y", 0)
//            .attr("width", function(d) { return x(d.location[0][1]) - x(d.location[0][0]); })
//            .attr("height", 20)
//            .style("fill", "steelblue");
//
//    
//       // show text only of longer things (that's why we have a tooltip)
//       fea.append("text")
//              .attr("x", function(d) { return 0.5 * (x(d.location[0][1]) - x(d.location[0][0])); })
//              .attr("y", 15)
//              .text(function(d) {
//                if ((d.location[0][1] - d.location[0][0]) > 350) {
//                  return d.name;
//                } else {
//                  return "";
//                }
//              });
//       }
//    
//      plot_group("fragment", 0);
//      plot_group("protein", 75);
//      plot_group("gene-single", 75);
//      plot_group("gene-exon1", 75);
//      plot_group("gene-exon2", 75);
//      plot_group("RNA_structure", 170);
//      plot_group("other", 170);
//    
//    }
 
}

