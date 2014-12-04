function get_ymax(data) {
 var ymax = d3.max(data[0][1]);
 var ytmp;
 for (i = 1; i < data.length; i++) {
   ytmp = d3.max(data[i][1]);
   if (ytmp > ymax) {
    ymax = ytmp;
   }
 }
 return ymax;
}

function update(data, id) {

 var div_width = $('.svg-container').width();

 var margin = {top: 10, right: 60, bottom: 150, left: 80},
     width = 0.9 * div_width - margin.left - margin.right,
     height_cov = 250, vpad = 20, height_genome = 200,
     height = $('.'+id).height() - margin.top - margin.bottom,
     height_cov = height - vpad - height_genome;

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

 var chart_ext = d3.select("."+id)
     .attr("width", width + margin.left + margin.right)
     .attr("height", height + margin.bottom + margin.top);

 var chart = chart_ext.append("g")
      .attr("id", "main-chart"+id)
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var charts = {"cov": chart.append("g")
                   .attr("id", "cov-chart"+id),
               "genome": chart.append("g")
                   .attr("class", "d3-chart genome-chart"+id)
                   .attr("transform", "translate(0," + (height_cov + vpad) + ")")};

 // Draw the grid lines (they should stay behind)
 charts.cov.append("g")
     .attr("class", "grid")
     .call(yAxisGrid);

 charts.cov.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 40)
      .style("text-anchor", "middle")
      .text("Position [bp]");
  
 charts.cov.append("g")
      .attr("class", "d3-axis")
      .call(yAxis)
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("dy", "-4.5em")
      .attr("x", -height / 2)
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
      .data(data.cov[0][1])
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

   var i_time = get_itime(value);
   handle.attr("cx", xsl(i_time + 1));
   // we could put a transition here, but it the code must be optimized first, or else it's
   // going to stutter a lot
   charts.cov.selectAll(".cov_point")
      .attr("cy", function(d, i) { return y(data.cov[i_time][1][i]); });

   // add template number line
   charts.cov.selectAll(".ntemplates").remove();
   charts.cov.selectAll(".ntemplates")
     .data(data.ntemplates)
     .enter()
     .append("line")
     .filter(function(d) { return d[0] == data.cov[i_time][0]; })
     .attr("class", "ntemplates")
     .attr("x1", x(0))
     .attr("x2", function(d) { return x(data.cov[i_time][1].length); })
     .attr("y1", function(d) { return y(d[1]); })
     .attr("y2", function(d) { return y(d[1]); })
     .attr("stroke", "darkred")
     .attr("stroke-width", 15)
     .attr("opacity", 0.3);

 }

  // update the genome
  update_genome(data.genome, id, x, height - height_cov, vpad);

}

function update_genome(data, id, x, height, vpad) {

 var tip = d3.tip()
      .attr('class', 'd3-tip')
      .html(function(d) { return d.name + ": " + (+d.location[0][0] + 1) + ", " + d.location[0][1]; });

  var chart = d3.select(".genome-chart"+id);
  var charty = +(chart.attr("transform").split(",")[1].split(")")[0]);

  var chart_main = d3.select("#main-chart"+id);
  
  // add tooltip
  chart.call(tip);
 
  var div_width = $('.svg-container').width();

  // vertical lines at the sides
  chart.selectAll(".side-line")
       .data([0, x.range()[1]])
       .enter()
       .append("line")
       .attr("class", "side-line")
       .attr("x1", function(d) { return d; })
       .attr("x2", function(d) { return d; })
       .attr("y1", -vpad)
       .attr("y2", height - vpad)
       .style("stroke", "black")
       .style("stroke-width", 1);

  var height_block = 20, vpad_block = 5;

  function get_n_feature_types(data) {
   features = [];
   var found = 0;
   for(i = 0; i < data.features.length; i++) {
    found = 0;
    for(j = 0; j < features.length; j++) {
     if (features[j] == data.features[i].type) {
      found = 1;
     }
    }
    if (found == 0) {
     features.push(data.features[i].type);
    }
   }
   return features.length; 
  }

  function get_feature_group(groupname) {
   var group = [];
   var gname = "";
   var tmpfea;
   for(i = 0; i < data.features.length; i++) {
    gname = data.features[i].type
    if (gname == groupname) {
     group.push(data.features[i]);
    } else if (gname + '-single' == groupname) {
     if (['vpu', 'vpr', 'nef', 'vif'].indexOf(data.features[i].name) != -1) {
      group.push(data.features[i]);
     }
    } else if (gname + '-exon1' == groupname) {
     if (['rev', 'tat'].indexOf(data.features[i].name) != -1) {
      tmpfea = {"name": data.features[i].name + ' exon 1',
                "type": "gene",
                "location": [data.features[i].location[0]]};
      group.push(tmpfea);
     }
    } else if (gname + '-exon2' == groupname) {
     if (['rev', 'tat'].indexOf(data.features[i].name) != -1) {
      tmpfea = {"name": data.features[i].name + ' exon 2',
                "type": "gene",
                "location": [data.features[i].location[1]]};
      group.push(tmpfea);
     }
    }
   }
   return group;
  }

  function height_feature(feature) {
   if (feature.type == "fragment") {
    return height_fragment(feature);
   } else if ((feature.type == "gene") | (feature.type == "protein")) {
    return height_inframe(feature); 
   } else {
    return 0;
   }
  }

  function height_fragment(fragment) {
   var fn = +(fragment.name[1]);
   var height = vpad_block;
   if ((fn % 2) == 0) {
    height += height_block + vpad_block;
   }
   return height;
  }

  function height_inframe(feature) {
   var frame;
   if (['tat exon 2', 'rev exon 2'].indexOf(feature.name) != -1) {
     frame = (+(feature.location[0][1]) - (+data.framestart)) % 3;
   } else {
     frame = (+(feature.location[0][0]) - (+data.framestart)) % 3;
   } 
   return frame * (height_block + vpad_block);
  }

  function moverFeature(d) {
    // find position of this box (FIXME: do it better)
    var boxy = +(d3.select(this).attr("transform").split(",")[1].split(")")[0]);

    // change color
    var rect = d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect");

    rect.style("fill", "darkred");  

    // show tooltip
    tip.show(d);

    // show line connecting the exons
    var exons = [];
    if (d.name.indexOf('exon 1') != -1) {
     exons = ["exon-1", "exon-2"];
    } else if (d.name.indexOf('exon 2') != -1) {
     exons = ["exon-2", "exon-1"];
    }
    if (exons.length > 0) {
      var rect2 = d3.select("#" + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "rect");
      var d2 = rect2.data()[0];
      var box2 = d3.select('#' + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "-box");
      // FIXME: rewrite (?): we extract the y coordinate from the affine transformation STRING!
      var box2y = +(box2.attr("transform").split(",")[1].split(")")[0]);

      chart.append('line')
       .attr("class", "exon-line")
       .attr("x1", 0.5 * (x(d.location[0][1]) + x(d.location[0][0])))
       .attr("x2", 0.5 * (x(d2.location[0][1]) + x(d2.location[0][0])))
       .attr("y1", boxy + 0.5 * height_block)
       .attr("y2", box2y + 0.5 * height_block)
       .style("opacity", 0.7)
       .style("stroke", "darkred")
       .style("stroke-width", 2);
    }

   // show vertical "grid" lines at the edges of the region
   chart_main.selectAll()
        .data([d.location[0][0], d.location[0][1]])
	.enter()
	.append("line")
        .attr("class", "grid feature-edge-line")
        .attr("x1", function(dd) { return x(dd); })
        .attr("x2", function(dd) { return x(dd); })
        .attr("y1", 0)
        .attr("y2", charty + boxy)
	.style("fill", "none")
	.style("stroke", "grey")
	.style("stroke-width", 1);

  }

  function moutFeature(d) {
    // change color back
    d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect")
     .style("fill", "steelblue");

    // hide tooltip
    tip.hide(d);   

    // hide exon line
    if ((d.name.indexOf('exon 1') != -1) | (d.name.indexOf('exon 2') != -1)) {
      chart.selectAll(".exon-line").remove();
    }

    // hide vertical "grid" lines
   chart_main.selectAll(".feature-edge-line")
          .remove();
  }

  function plot_group(groupname, dy) {
   var group = get_feature_group(groupname);
   var fea = chart.selectAll()
                  .data(group)
        	   .enter()
        	   .append("g")
                  .attr("class", "featurebox " + groupname)
       	   .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "-box"; })
        	  .attr("transform", function(d) { return "translate(" + x(d.location[0][0]) + "," + (dy + height_feature(d)) + ")";})
        .on('mouseover', moverFeature)
        .on('mouseout', moutFeature);

   var fearect = fea.append("rect")
        .attr("class", "featurerect")
        .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "rect"; })
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", function(d) { return x(d.location[0][1]) - x(d.location[0][0]); })
        .attr("height", 20)
        .style("fill", "steelblue")
	.style("fill-opacity", 0.5)
	.style("stroke-width", 1)
	.style("stroke", "black");

   // show text only of longer things (that's why we have a tooltip)
   fea.append("text")
          .attr("x", function(d) { return 0.5 * (x(d.location[0][1]) - x(d.location[0][0])); })
          .attr("y", 15)
          .text(function(d) {
            if ((d.location[0][1] - d.location[0][0]) > 350) {
              return d.name;
            } else {
              return "";
            }
          })
          .style("text-anchor", "middle");
   }

  plot_group("fragment", 0);
  plot_group("protein", 75);
  plot_group("gene-single", 75);
  plot_group("gene-exon1", 75);
  plot_group("gene-exon2", 75);
  plot_group("RNA_structure", 170);
  plot_group("other", 170);

}
 
