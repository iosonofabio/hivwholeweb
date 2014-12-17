function update(data, id) {

 var datal = data.len;

 var div_width = $('.svg-container').width();

 var margin = {top: 45, right: 50, bottom: 50, left: 80},
     width = div_width - margin.left - margin.right,
     height = $('.'+id).height() - margin.top - margin.bottom,
     vpad = 50, height_genome=200,
     height_aft = (height - height_genome - vpad - vpad) / 2;

 var chart_ext = d3.select("."+id)
     .attr("width", width + margin.left + margin.right)
     .attr("height", height + margin.bottom + margin.top);

 var vis = chart_ext.append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var visAft = vis.append("g")
      .attr("class", "d3-chart aft-chart"+id);

 var visAfg = vis.append("g")
      .attr("class", "d3-chart afg-chart"+id)
      .attr("transform", "translate(0," + (height_aft + vpad) + ")");

 var visGenome = vis.append("g")
          .attr("class", "d3-chart genome-chart"+id)
          .attr("transform", "translate(0," + (2 * height_aft + 2 * vpad) + ")");

 var y = d3.scale.log()
      .domain([0.0009, 0.9991])
      .range([height_aft, 0])
      .clamp(true);
 
 var x = d3.scale.linear()
      .domain([-10, 1.05 * data.tmax])
      .range([0, width])
      .clamp(true);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");
 
 var yAxis = d3.svg.axis()
     .scale(y)
     .orient("left");

 var yAxisRight = d3.svg.axis()
     .scale(y)
     .orient("right");

 visAft.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(0," + height_aft + ")")
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
      .attr("x", -(height_aft) / 2)
      .style("text-anchor", "middle")
      .text("Frequency");

 visAft.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(" + width + ",0)")
      .call(yAxisRight);


 var lineFunc = d3.svg.line()
                      .x(function(d) { return x(d[0]); })
       	       .y(function(d) { return y(d[1]); })
       	       .interpolate('monotone');

 var colors = d3.scale.linear()
 .domain([0, data.len / 4, data.len / 3, data.len / 2, 2 * data.len / 3, 3 * data.len / 4, data.len])
 .interpolate(d3.interpolateRgb)
 .range(["darkblue", "blue", "cyan", "green", "yellow", "orange", "red"]);

 visAft.selectAll(".aft")
      .data(data.aft)
      .enter()
      .append("svg:path")
      .attr("class", "aft")
      .attr("id", function(d, i) { return "aft-" + (i+1); })
      .attr("d", function(d) { return lineFunc(d[2]); })
      .attr("stroke", function(d) { return colors(d[0]); })
      .attr("stroke-width", 2)
      .attr("fill", "none")
      .attr("opacity", 0.5);

 // Add number of templates line if present
 visAft.append("svg:path")
   .attr("class", "ntemplates")
   .attr("d", d3.svg.line()
                .x(function(d) { return x(d[0]); })
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
      .attr("x1", x(0))
      .attr("x2", x.range()[1])
      .attr("y1", y(2e-3))
      .attr("y2", y(2e-3))
      .attr("stroke", "steelblue")
      .attr("stroke-width", 15)
      .attr("opacity", 0.4);


  // update the allele frequencies along the genome
  updateAfg(data.aft, id, width, height_aft, vpad, data.genome.len, colors, x, data.times, data.ntemplates);

  // update the genome
  update_genome(data.genome, id, width, height_genome, vpad, height_aft);

}

function updateAfg (data, id, width, height, vpad, len, colors, xsl, times, ntemplates) {
 var x = d3.scale.linear()
      .domain([-0.05 * len, 1.05 * len])
      .range([0, width]);

 var y = d3.scale.log()
      .domain([0.0009, 0.9991])
      .range([height, 0])
      .clamp(true);

  var xt = d3.scale.linear()
      .domain([times[0], times[times.length - 1]])
      .range([xsl(times[0]), xsl(times[times.length - 1])])
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

 var vis = d3.select(".afg-chart"+id);

 // Slider
 vis.append("g")
     .attr("class", "x axis")
     .attr("transform", "translate(0," + (-height -vpad - 20) + ")")
     .call(d3.svg.axis()
       .scale(xt)
       .orient("top"))
   .select(".domain")
   .select(function() { return this.parentNode.appendChild(this.cloneNode(true)); })
     .attr("class", "halo");

 var brush = d3.svg.brush()
     .x(xt)
     .extent([1, 1])
     .on("brush", brushed);

  var slider = vis.append("g")
     .attr("class", "slider")
     .call(brush);
 
 slider.selectAll(".extent,.resize")
     .remove();

 var handle = slider.append("circle")
     .attr("class", "handle")
     .attr("transform", "translate(0," + (-vpad - height - 20) + ")")
     .attr("r", 5);

 slider.call(brush.event)
     .call(brush.extent([1, 1]))
     .call(brush.event);

 // NOTE: we should make transitions here, but this involves enumerating the cases and is therefore very tedious
 function brushed() {
   var value = brush.extent()[0];
 
   if (d3.event.sourceEvent) { // not a programmatic event
     value = xt.invert(d3.mouse(this)[0]);
     brush.extent([value, value]);
   }

   var time = getTime(times, value);
   handle.attr("cx", xt(time));

   // update time tracker
   vis.select("#aftTimeTracker")
    .attr("x1", xt(time))
    .attr("x2", xt(time));

   // update data
   vis.selectAll(".afg_point")
       .remove();

   vis.select(".afg")
        .selectAll(".afg_point")
        .data(data)
        .enter()
        .append("circle")
        .filter(function(d) { return hasTime(d[2], time); })
        .attr("class", "circle afg_point")
        .attr("cx", function(d) { return x(d[0]); })
        .attr("cy", function(d) { return y(getFreq(d[2], time)); })
        .attr("r", 3)
        .attr("fill", function(d) { return colors(d[0]); })
        .attr("stroke", "grey")
        .attr("stroke-width", 0.3);

  // update ntemplates
  vis.select("#afgTemplates")
   .remove();
  
  if (hasTime(ntemplates, time)) {
   vis.append("svg:line")
      .attr("id", "afgTemplates")
      .datum(1.0 / getFreq(ntemplates, time))
      .attr("x1", xsl.range()[0])
      .attr("x2", xsl.range()[1])
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
      .attr("x", -(height) / 2)
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
      .data(data)
      .enter()
      .append("circle")
      .filter(function(d) { return hasTime(d[2], times[0]); })
      .attr("class", "circle afg_point")
      .attr("cx", function(d) { return x(d[0]); })
      .attr("cy", function(d) { return y(d[2][0][1]); })
      .attr("r", 3)
      .attr("fill", function(d) { return colors(d[0]); })
      .attr("stroke", "grey")
      .attr("stroke-width", 0.3);

 // number of templates if appropriate
 if (hasTime(ntemplates, times[0])) {
  vis.append("svg:line")
     .attr("id", "afgTemplates")
     .datum(1.0 / getFreq(ntemplates, times[0]))
     .attr("x1", xsl.range()[0])
     .attr("x2", xsl.range()[1])
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
     .attr("x1", xsl.range()[0])
     .attr("x2", xsl.range()[1])
     .attr("y1", function(d) { return y(d); })
     .attr("y2", function(d) { return y(d); })
     .attr("stroke", "steelblue")
     .attr("stroke-width", 15)
     .attr("fill", "none")
     .attr("opacity", 0.4);


 // vertical line for time tracking
 vis.append("svg:line")
  .attr("id", "aftTimeTracker")
  .attr("x1", xt(0))
  .attr("x2", xt(0))
  .attr("y1", -height - vpad - 15)
  .attr("y2", -vpad)
  .attr("stroke", "grey")
  .attr("stroke-width", 5)
  .attr("opacity", 0.5);

}

function update_genome(data, id, width, height, vpad, height_aft) {

 var tip = d3.tip()
      .attr('class', 'd3-tip')
      .html(function(d) { return d.name + ": " + (+d.location[0][0] + 1) + ", " + d.location[0][1]; });

 var chart = d3.select(".genome-chart"+id);
 var charty = +(chart.attr("transform").split(",")[1].split(")")[0]);

 // add tooltip
 chart.call(tip);
 
 var x = d3.scale.linear()
      .domain([-0.05 * data.len, 1.05 * data.len])
      .range([0, width]);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");

 chart.append("g")
      .attr("class", "d3-axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 40)
      .style("text-anchor", "middle")
      .text("Position [bp]");
 
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

   // thicken lines in this feature
   d3.selectAll(".aft")
     .attr("opacity", function(daft) {
       if ((daft[0] >= d.location[0][0]) & (daft[0] < d.location[0][1])) {
	return 1;
       } else {
        return 0.05;
       }
     });

   // show vertical "grid" lines at the edges of the region
   chart.selectAll(".feature-edge-line")
        .data([d.location[0][0], d.location[0][1]])
	.enter()
	.append("line")
        .attr("class", "grid feature-edge-line")
        .attr("x1", function(dd) { return x(dd); })
        .attr("x2", function(dd) { return x(dd); })
        .attr("y1", -(vpad + height_aft))
        .attr("y2", boxy)
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

   // restore opacity of lines
   d3.selectAll(".aft")
     .attr("opacity", 0.5);   

    // hide vertical "grid" lines
   chart.selectAll(".feature-edge-line")
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
 
function getTime(times, value) {
  var t = times[0], d = Math.abs(value - times[0]), dnew;
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
 for (i = 0; i < arr.length; i++) {
  if (arr[i][0] == time) return true;
 }
 return false;
}

function getFreq(arr, time) {
 for (i = 0; i < arr.length; i++) {
  if (arr[i][0] == time) return arr[i][1];
 }
 return undefined;
}
