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

 function update(data, id) {

  var tip = d3.tip()
       .attr('class', 'd3-tip')
       .html(function(d) { return d.name + ": " + (+d.location[0][0] + 1) + ", " + d.location[0][1]; });

  var n_feature_types = get_n_feature_types(data);

  var div_width = $('.svg-container').width();

  var margin = {top: 3, right: 60, bottom: 42, left: 60},
      width = div_width - margin.left - margin.right,
      height = 200, height_block = 20, vpad_block = 5;

   var x = d3.scale.linear()
        .domain([-50, data.len + 50])
        .range([20, width - 20]);

   var xAxis = d3.svg.axis()
       .scale(x)
       .orient("bottom");

   // external chart, including axis
   var chart_ext = d3.select("."+id)
       .attr("width", width + margin.left + margin.right)
       .attr("height", height + margin.bottom + margin.top);

   // internal chart, only the bars
   var chart = chart_ext.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

   // add tooltip
   chart.call(tip);

   // add axis to the chart
   var xAxisObj = chart.append("g")
                       .attr("class", "d3-axis")
                       .attr("transform", "translate(0," + height + ")")
                       .call(xAxis)

   xAxisObj.append("text")
           .attr("x", width / 2)
           .attr("y", 40)
           .style("text-anchor", "middle")
	   .text("Position [bp]");

   // add the rest of the rectangle box
   xAxisObj.append("line")
    .attr({"class": "bBox", "x1": x(-50), "x2": x(-50), "y1": 0, "y2": -height});
   xAxisObj.append("line")
    .attr({"class": "bBox", "x1": x(data.len + 50), "x2": x(data.len + 50), "y1": 0, "y2": -height});
   xAxisObj.append("line")
    .attr({"class": "bBox", "x1": x(-50), "x2": x(data.len + 50), "y1": -height, "y2": -height});

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
    } else if ((feature.type == "protein") | (feature.type == "gene")) {
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
       var boxy = +(d3.select(this).attr("transform").split(",")[1].split(")")[0]);
       var box2y = +(box2.attr("transform").split(",")[1].split(")")[0]);

       chart.append('line')
	.attr("class", "exon-line")
	.attr("x1", 0.5 * (x(d.location[0][1]) + x(d.location[0][0])))
	.attr("x2", 0.5 * (x(d2.location[0][1]) + x(d2.location[0][0])))
	.attr("y1", boxy + 0.5 * height_block)
	.attr("y2", box2y + 0.5 * height_block)
        .style("stroke", "darkred")
        .style("stroke-width", 2);
     }
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
         .attr("style", "fill:steelblue;fill-opacity:0.5;stroke-width:1;stroke:rgb(0,0,0)");

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
