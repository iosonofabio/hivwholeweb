/* 
 * Plot phylogenetic trees with radial or rectangular representation
 * Author: Fabio Zanini
 */
function maxDepth(n, maxOffset) {
  if (n.parent == null) n.depthScaled = n.length;
  else n.depthScaled = n.parent.depthScaled + n.length;

  // If an internal node, maxOffset is given by children
  if (n.children) {
    var i;
    for(i=0; i < n.children.length; i++) {
      var tmp = maxDepth(n.children[i], maxOffset);
      if (tmp > maxOffset) maxOffset = tmp;
    }

  // if terminal node, just checks itself
  } else {
   if (n.depthScaled > maxOffset) {
    maxOffset = n.depthScaled;
   }
  }
  return maxOffset;
}

function phyloScale(n, treeScale) {
  n.y = n.depthScaled * treeScale;
  if (n.children)
    n.children.forEach(function(n) {
      phyloScale(n, treeScale);
    });
}

function getNumberTerminals(n, number) {
    if (n.branchset) {
        var i;
        for(i=0; i < n.branchset.length; i++) {
            number += getNumberTerminals(n.branchset[i], 0);
        }
    } else
        number += 1;
    return number;
}

// NOTE: this function is ready for callbacks, they have to set data = {id: <id>, chartType: <chartType>}
function update(data) {
    var div = d3.select("#phylogram_"+data.id),
        divWidth = +($("#phylogram_"+data.id).width());

    // if this function is called with some useful data, bind it to the DOM
    if ("tree" in data)
        div.datum(data)

    var chart = treeChart().svgWidth(0.9 * divWidth)
                           .chartType(data.chartType);

    if (data.chartType == "radial")
        chart.svgHeight(0.9 * divWidth);
    else
        chart.svgHeight(15 * getNumberTerminals(div.datum().tree, 0));

    div.call(chart);

}

/* tree chart closure as of Mark Bostock: http://bost.ocks.org/mike/chart/ */
function treeChart() {

    var svgWidth = 400,
        svgHeight = 400,
        margin = {top: 5, bottom: 5, left: 5, right: 5},
        width = svgWidth - margin.left - margin.right,
        height = svgHeight - margin.top - margin.bottom,
        chartType = "radial";

    function chart(selection) {
        selection.each(function (data) {
        
            // TREE CHART FUNCTION
            var tree = data.tree,
                id = data.id,
                pname = data.pname,
                region = data.region;
           
            // Create outer chart (SVG) and make sure there are no other ones
            // TODO: recycle charts
            var div = d3.select(this);
            div.selectAll("svg").remove();
            var svg = div.append("svg");

            // Set the outer dimensions.
            //responsive SVG needs these 2 attributes and no width and height attr
            svg.attr("preserveAspectRatio", "xMinYMin meet")
               .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom));
            //svg.attr("width", width + margin.left + margin.right)
            //   .attr("height", height + margin.top + margin.bottom)

            if (chartType == "radial")
                makeRadial();
            else
                makeRectangular();

            // RADIAL CHART
            function makeRadial() {
            
                // Create inner chart, centered in the center of the circle
                var maxAngle = 360,
                    r = width / 2,
                    rInternal = r - 170,
                    rLeaves = rInternal - 100,
                    vis = svg.append("g")
                             .attr("transform", "translate(" + (margin.top + r) + "," + (margin.left + r) + ")");
           
                var cluster = d3.layout.cluster()
                   .size([maxAngle, 1])
                   .sort(null)
                   .value(function(d) { return d.length; })
                   .children(function(d) { return d.branchset; })
                   .separation(function(a, b) { return 1; });
           
                var nodes = cluster.nodes(tree);
                var depth = maxDepth(nodes[0], 0),
                    treeScale = 0.9 * rInternal / depth;
                
                // adjust the bar to calibration
                var barLengthData = (30.0 / treeScale).toPrecision(1),
                    barLength = treeScale * barLengthData;
           
                // add scaled depths to the tree
                phyloScale(nodes[0], treeScale);
           
                var link = vis.selectAll("path.link")
                     .data(cluster.links(nodes))
                     .enter()
                     .append("path")
                     .attr("class", "link")
                     .attr("d", stepRadial)
                     .attr("fill", "none")
                     .attr("stroke", "black")
                     .attr("stroke-width", 2);
           
                // scale bar
                var bar = vis.append("g")
                             .attr("class", "lengthbar")
                          .attr("transform", "translate(" + (r - 50) + "," + (r - 20) + ")");
           
                bar.selectAll(".lengthbar")
                   .data([[-barLength, 0, 0, 0], [-barLength, -barLength, -7, 7], [0, 0, -7, 7]])
                   .enter()
                   .append("svg:line")
                   .attr("x1", function(d) { return d[0]; })
                   .attr("x2", function(d) { return d[1]; })
                   .attr("y1", function(d) { return d[2]; })
                   .attr("y2", function(d) { return d[3]; })
                   .style("stroke", "black")
                   .style("stroke-width", 2);
           
                bar.append("text")
                 .attr("x", -barLength / 2)
                 .attr("y", 25)
                 .text(barLengthData)
                 .attr("text-anchor", "middle");
           
                // line connecting the leaves to their labels
                var label = vis.selectAll(".anno")
                     .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
                     .enter()
                     .append("g")
                     .attr("class", "anno");
           
                 label.append("path")
                     .attr("class", "anno")
                     .attr("d", function(d) { return stepAnnoRadial(d, rInternal); })
                     .attr("fill", "none")
                     .attr("stroke", "lightgrey")
                     .style("stroke-dasharray", ("3, 3"))
                     .attr("stroke-width", 2);
           
                 // leaf labels
                 label.append("text")
                     .attr("dy", ".31em")
                     .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
                     .attr("transform", function(d) {
                         return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 8) + 
                             ")rotate(" + (d.x < 180 ? 0 : 180) + ")"; })
                     .text(function(d) {
                      // local trees have the sequences attached
                      if (region.indexOf('minor') == -1) {
                       return d.name.replace(/_/g, ' ');
                      } else {
                    return d.name.split('_').slice(1).join(" ");
                      }
                     })
                     .on("mouseover", moverRadial)
                     .on("mouseout", moutRadial);
           
                function moverRadial(d) {
                  var t = projectRadial(d);
                  vis.append("circle")
                      .attr("class", "highlight")
                      .attr("cx", t[0])
                      .attr("cy", t[1])
                      .attr("r", 8)
                      .style("stroke", "steelblue")
                      .style("stroke-width", 3)
                      .style("fill", "none");
                
                }
           
                function moutRadial(d) {
                  vis.selectAll(".highlight")
                     .remove();
                }
           
                function projectRadial(d) {
                  var r = d.y, a = (d.x - 90) / 180 * Math.PI;
                  return [r * Math.cos(a), r * Math.sin(a)];
                }
                
                function stepRadial(d) {
                  var s = projectRadial(d.source),
                      m = projectRadial({x: d.target.x, y: d.source.y}),
                      t = projectRadial(d.target),
                      r = d.source.y,
                      sweep = d.target.x > d.source.x ? 1 : 0;
                  return (
                    "M" + s[0] + "," + s[1] +
                    "A" + r + "," + r + " 0 0," + sweep + " " + m[0] + "," + m[1] +
                    "L" + t[0] + "," + t[1]);
                }
                
                function stepAnnoRadial(d, r) {
                  var s = projectRadial({x: d.x, y: d.y + 5}),
                      t = projectRadial({x: d.x, y: r});
                  return (
                    "M" + s[0] + "," + s[1] +
                    "L" + t[0] + "," + t[1]);
                }


            }

            // RECTANGULAR CHART
            function makeRectangular() {

                var vis = svg.append("g")
                             .attr("transform", "translate(" + (margin.top) + "," + (margin.left) + ")");

                // note: at present, x and y are swapped to keep consistency with the radial layout
                var cluster = d3.layout.cluster()
                   .size([height, 0.85 * width])
                   .sort(null)
                   .value(function(d) { return d.length; })
                   .children(function(d) { return d.branchset; })
                   .separation(function(a, b) { return 1; });
           
                var nodes = cluster.nodes(tree);
                var depth = maxDepth(nodes[0], 0),
                    treeScale = cluster.size()[1] / depth;
                
                // adjust the bar to calibration
                var barLengthData = (30.0 / treeScale).toPrecision(1),
                    barLength = treeScale * barLengthData;
           
                // add scaled depths to the tree
                phyloScale(nodes[0], treeScale);
                
                var link = vis.selectAll("path.link")
                     .data(cluster.links(nodes))
                     .enter()
                     .append("path")
                     .attr("class", "link")
                     .attr("d", stepRectangular)
                     .attr("fill", "none")
                     .attr("stroke", "black")
                     .attr("stroke-width", 2);

                // scale bar
                var bar = vis.append("g")
                             .attr("class", "lengthbar")
                             .attr("transform", "translate(" + 50 + "," + (height - 20) + ")");
           
                bar.selectAll(".lengthbar")
                   .data([[-barLength, 0, 0, 0], [-barLength, -barLength, -7, 7], [0, 0, -7, 7]])
                   .enter()
                   .append("svg:line")
                   .attr("x1", function(d) { return d[0]; })
                   .attr("x2", function(d) { return d[1]; })
                   .attr("y1", function(d) { return d[2]; })
                   .attr("y2", function(d) { return d[3]; })
                   .style("stroke", "black")
                   .style("stroke-width", 2);
           
                bar.append("text")
                   .attr("x", -barLength / 2)
                   .attr("y", 25)
                   .text(barLengthData)
                   .attr("text-anchor", "middle");

                // dashed lines going from leaf to label
                var label = vis.selectAll(".anno")
                     .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
                     .enter()
                     .append("g")
                     .attr("class", "anno");
           
                 label.append("path")
                     .attr("class", "anno")
                     .attr("d", function(d) { return stepAnnoRectangular(d); })
                     .attr("fill", "none")
                     .attr("stroke", "lightgrey")
                     .style("stroke-dasharray", ("3, 3"))
                     .attr("stroke-width", 2);

                 // leaf labels
                 label.append("text")
                     .attr("dy", ".31em")
                     .attr("text-anchor", "end")
                     .attr("transform", function(d) { return "translate(" + width + "," + d.x + ")"; })
                     .text(function(d) {
                        // local trees have the sequences attached
                        if (region.indexOf('minor') == -1)
                            return d.name.replace(/_/g, ' ');
                        else
                            return d.name.split('_').slice(1).join(" ");
                     })
                     .on("mouseover", moverRectangular)
                     .on("mouseout", moutRectangular);
           
                function moverRectangular(d) {
                  vis.append("circle")
                      .attr("class", "highlight")
                      .attr("cx", d.y)
                      .attr("cy", d.x)
                      .attr("r", 8)
                      .style("stroke", "steelblue")
                      .style("stroke-width", 3)
                      .style("fill", "none");
                
                }
           
                function moutRectangular(d) {
                  vis.selectAll(".highlight")
                     .remove();
                }


                function stepRectangular(d) {
                    return (
                        "M" + d.source.y + "," + d.source.x +
                        "V" + d.target.x +
                        "H" + d.target.y);
                }

                function stepAnnoRectangular(d) {
                  return (
                    "M" + d.y + "," + d.x +
                    "H" + (10 + depth * treeScale));
                }

            }

        });
    }

    chart.margin = function (_) {
        if (!arguments.length) return margin;
        margin = _;
        width = svgWidth - margin.left - margin.right;
        height = svgHeight - margin.top - margin.bottom;
        return chart;
    };

    chart.svgWidth = function (_) {
        if (!arguments.length) return width;
        svgWidth = _;
        width = svgWidth - margin.left - margin.right;
        return chart;
    };

    chart.svgHeight = function (_) {
        if (!arguments.length) return height;
        svgHeight = _;
        height = svgHeight - margin.top - margin.bottom;
        return chart;
    };

    chart.chartType = function (_) {
        if (!arguments.length) return chartType;
        chartType = _;
        return chart;
    };

    return chart;
}
