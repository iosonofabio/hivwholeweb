function cross(a, b) { return a[0] * b[1] - a[1] * b[0]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1]; }

function project(d) {
  var r = d.y, a = (d.x - 90) / 180 * Math.PI;
  return [r * Math.cos(a), r * Math.sin(a)];
}

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

function step(d) {
  var s = project(d.source),
      m = project({x: d.target.x, y: d.source.y}),
      t = project(d.target),
      r = d.source.y,
      sweep = d.target.x > d.source.x ? 1 : 0;
  return (
    "M" + s[0] + "," + s[1] +
    "A" + r + "," + r + " 0 0," + sweep + " " + m[0] + "," + m[1] +
    "L" + t[0] + "," + t[1]);
}

function stepAnno(d, r) {
  var s = project({x: d.x, y: d.y + 5}),
      t = project({x: d.x, y: r});
  return (
    "M" + s[0] + "," + s[1] +
    "L" + t[0] + "," + t[1]);
}

function update(text, id, pname, region) {

 var maxAngle = 360;
 var div_width = +($("#phylogram_"+id).width()),
     svg_width = 0.9 * div_width,
     margins = {"top": 5, "bottom": 5, "left": 5, "right": 5},
     width = svg_width - margins.left - margins.right,
     height = svg_width - margins.top - margins.bottom,
     r = width / 2,
     rInternal = r - 170,
     rLeaves = rInternal - 100;

 var chart_ext = d3.select("#phylogram_"+id).append("svg")
    .attr("width", width + margins.left + margins.right)
    .attr("height", height + margins.top + margins.bottom);
 
 var vis = chart_ext.append("g")
    .attr("transform", "translate(" + (margins.top + r) + "," + (margins.left + r) + ")");

 var cluster = d3.layout.cluster()
    .size([maxAngle, 1])
    .sort(null)
    .value(function(d) { return d.length; })
    .children(function(d) { return d.branchset; })
    .separation(function(a, b) { return 1; });

 var tree = Newick.parse(text);
 var nodes = cluster.nodes(tree);
 var depth = maxDepth(nodes[0], 0),
     treeScale = 0.9 * rInternal / depth,
     barLength = 30.0,
     barLengthData = (30.0 / treeScale).toPrecision(1);
 
 // adjust the bar to calibration
 barLength = treeScale * barLengthData;

 // add scaled depths to the tree
 phyloScale(nodes[0], treeScale);

 var link = vis.selectAll("path.link")
      .data(cluster.links(nodes))
      .enter()
      .append("path")
      .attr("class", "link")
      .attr("d", step)
      .attr("fill", "none")
      .attr("stroke", "black")
      .attr("stroke-width", 2);

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

 var label = vis.selectAll(".anno")
      .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
      .enter()
      .append("g")
      .attr("class", "anno");

  label.append("path")
      .attr("class", "anno")
      .attr("d", function(d) { return stepAnno(d, rInternal); })
      .attr("fill", "none")
      .attr("stroke", "lightgrey")
      .style("stroke-dasharray", ("3, 3"))
      .attr("stroke-width", 2);

  // leaf labels
  label.append("text")
      .attr("dy", ".31em")
      .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 8) + ")rotate(" + (d.x < 180 ? 0 : 180) + ")"; })
      .text(function(d) {
       // local trees have the sequences attached
       if (region.indexOf('minor') == -1) {
        return d.name.replace(/_/g, ' ');
       } else {
	return d.name.split('_').slice(1).join(" ");
       }
      })
      .on("mouseover", mover)
      .on("mouseout", mout);

function mover(d) {
  var t = project(d);
  vis.append("circle")
      .attr("class", "highlight")
      .attr("cx", t[0])
      .attr("cy", t[1])
      .attr("r", 8)
      .style("stroke", "steelblue")
      .style("stroke-width", 3)
      .style("fill", "none");

}

function mout(d) {
  vis.selectAll(".highlight")
     .remove();
}

}
