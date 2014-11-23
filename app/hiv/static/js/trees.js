function cross(a, b) { return a[0] * b[1] - a[1] * b[0]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1]; }

function project(d) {
  var r = d.y, a = (d.x - 90) / 180 * Math.PI;
  return [r * Math.cos(a), r * Math.sin(a)];
}

function maxDepth(n, offset, maxOffset) {
  var tmp;
  if (n.length != null) offset += n.length;
  if (offset > maxOffset) maxOffset = offset;
  if (n.children) {
    for(i=0; i < n.children.length; i++) {
      tmp = maxDepth(n.children[i], offset, maxOffset);
      if (tmp > maxOffset) maxOffset = tmp;
    }
  }
  return maxOffset;
}

function phylo(n, offset, treeScale) {
  if (n.length != null) offset += treeScale * n.length;
  n.y = offset;
  if (n.children)
    n.children.forEach(function(n) {
      phylo(n, offset, treeScale);
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

function update(text, id, pname) {

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
 var depth = maxDepth(nodes[0], 0, 0),
     treeScale = 100 / depth;
 phylo(nodes[0], 0, treeScale);

 var link = vis.selectAll("path.link")
      .data(cluster.links(nodes))
      .enter()
      .append("path")
      .attr("class", "link")
      .attr("d", step)
      .attr("fill", "none")
      .attr("stroke", "black")
      .attr("stroke-width", 2);

 /* label nodes with circles
 var node = vis.selectAll("g.node")
      .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })

  node.append("circle")
      .attr("r", 5)
      .style("stroke", "steelblue")
      .style("stroke-width", 2)
      .style("fill", "none");
  */

 var annoLine = vis.selectAll("path.anno")
      .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
      .enter()
      .append("path")
      .attr("class", "anno")
      .attr("d", function(d) { return stepAnno(d, rInternal); })
      .attr("fill", "none")
      .attr("stroke", "lightgrey")
      .style("stroke-dasharray", ("3, 3"))
      .attr("stroke-width", 2);

 var label = vis.selectAll("text")
      .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
    .enter().append("text")
      .attr("dy", ".31em")
      .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 8) + ")rotate(" + (d.x < 180 ? 0 : 180) + ")"; })
      .text(function(d) { return d.name.replace(/_/g, ' '); });
}
