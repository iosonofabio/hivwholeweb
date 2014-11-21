function update(treestr, id, pname) {
  var newick = Newick.parse(treestr)

 
  var newickNodes = []
  function buildNewickNodes(node, callback) {
    newickNodes.push(node)
    if (node.branchset) {
      for (var i=0; i < node.branchset.length; i++) {
        buildNewickNodes(node.branchset[i])
      }
    }
  }
  buildNewickNodes(newick)

  if (pname == "all") {
    var width = 600
    var height = 1000	
  } else {
    var width = 400
    var height = 300
  }

  var dict = d3.phylogram.build("#phylogram_"+id, newick, {
    width: width,
    height: height
  });

  dict.vis.selectAll("circle")
   .attr('stroke',  'steelblue')
   .attr('fill', 'steelblue');

}
