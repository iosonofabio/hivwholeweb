function update(data, id) {

 var margin = {'top': 20, 'bottom': 10, 'left': 20, 'right': 20},
     width = $('#haploThumbnail').width() - margin.left - margin.right,
     height = $('#haploThumbnail').height() - margin.top - margin.bottom;

  var vis = d3.select("#haploThumbnail").append("g")
       .attr("id", "haploThumbnailVis")
       .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
       .datum(data);

  var x = d3.scale.linear()
           .domain([0, data.len + 50])
	   .range([0, width])
 
  vis.append("line")
   .attr("x1", x.range()[0])
   .attr("x2", x.range()[1])
   .attr("y1", 0)
   .attr("y2", 0)
   .style("stroke", "black")
   .style("stroke-width", 5);

  addRegion($("#region option:selected").text());

}


function addRegion(name) {
 var margin = {'top': 10, 'bottom': 10, 'left': 20, 'right': 20},
     width = $('#haploThumbnail').width() - margin.left - margin.right,
     height = $('#haploThumbnail').height() - margin.top - margin.bottom;

 var vis = d3.select("#haploThumbnailVis");
 var data = vis.datum();

 var x = d3.scale.linear()
          .domain([0, data.len + 50])
          .range([0, width])

 vis.selectAll(".regionGroup")
  .transition()
  .duration(300)
  .style("opacity", 0)
  .remove();

 var regionGroup = vis.selectAll(".regionGroup")
  .data(data.features)
  .enter()
  .append("g")
  .filter(function(d) { return d.name == name; })
  .attr("class", "regionGroup")
  .attr("transform", function(d) { return "translate(" + x(d.location[0][0]) + "," + 20 + ")"; })
  .style("opacity", 0);

 regionGroup.append("rect")
  .attr("x", 0)
  .attr("y", -20)
  .attr("width", function(d) { return x(d.location[0][1]) - x(d.location[0][0]); })
  .attr("height", 20)
  .style("fill", "steelblue")
  .style("opacity", 0.5)
  .style("stroke", "none");

 regionGroup.append("text")
  .text(function(d) { return "HXB2: " + d.location[0][0] + " - " + d.location[0][1]; })
  .attr("x", function(d) { return x(0.5 * (d.location[0][1] - d.location[0][0])); })
  .attr("y", 20)
  .style("text-anchor", "middle");

 regionGroup.transition()
  .duration(300)
  .style("opacity", 1);
  
}
