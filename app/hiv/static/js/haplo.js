function updateHaplo(id, data) {
    var margin = {'top': 20, 'bottom': 10, 'left': 40, 'right': 40},
        width = $('#'+id).width() - margin.left - margin.right,
        height = $('#'+id).height() - margin.top - margin.bottom;
  
    var svg = d3.select("#"+id),
        vis = svg.append("g")
         .attr("id", "haploThumbnailVis")
         .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
         .datum(data);
  
    var x = d3.scale.linear()
        .domain([1, data.len + 50])
  	.range([0, width]);
   
    vis.append("line")
        .attr("x1", x.range()[0])
        .attr("x2", x.range()[1])
        .attr("y1", 0)
        .attr("y2", 0)
        .style("stroke", "black")
        .style("stroke-width", 5);
  
    addRegion($("#region option:selected").text());
}


function moveRegion(name) {

 var margin = {'top': 10, 'bottom': 10, 'left': 40, 'right': 40},
     width = $('#haploThumbnail').width() - margin.left - margin.right,
     height = $('#haploThumbnail').height() - margin.top - margin.bottom;

 var vis = d3.select("#haploThumbnailVis");
 var data = vis.datum();

 var x = d3.scale.linear()
          .domain([0, data.len + 50])
          .range([0, width]);

 var datum = getFeature(data, name);
 var regionGroup = vis.select("g.regionGroup");

 // NOTE: these transitions should take place at the same time, but that's laborious to code
 regionGroup.select("rect")
  .transition()
  .duration(100)
  .attr("width", x(datum.location[0][1]) - x(datum.location[0][0] + 1));

 regionGroup.select("text")
  .text("HXB2: " + (datum.location[0][0] + 1) + " - " + datum.location[0][1])
  .attr("x", x(0.5 * (datum.location[0][1] - datum.location[0][0])));

 regionGroup.transition()
  .duration(900)
  .attr("transform", "translate(" + x(datum.location[0][0] + 1) + "," + 20 + ")");

}


function addRegion(name) {
 var margin = {'top': 10, 'bottom': 10, 'left': 40, 'right': 40},
     width = $('#haploThumbnail').width() - margin.left - margin.right,
     height = $('#haploThumbnail').height() - margin.top - margin.bottom;

 var vis = d3.select("#haploThumbnailVis");
 var data = vis.datum();

 var x = d3.scale.linear()
          .domain([0, data.len + 50])
          .range([0, width]);

 // NOTE: we get regions in 0 - (max+1) coordinates a la Python, we convert them 1+ both included.
 var datum = getFeature(data, name);
 var regionGroup = vis.append("g")
  .attr("class", "regionGroup")
  .attr("transform", "translate(" + x(datum.location[0][0] + 1) + "," + 20 + ")")
  .style("opacity", 0);

 regionGroup.append("rect")
  .attr("x", 0)
  .attr("y", -20)
  .attr("width", x(datum.location[0][1]) - x(datum.location[0][0] + 1))
  .attr("height", 20)
  .style("fill", "steelblue")
  .style("opacity", 0.5)
  .style("stroke", "none");

 regionGroup.append("text")
  .text("HXB2: " + (datum.location[0][0] + 1) + " - " + datum.location[0][1])
  .attr("x", x(0.5 * (datum.location[0][1] - datum.location[0][0])))
  .attr("y", 20)
  .style("text-anchor", "middle");

 regionGroup.transition()
  .duration(300)
  .style("opacity", 1);
  
}

function getFeature(data, name) {
 var i;
 for(i=0; i < data.features.length; i++) {
  if(data.features[i].name == name) {
   return data.features[i];
  }
 }
}

