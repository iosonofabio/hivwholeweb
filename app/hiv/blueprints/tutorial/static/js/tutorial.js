// Full list of configuration options available here:
// https://github.com/hakimel/reveal.js#configuration
Reveal.initialize({
	controls: true,
	progress: true,
	history: true,
	center: true,

	theme: Reveal.getQueryHash().theme, // available themes are in /css/theme
	transition: Reveal.getQueryHash().transition || 'fade', // default/cube/page/concave/zoom/linear/fade/none

	// Parallax scrolling
	// parallaxBackgroundImage: 'https://s3.amazonaws.com/hakim-static/reveal-js/reveal-parallax-1.jpg',
	// parallaxBackgroundSize: '2100px 900px',

	// Optional libraries used to extend on reveal.js
	dependencies: [
		{ src: 'lib/js/classList.js', condition: function() { return !document.body.classList; } },
		{ src: 'plugin/markdown/marked.js', condition: function() { return !!document.querySelector( '[data-markdown]' ); } },
		{ src: 'plugin/markdown/markdown.js', condition: function() { return !!document.querySelector( '[data-markdown]' ); } },
					{ src: 'plugin/highlight/highlight.js', async: true, callback: function() { hljs.initHighlightingOnLoad(); } },
					{ src: 'plugin/zoom-js/zoom.js', async: true, condition: function() { return !!document.body.classList; } },
					{ src: 'plugin/notes/notes.js', async: true, condition: function() { return !!document.body.classList; } }
				]
});

// d3 part below

function createPhysio1() {

 var secWidth = $("section").width();

 var margin = {left: 200, right:150, top:50, bottom:80},
    width = secWidth - margin.left - margin.right,
    height = 400;

 var vis = d3.select("#d3-plot-time1")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var x = d3.scale.linear()
      .domain([0, 3000])
      .range([0, width]);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");

 vis.append("g")
      .attr("class", "d3-axis-white")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 65)
      .style("text-anchor", "middle")
      .text("Time from infection [days]");


vis.append("defs").append("marker")
    .attr("id", "arrowhead")
    .attr("refX", 1) /*must be smarter way to calculate shift*/
    .attr("refY", 2)
    .attr("markerWidth", 6)
    .attr("markerHeight", 4)
    .attr("orient", "auto")
    .append("path")
        .attr("d", "M 0,0 V4 L6,2 Z")
	.attr("fill", "white");


// FIXME: here an animation would be good, but it's buggy
vis.append("line")
   .attr("marker-end", "url(#arrowhead)")
   .attr("x1", 0)
   .attr("x2", x(2000))
   .attr("y1", height / 2)
   .attr("y2", height / 2)
   .attr("stroke-width", 20)
   .attr("stroke", "white");

}

function createPhysio2() {

 var secWidth = $("section").width();

 var margin = {left: 200, right:150, top:50, bottom:80},
    width = secWidth - margin.left - margin.right,
    height = 400;

 var vis = d3.select("#d3-plot-time2")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var x = d3.scale.linear()
      .domain([0, 3000])
      .range([0, width]);

 var y = d3.scale.log()
      .domain([50, 100000])
      .range([height, 0]);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");

 var yAxis = d3.svg.axis()
     .scale(y)
     .orient("left");

 vis.append("g")
      .attr("class", "d3-axis-white")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 65)
      .style("text-anchor", "middle")
      .text("Time from infection [days]");

 vis.append("g")
      .attr("class", "d3-axis-white")
      .call(yAxis)
      .append("g")
      .attr("class", "d3-axis-textbox")
      .attr("transform", "rotate(-90)")
      .append("text")
      .attr("dy", "-4.8em")
      .attr("x", -height / 2)
      .style("text-anchor", "middle")
      .text("Viral load [copies/ml]");


vis.append("defs").append("marker")
    .attr("id", "arrowhead_red")
    .attr("refX", 1) /*must be smarter way to calculate shift*/
    .attr("refY", 2)
    .attr("markerWidth", 6)
    .attr("markerHeight", 4)
    .attr("orient", "auto")
    .append("path")
        .attr("d", "M 0,0 V4 L6,2 Z")
	.attr("fill", "darkred");


// FIXME: here an animation would be good, but it's buggy
var data = [[150, 1500], [300, 1200], [500, 1500],
            [800, 1800], [1100, 1100], [1500, 1750],
            [1900, 2500], [2100, 5000], [2500, 20000]];
vis.append("path")
   .attr("marker-end", "url(#arrowhead_red)")
   .attr("d", d3.svg.line()
     .x(function(d) { return x(d[0]); })
     .y(function(d) { return y(d[1]); })
     .interpolate("monotone")(data))
   .attr("stroke-width", 5)
   .attr("stroke", "darkred")
   .attr("fill", "none");

}

function createPhysio3() {

 var secWidth = $("section").width();

 var margin = {left: 200, right:150, top:50, bottom:80},
    width = secWidth - margin.left - margin.right,
    height = 400;

 var vis = d3.selectAll("#d3-plot-time3")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var x = d3.scale.linear()
      .domain([0, 3000])
      .range([0, width]);

 var y = d3.scale.log()
      .domain([50, 100000])
      .range([height, 0]);

 var y_cc = d3.scale.linear()
      .domain([0, 900])
      .range([height, 0]);

 var xAxis = d3.svg.axis()
     .scale(x)
     .orient("bottom");

 var yAxis = d3.svg.axis()
     .scale(y)
     .orient("left");

 var yAxis_cc = d3.svg.axis()
     .scale(y_cc)
     .orient("right");

 vis.append("g")
      .attr("class", "d3-axis-white")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
      .attr("x", width / 2)
      .attr("y", 65)
      .style("text-anchor", "middle")
      .text("Time from infection [days]");

 vis.append("g")
      .attr("class", "d3-axis-white")
      .call(yAxis)
      .append("g")
      .attr("class", "d3-axis-textbox")
      .attr("transform", "rotate(-90)")
      .append("text")
      .attr("dy", "-4.8em")
      .attr("x", -height / 2)
      .style("text-anchor", "middle")
      .text("Viral load [copies/ml]");

 var ccTextBox = vis.append("g")
      .attr("class", "d3-axis-white")
      .attr("transform", "translate(" + width + " ,0)")
      .call(yAxis_cc)
      .append("g")
      .attr("class", "d3-axis-textbox")
      .attr("transform", "rotate(+90)");

 ccTextBox.append("text")
      .attr("dy", "-4em")
      .attr("x", height / 2)
      .style("text-anchor", "middle")
      .text("CD4+ cells [cells/ml]");

 ccTextBox.append("rect")
      .attr("class", "circles ccAxis")
      .attr("y", "-4.8em")
      .attr("x", 0.2 * height)
      .attr("width", 30)
      .attr("height", 14)
      .style("fill", "steelblue");


vis.append("defs").append("marker")
    .attr("id", "arrowhead_red")
    .attr("refX", 1) /*must be smarter way to calculate shift*/
    .attr("refY", 2)
    .attr("markerWidth", 6)
    .attr("markerHeight", 4)
    .attr("orient", "auto")
    .append("path")
        .attr("d", "M 0,0 V4 L6,2 Z")
	.attr("fill", "darkred");

vis.append("defs").append("marker")
    .attr("id", "arrowhead_blue")
    .attr("refX", 1) /*must be smarter way to calculate shift*/
    .attr("refY", 2)
    .attr("markerWidth", 6)
    .attr("markerHeight", 4)
    .attr("orient", "auto")
    .append("path")
        .attr("d", "M 0,0 V4 L6,2 Z")
	.attr("fill", "steelblue");


// FIXME: here an animation would be good, but it's buggy
var data = [[150, 1500], [300, 1200], [500, 1500],
            [800, 1800], [1100, 1100], [1500, 1750],
            [1900, 2500], [2100, 5000], [2500, 20000]];
vis.append("path")
   .attr("marker-end", "url(#arrowhead_red)")
   .attr("d", d3.svg.line()
     .x(function(d) { return x(d[0]); })
     .y(function(d) { return y(d[1]); })
     .interpolate("monotone")(data))
   .attr("stroke-width", 5)
   .attr("stroke", "darkred")
   .attr("fill", "none");

var data_cc = [[150, 800], [400, 750], [800, 700],
               [1100, 600], [1500, 650],
               [2100, 500], [2500, 200]];
vis.append("path")
   .attr("marker-end", "url(#arrowhead_blue)")
   .attr("d", d3.svg.line()
     .x(function(d) { return x(d[0]); })
     .y(function(d) { return y_cc(d[1]); })
     .interpolate("monotone")(data_cc))
   .attr("stroke-width", 5)
   .attr("stroke", "steelblue")
   .attr("fill", "none");


}

function createVirion() {

 var secWidth = $("section").width();

 var margin = {left: 30, right:30, top:30, bottom:30},
     width = secWidth - margin.left - margin.right,
     height = 340,
     radius = 80,
     radiusCapsid = radius - 20;

 var vis = d3.selectAll("#virion-zoom")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var vir = vis.append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + (width / 4) + "," + (height / 2) + ")");


 // draw virus circle
 vir.append("circle")
   .attr("cx", 0)
   .attr("cy", 0)
   .attr("r", radius)
   .style("fill", "steelblue")
   .style("fill-opacity", 0.3)
   .style("stroke", "steelblue")
   .style("stroke-width", 8);
 
 
 // draw capsid
 var hexagonData = [
       { "x": radiusCapsid/7,  "y": -radiusCapsid},
       { "x": radiusCapsid/3.5,  "y": -0.9 * radiusCapsid},
       { "x": radiusCapsid/2,   "y": 0.8 * radiusCapsid}, 
       { "x": radiusCapsid/3,   "y": radiusCapsid}, 
       { "x": -radiusCapsid/3,  "y": radiusCapsid},
       { "x": -radiusCapsid/2,  "y": 0.8 * radiusCapsid},
       { "x": -radiusCapsid/3.5,  "y": -0.9 * radiusCapsid},
       { "x": -radiusCapsid/7,  "y": -radiusCapsid},
     ];
 
 drawHexagon = d3.svg.line()
       .x(function(d) { return d.x; })
       .y(function(d) { return d.y; })
       .interpolate("linear-closed");
 
 vir.append("path")
    .attr("transform", "rotate(100)")
    .attr("d", drawHexagon(hexagonData))
    .attr("stroke", "darkred")
    .attr("stroke-width", 5)
    .style("fill", "darkred")
    .style("fill-opacity", 0.7);

  // HIV membrane protein
  var angleGp120 = [0, 25, 50, 75, 100, 125, 150, 175, 200, 226, 253,
                    281, 308, 334];
  var dataGp120 = [
    { "x": -4, "y": 6} ,
    { "x": -4, "y": -8} ,
    { "x": -8, "y": -12} ,
    { "x": 0, "y": -16} ,
    { "x": 8, "y": -12} ,
    { "x": 4, "y": -8} ,
    { "x": 4, "y": 6} ,
  ];
  vir.selectAll(".spike")
     .data(angleGp120)
     .enter()
     .append("path")
     .attr("transform", function(d) {
       return "translate(0," + (-radius) + "), rotate(" + d + ", 0, " + radius + ")";
     })
     .attr("d", d3.svg.line()
                  .x(function(d) { return d.x; })
                  .y(function(d) { return d.y; })
		  .interpolate("basis")(dataGp120)
		  )
    .attr("stroke", "darkred")
    .attr("stroke-width", 5)
    .attr("stroke-opacity", 0.7)
    .style("fill", "none");

 var genomeData = [
  {"x": -50, "y": 0},
  {"x": -20, "y": 5},
  {"x": 10, "y": -5},
  {"x": 40, "y": 0},
 ];
 vir.append("path")
    .attr("transform", "rotate(15), scale(0.8)")
    .attr("d", d3.svg.line()
                 .x(function(d) { return d.x; })
		 .y(function(d) { return d.y; })
		 .interpolate("basis")(genomeData))
    .attr("stroke", "white")
    .attr("stroke-width", 4);
  
 // draw genome
 vis.append("line")
    .attr("x1", width / 4 - 30)
    .attr("x2", 0.67 * width - 10)
    .attr("y1", height / 2 - 15)
    .attr("y2", 0.4 * height)
    .attr("stroke", "white")
    .attr("opacity", 0.5)
    .attr("stroke-dasharray", "7,3")
    .attr("stroke-width", 2);

 vis.append("line")
    .attr("x1", width / 4 + 40)
    .attr("x2", 0.67 * width + 300)
    .attr("y1", height / 2 + 10)
    .attr("y2", 0.4 * height + 5)
    .attr("stroke", "white")
    .attr("opacity", 0.4)
    .attr("stroke-dasharray", "7,3")
    .attr("stroke-width", 2);

 var gen = vis.append("g")
   .attr("class", "genome")
   .attr("transform", "translate(" + (2 * width / 3) + "," + (0.4 * height) + ")");

 gen.append("line")
    .attr("x1", 0)
    .attr("x2", 300)
    .attr("y1", 0)
    .attr("y2", 0)
    .attr("stroke", "white")
    .attr("stroke-width", 8);

 var dataMuts = [50, 250];
 gen.selectAll(".escapeMut")
    .data(dataMuts)
    .enter()
    .append("circle")
    .attr("class", "escapeMut")
    .attr("cx", function(d) { return d; })
    .attr("cy", 0)
    .attr("r", 15)
    .attr("fill", "steelblue");



}

function updateVirion(active) {
 var gen = d3.select(".genome");
 console.log(gen);

 if (active) {

 var dataMuts = [20, 80, 100, 130, 180, 200, 220, 280];
 gen.selectAll(".otherMut")
    .data(dataMuts)
    .enter()
    .append("circle")
    .attr("class", "otherMut")
    .attr("cx", function(d) { return d; })
    .attr("cy", 0)
    .attr("r", 8)
    .attr("fill", "darkred")
    .style("opacity", 0)
    .transition()
    .duration(400)
    .style("opacity", 1);

 } else {
  gen.selectAll(".otherMut")
   .transition()
   .duration(400)
   .style("opacity", 0)
   .remove();
 }

}

function createMutSpread() {

 var secWidth = $("section").width();

 var margin = {left: 30, right:30, top:30, bottom:120},
     width = secWidth - margin.left - margin.right,
     height = 340,
     radius = 80,
     radiusCapsid = radius - 20;

 var vis = d3.selectAll("#mutSpread")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

 var dataVirion = [
  { "x": 0.25 * width, "y": 0.5 * height, "mut": true}, 
  { "x": 0.15 * width, "y": 0.82 * height, "mut": false}, 
  { "x": 0.07 * width, "y": 0.28 * height, "mut": false}, 
  { "x": 0.85 * width, "y": 0.5 * height, "mut": true}, 
  { "x": 0.75 * width, "y": 0.82 * height, "mut": true}, 
  { "x": 0.67 * width, "y": 0.28 * height, "mut": true}, 
 ];
 var vir = vis.selectAll("g.smallVirion")
   .data(dataVirion)
   .enter()
   .append("g")
   .attr("class", "smallVirion")
    .attr("id", function(d, i) { return "smallVirion" + (i+1); })
   .attr("transform", function(d) {
    return "translate(" + d.x + "," + d.y + "), scale(0.6)";
   });

 // draw virus circle
 vir.append("circle")
   .attr("cx", 0)
   .attr("cy", 0)
   .attr("r", radius)
   .style("fill", "steelblue")
   .style("fill-opacity", 0.3)
   .style("stroke", "steelblue")
   .style("stroke-width", 8);
 
 
 // draw capsid
 var hexagonData = [
       { "x": radiusCapsid/7,  "y": -radiusCapsid},
       { "x": radiusCapsid/3.5,  "y": -0.9 * radiusCapsid},
       { "x": radiusCapsid/2,   "y": 0.8 * radiusCapsid}, 
       { "x": radiusCapsid/3,   "y": radiusCapsid}, 
       { "x": -radiusCapsid/3,  "y": radiusCapsid},
       { "x": -radiusCapsid/2,  "y": 0.8 * radiusCapsid},
       { "x": -radiusCapsid/3.5,  "y": -0.9 * radiusCapsid},
       { "x": -radiusCapsid/7,  "y": -radiusCapsid},
     ];
 
 drawHexagon = d3.svg.line()
       .x(function(d) { return d.x; })
       .y(function(d) { return d.y; })
       .interpolate("linear-closed");
 
 vir.append("path")
    .attr("transform", "rotate(100)")
    .attr("d", drawHexagon(hexagonData))
    .attr("stroke", "darkred")
    .attr("stroke-width", 5)
    .style("fill", "darkred")
    .style("fill-opacity", 0.7);

  // HIV membrane protein
  var angleGp120 = [0, 25, 50, 75, 100, 125, 150, 175, 200, 226, 253,
                    281, 308, 334];
  var dataGp120 = [
    { "x": -4, "y": 6} ,
    { "x": -4, "y": -8} ,
    { "x": -8, "y": -12} ,
    { "x": 0, "y": -16} ,
    { "x": 8, "y": -12} ,
    { "x": 4, "y": -8} ,
    { "x": 4, "y": 6} ,
  ];
  vir.selectAll(".spike")
     .data(angleGp120)
     .enter()
     .append("path")
     .attr("transform", function(d) {
       return "translate(0," + (-radius) + "), rotate(" + d + ", 0, " + radius + ")";
     })
     .attr("d", d3.svg.line()
                  .x(function(d) { return d.x; })
                  .y(function(d) { return d.y; })
		  .interpolate("basis")(dataGp120)
		  )
    .attr("stroke", "darkred")
    .attr("stroke-width", 5)
    .attr("stroke-opacity", 0.7)
    .style("fill", "none");

 var genomeData = [
  {"x": -50, "y": 0},
  {"x": -20, "y": 5},
  {"x": 10, "y": -5},
  {"x": 40, "y": 0},
 ];
 vir.append("path")
    .attr("class", "littleGenome")
    .attr("id", function(d, i) { return "littleGenome" + (i+1); })
    .attr("transform", "rotate(15), scale(0.8)")
    .attr("d", d3.svg.line()
                 .x(function(d) { return d.x; })
		 .y(function(d) { return d.y; })
		 .interpolate("basis")(genomeData))
    .attr("stroke", "white")
    .attr("stroke-width", 7);

 vir.filter(function(d) { return d.mut == true; })
  .append("circle")
  .attr("cx", 10)
  .attr("cy", 1)
  .attr("r", 10)
  .attr("fill", "steelblue");

// line of evolution
vis.append("defs").append("marker")
    .attr("id", "arrowheadEvo")
    .attr("refX", 1) /*must be smarter way to calculate shift*/
    .attr("refY", 2)
    .attr("markerWidth", 6)
    .attr("markerHeight", 4)
    .attr("orient", "auto")
    .append("path")
        .attr("d", "M 0,0 V4 L6,2 Z")
	.attr("fill", "white");

 var dataLines = [
 {'x1': dataVirion[0].x + 1.2 * radius,
  'x2': dataVirion[3].x - 1.2 * radius,
  'y1': dataVirion[0].y,
  'y2': dataVirion[3].y},
 {'x1': dataVirion[0].x + 1.2 * radius,
  'x2': dataVirion[4].x - 1.2 * radius,
  'y1': dataVirion[0].y,
  'y2': dataVirion[4].y - 33},
 {'x1': dataVirion[0].x + 1.2 * radius,
  'x2': dataVirion[5].x - 1.2 * radius,
  'y1': dataVirion[0].y,
  'y2': dataVirion[5].y + 27},
 ];
 vis.selectAll(".evoLine")
    .data(dataLines)
    .enter()
    .append("line")
    .attr("marker-end", "url(#arrowheadEvo)")
    .attr("x1", function(d) { return d.x1; })
    .attr("x2", function(d) { return d.x2; })
    .attr("y1", function(d) { return d.y1; })
    .attr("y2", function(d) { return d.y2; })
    .attr("stroke-width", 4)
    .attr("stroke", "white");

 var dataFracs = [
 {'x': 0.12 * width, 'text': '33%'},
 {'x': 0.72 * width, 'text': '100%'},
 ];
 vis.selectAll(".fracs")
    .data(dataFracs)
    .enter()
    .append("text")
    .attr("x", function(d) { return d.x; })
    .attr("y", height + 60)
    .style("fontsize", "30px")
    .style("fill", "white")
    .text(function(d) { return d.text; });

}


function createDiversity() {

 var secWidth = $("section").width();

 var margin = {left: 180, right:30, top:80, bottom:120},
     width = secWidth - margin.left - margin.right,
     height = 240, roi = [6500, 7500];

 var vis = d3.selectAll("#divDivEnv")
   .attr("width", width + margin.left + margin.right)
   .attr("height", height + margin.top + margin.bottom)
   .append("g")
   .attr("class", "virion")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var y = d3.scale.log()
       .domain([0.001, 0.1])
       .range([height, 0]);
  
  var x = d3.scale.linear()
       .domain(roi)
       .range([20, width - 20]);
  
  // Axis stuff
  var xAxis = d3.svg.axis()
      .scale(x)
      .orient("bottom");
  
  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left");

  // X axis
  vis.append("g")
       .attr("class", "d3-axis-white")
       .attr("transform", "translate(0," + height + ")")
       .call(xAxis)
       .append("text")
       .attr("x", width / 2)
       .attr("y", 80)
       .style("text-anchor", "middle")
       .text("Position  in HIV genome [bp]");
   
  // Y axis
  vis.append("g")
       .attr("class", "d3-axis-white")
       .call(yAxis)
       .append("text")
       .attr("transform", "rotate(-90)")
       .attr("dy", "-4em")
       .attr("x", -height / 2)
       .style("text-anchor", "middle")
       .text("Diversity");

  // Diversity
  d3.xhr("/divdiv_local/")
  .header("Content-Type", "application/json")
  .post(
    JSON.stringify({patient: "p1",
                    observables: ["ds"],
		    itimes: [4],
		    roi: ["genomewide", roi[0], roi[1]]}),
    function(error, request) {
     var dataDiv = JSON.parse(request.responseText).data;

     vis.append("path")
        .attr("d", d3.svg.line()
                 .x(function(d, i) { return x(roi[0] + i * dataDiv.dx); })
           	 .y(function(d) { return y(d); })
           	 .interpolate("basis")
           	 (dataDiv.ds[0][1])
                    )
        .attr("stroke", "steelblue")
        .attr("stroke-width", 5)
        .attr("fill", "none");

     var barPos = [6760, 7060];
     var divBarGroup = vis.append("g")
	.datum(barPos)
        .attr("class", "divBarGroup")
	.attr("transform", function(d) { return "translate(" + x(d[0]) + ", " + (height - 25) + ")"; });

     divBarGroup.append("rect")
       .attr("id", "diversityBar")
       .attr("x", 0)
       .attr("y", 0)
       .attr("width", function(d) {return x(d[1]) - x(d[0]); })
       .attr("height", 10)
       .attr("fill", "darkred");
    
     divBarGroup.append("line")
      .attr("x1", 2)
      .attr("x2", 2)
      .attr("y1", 0)
      .attr("y2", 25 - height)
      .attr("stroke", "darkred")
      .attr("stroke-width", 4)
      .attr("opacity", 0.3)
      .style("stroke-dasharray", "7, 3");
      
     divBarGroup.append("line")
      .attr("x1", function(d) { return x(d[1]) - x(d[0]) - 2; })
      .attr("x2", function(d) { return x(d[1]) - x(d[0]) - 2; })
      .attr("y1", 0)
      .attr("y2", 25 - height)
      .attr("stroke", "darkred")
      .attr("stroke-width", 4)
      .attr("opacity", 0.3)
      .style("stroke-dasharray", "7, 3");

    });

}

function updateDiversity(active) {

 var secWidth = $("section").width();

 var margin = {left: 180, right:30, top:80, bottom:120},
     width = secWidth - margin.left - margin.right,
     height = 240, roi = [6500, 7500];

 var x = d3.scale.linear()
       .domain(roi)
       .range([20, width - 20]);

 if (active) {
  d3.select("g.divBarGroup")
	.datum([7140])
	.transition()
	.duration(1000)
	.attr("transform", function(d) { return "translate(" + x(d[0]) + ", " + (height - 25) + ")"; });

 } else {
  d3.select("g.divBarGroup")
	.datum([6760])
	.transition()
	.duration(1000)
	.attr("transform", function(d) { return "translate(" + x(d[0]) + ", " + (height - 25) + ")"; });

 }

}


/* Fragment events */
Reveal.addEventListener('fragmentshown', function(event) {
   // event.fragment = the fragment DOM element
   if (event.fragment.id == 'otherMuts') updateVirion(true);
   else if (event.fragment.id == 'divRRE') updateDiversity(true);
});
Reveal.addEventListener('fragmenthidden', function(event) {
   if (event.fragment.id == 'otherMuts') updateVirion(false);
   else if (event.fragment.id == 'divRRE') updateDiversity(false);
});


/* Load everything for now */
createPhysio1();
createPhysio2();
createPhysio3();
createVirion();
createMutSpread();
createDiversity();
