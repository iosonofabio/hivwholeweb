/* 
 * Plot phylogenetic trees with radial or rectangular representation
 * Author: Fabio Zanini
 */
// NOTE: this function is ready for callbacks, they have to set data = {id: <id>, chartType: <chartType>}
function update(data) {
    var svg = d3.select("#treeSvg"+data.id),
        divWidth = +($("#treeDiv"+data.id).width());

    // if this function is called with some useful data, bind it to the DOM
    if ("tree" in data)
        svg.datum(data)

    var chart = treeChart().svgWidth(0.9 * divWidth);

    // Figure out settings from the DOM, e.g. coloring and chart type
    if ($("#switchRadial").hasClass("active")) {
        chart.chartType("radial");
        chart.svgHeight(0.9 * divWidth);
    } else
        chart.svgHeight(15 * getNumberTerminals(svg.datum().tree, 0));

    if ($("#switchColorLinkDate").hasClass("active"))
        chart.colorLinkType("date");


    svg.call(chart);

    function getNumberTerminals(n, number) {
        if (n.children)
            for(var i=0; i < n.children.length; i++)
                number += getNumberTerminals(n.children[i], 0);
        else
            number += 1;
        return number;
    }

}

/* tree chart closure as of Mark Bostock: http://bost.ocks.org/mike/chart/ */
function treeChart() {

    var svgWidth = 400,
        svgHeight = 400,
        margin = {top: 5, bottom: 5, left: 5, right: 5},
        width = svgWidth - margin.left - margin.right,
        height = svgHeight - margin.top - margin.bottom,
        chartType = "rectangular",
        colorLinkType = "black";

    // TREE CHART FUNCTION
    function chart(selection) {
        selection.each(function (data) {

            var colorMap = d3.scale.linear()
                .interpolate(d3.interpolateRgb)
                .domain([0, 0.25, 0.33, 0.5, 0.67, 0.75, 1])
                .range(["darkblue", "blue", "cyan", "green",
                        "yellow", "orange", "red"]);

            var tree = data.tree,
                id = data.id,
                pname = data.pname,
                region = data.region;
           
            // Create outer chart (SVG) and make sure there are no other ones
            var svg = d3.select(this);
            svg.selectAll("*").remove();

            // Set the outer dimensions.
            //svg.attr("width", width + margin.left + margin.right)
            //   .attr("height", height + margin.top + margin.bottom)
            //responsive SVG needs these 2 attributes and no width and height attr
            svg.attr("preserveAspectRatio", "xMinYMin meet")
               .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom));


            // calculate various properties of all nodes
            setDepths(tree, 0);
            setNTerminals(tree);
            setDSI(tree);

            var depth = getMaxTree(tree, function(d) { return d.depthScaled; });

            // set coloring functions
            if (colorLinkType == "black") {
                var colorLinkFunc = function(d) { return "black"; }
            
            } else {
                var dsiMax = getMaxTree(tree, function(d) { return d.DSI; }),
                    colorLinkFunc = function(d) {
                        var dsi = d.target.DSI;
                        if (dsi != "undefined")
                            return colorMap(dsi / dsiMax);
                        else
                            return "grey";
                    };
            }

            // tooltip
            var tip = d3.tip()
                .attr('class', 'd3-tip')
                .html(tooltipFunc);

            // prepare cluster representation
            var cluster = d3.layout.cluster()
               .sort(null)
               .value(function(d) { return d.branch_length; })
               .separation(function(a, b) { return 1; });

            // plot the chart
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


                // activate tip on visualized svg
                vis.call(tip);

                cluster.size([maxAngle, 1]);

                var nodes = cluster.nodes(tree),
                    treeScale = 0.9 * rInternal / depth;

                // adjust the bar to calibration
                var barLengthData = (30.0 / treeScale).toPrecision(1),
                    barLength = treeScale * barLengthData;

                // add scaled depths (in pixels) to the tree
                nodes.map(function(n) {n.y = n.depthScaled * treeScale; });
           
                // links
                var link = vis.selectAll("path.link")
                     .data(cluster.links(nodes))
                     .enter()
                     .append("path")
                     .attr("class", "link")
                     .attr("d", stepRadial)
                     .attr("fill", "none")
                     .attr("stroke", colorLinkFunc)
                     .attr("stroke-width", 2)
                     .on("mouseover", moverLinksRadial)
                     .on("mouseout", moutLinksRadial);
           
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
                      d.name = "" + d.name;
                      // local trees have the sequences attached
                      if (region.indexOf('minor') == -1) {
                       return d.name.replace(/_/g, ' ');
                      } else {
                    return d.name.split('_').slice(1).join(" ");
                      }
                     })
                     .on("mouseover", moverRadial)
                     .on("mouseout", moutRadial);
           
                function moverLinksRadial(d) {
                    tip.show(d);
                    moverRadial(d.target);
                }
           
                function moutLinksRadial(d) {
                    tip.hide(d);   
                    moutRadial(d.target);
                }

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

                // activate tip on visualized svg
                vis.call(tip);

                // note: at present, x and y are swapped to keep consistency with the radial layout
                cluster.size([height, 0.85 * width]);
           
                var nodes = cluster.nodes(tree),
                    treeScale = cluster.size()[1] / depth;

                // add scaled depths (in pixels) to the tree
                nodes.map(function(n) {n.y = n.depthScaled * treeScale; });

                // adjust the bar to calibration
                var barLengthData = (30.0 / treeScale).toPrecision(1),
                    barLength = treeScale * barLengthData;

                // links
                var link = vis.selectAll("path.link")
                     .data(cluster.links(nodes))
                     .enter()
                     .append("path")
                     .attr("class", "link")
                     .attr("d", stepRectangular)
                     .attr("fill", "none")
                     .attr("stroke", colorLinkFunc)
                     .attr("stroke-width", 2)
                     .on("mouseover", moverLinksRectangular)
                     .on("mouseout", moutLinksRectangular);

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

                // line connecting the leaves to their labels
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
           
                function moverLinksRectangular(d) {
                    tip.show(d);
                    moverRectangular(d.target);
                }
           
                function moutLinksRectangular(d) {
                    tip.hide(d);
                    moutRectangular(d.target);
                }


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

        // NOTE: node.depth is a d3 construct to get the integer depth
        // so we use depthScaled
        function setDepths(n, parentDepth) {
            n.depthScaled = parentDepth + n.branch_length;
            if (n.children)
                for (var i=0; i < n.children.length; i++)
                    setDepths(n.children[i], n.depthScaled);        
        }

        function getMaxTree(n, accessor) {
            if (n.children)
                return d3.max(n.children, function (d) { return getMaxTree(d, accessor); });
            else
                return accessor(n);
        }

        function setNTerminals(n) {
            if (n.children) {
                n.children.map(setNTerminals);
                n.nTerminals = d3.sum(n.children, function(d) { return d.nTerminals; });
            } else
                n.nTerminals = 1;
        }

        // Get weighted mean of DSIs of terminal nodes
        function setDSI(n) {
            if (n.children) {
                n.children.map(setDSI);
                var dsicum = d3.sum(n.children, function (d) {
                    var dsi = d.DSI;
                    if (dsi == "undefined")
                        return 0;
                    return dsi * d.nTerminals;
                });
                var childrencum = d3.sum(n.children, function(d) {
                    if (d.DSI == "undefined")
                        return 0;
                    return d.nTerminals;
                });
                if (childrencum == 0)
                    n.DSI = "undefined";
                else
                    n.DSI = dsicum / childrencum;
            }
        }

        // Tooltip function
        function tooltipFunc(d) {
            var n = d.target,
                msg = "";

            if (!(n.children)) {
                var pname =  n.name.split('_')[0];
                if (pname[0] == "p")
                    msg = msg + "Patient: " + pname + "</br>";
            
                if (isNumeric(n.CD4))
                    msg = msg + "CD4+ cell count [cells/ml]: " + n.CD4 + "</br>";

                if (isNumeric(n.VL))
                    msg = msg + "Viral load [virions/ml]: " + n.VL + "</br>";

            }   

            if (isNumeric(n.DSI))
                msg = msg + "Day since infection: " + n.DSI.toFixed(0) + "</br>";

            msg = msg + "Mutations on this branch: ";
            if (n.muts.length > 0) {
                var muts = n.muts.split(" "),
                    nMuts = muts.length,
                    nMutsPerLine = 10,
                    nMutLines = Math.ceil(nMuts / nMutsPerLine);
                
                if (nMutLines == 1)
                    msg = msg + n.muts;
                else {
                    msg = msg + "</br>";
                    for(var i=0; i < nMutLines - 1; i++)
                        msg = msg + muts.slice(i * nMutsPerLine, (i+1) * nMutsPerLine).join(" ") + "</br>";
                    msg = msg + muts.slice((nMutLines - 1) * nMutsPerLine, nMuts).join(" ");
                
                }
            } else
                msg = msg + "(none)"

            return msg;
        }

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

    chart.colorLinkType = function (_) {
        if (!arguments.length) return colorLinkType;
        colorLinkType = _;
        return chart;
    };

    return chart;
}
