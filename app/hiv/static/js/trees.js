function emptyTree(id, keepData) {

 var svg = d3.select('#'+id);
 svg.selectAll("*").remove();

 if ((typeof(keepData) == "undefined") | (!keepData)) svg.datum(null);

}

/* 
 * Plot phylogenetic trees with radial or rectangular representation
 * Author: Fabio Zanini
 */
// NOTE: this function is ready for callbacks, they have to set data = {chartType: <chartType>}
function updateTree(id, data) {
    if (typeof(data)==='undefined') data = {};

    var svg = d3.select('#'+id),
        divWidth = $('#'+id).parent().width();

    // if this function is called with some useful data, bind it to the DOM
    if ("tree" in data)
        svg.datum(data)

    var chart = treeChart().svgWidth(0.9 * divWidth);

    // Figure out settings from the DOM, e.g. coloring and chart type
    if ($("#switchRectangular").hasClass("active")) {
        chart.chartType("rectangular");
        chart.svgHeight(15 * getNumberTerminals(svg.datum().tree, 0));
    } else
        chart.svgHeight(0.9 * divWidth);

    if ($("#switchColorLinkDate").hasClass("active"))
        chart.colorLinkType("date");
    else if ($("#switchColorLinkSubtype").hasClass("active"))
        chart.colorLinkType("subtype");
    else if ($("#switchColorLinkPatient").hasClass("active"))
        chart.colorLinkType("patient");

    // Leaf labels defaults to false for minor variant trees
    if ((data.leafLabels === false) || ((typeof(data.leafLabels) == "undefined") && (svg.datum().region.indexOf("minor") != -1)))
        chart.leafLabels(false);

    if (data.optimizeSpace === true)
        chart.optimizeSpace(true);

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
        chartType = "radial",
        colorLinkType = "black",
        leafLabels = true,
        optimizeSpace = false;

    // TREE CHART FUNCTION
    function chart(selection) {
        selection.each(function (data) {

            var colorMap = d3.scale.linear()
                .interpolate(d3.interpolateRgb)
                .domain([0, 0.25, 0.33, 0.5, 0.67, 0.75, 1])
                .range(["darkblue", "blue", "cyan", "green",
                        "yellow", "orange", "red"]);
            
            var tree = data.tree,
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
            setPatientInternal(tree);

            var depth = Math.max(1e-7, getMaxTree(tree, function(d) { return d.depthScaled; })),
                dsiMax = getMaxTree(tree, function(d) { return d.DSI; });

            // set link coloring function
            if (colorLinkType == "date")
                var colorLinkFunc = function(d) {
                    var dsi = d.target.DSI;
                    if (dsi != "undefined")
                        return colorMap(dsi / dsiMax);
                    else
                        return "grey";
                }
            else if (colorLinkType == "subtype")
                var colorLinkFunc = function(d) {
                    var n = d.target;
                    // NOTE: colors from colorbrewer dark2
                    if (n.subtype === 'B')
                        return "#7570b3";
                    else if (n.subtype === 'C')
                        return "#d95f02";
                    else if (n.subtype === 'AE')
                        return "#1b9e77"
                    else
                        return "grey";
                }
            else if (colorLinkType == "patient")
                var colorLinkFunc = function(d) {
                    var cscale = d3.scale.category20()
                        .domain(['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11']);
                    if (d.target.patient != "undefined")
                        return cscale(d.target.patient);
                    else if (d.target.name != "undefined")
                        return "black";
                    else
                        return "grey";
                }
            else
                var colorLinkFunc = function(d) { return "black"; }

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
                    treeCenter = {'cx': r, 'cy': r};

                // if the leaf labels are not shown, center the tree differently
                if (leafLabels === true) {
                    var rInternal = r - 170;
                } else {
                    var rInternal = r;
                }

                // adjust the bar to calibration
                var treeScale = 0.9 * rInternal / depth;

                // set up the d3 cluster
                cluster.size([maxAngle, 1]);
                var nodes = cluster.nodes(tree);

                // add position attributes to the tree:
                // - angle
                // - radius from the coordinate center (treeCenter)
                // - x coordinate from the treeCenter
                // - y coordinate from the treeCenter
                nodes.map(function(n) {
                    n.radius = n.depthScaled * treeScale;
                    n.angle = n.x;
                    setProjectRadial(n);
                });

                // if requested, optimize space by resetting the viewbox
                if (optimizeSpace === true) {
                    // get max and min x and y
                    var xMax = d3.max(nodes, function(d){ return d.x; }),
                        xMin = d3.min(nodes, function(d){ return d.x; }),
                        yMax = d3.max(nodes, function(d){ return d.y; }),
                        yMin = d3.min(nodes, function(d){ return d.y; });

                    svg.attr("viewBox", ((r + xMin - 20) + " " + 
                                         (r + yMin - 20) + " " + 
                                         (xMax - xMin + 40 + margin.left + margin.right) + " " +
                                         (yMax - yMin + 40 + margin.top + margin.bottom))
                            );

                }

                // SVG group to render tree in
                var vis = svg.append("g")
                             .attr("transform", ("translate(" +
                                                 (margin.top + treeCenter.cy) + "," +
                                                 (margin.left + treeCenter.cx) + ")")
                                  );

                //// Test dot in the treeCenter
                //svg.append("circle")
                //    .attr("cx", margin.left + r)
                //    .attr("cy", margin.top + r)
                //    .attr("r", 10)
                //    .style("fill", "red");

                //svg.append("rect")
                //    .attr("x", margin.left)
                //    .attr("y", margin.top)
                //    .attr("width", 2 * r)
                //    .attr("height", 2 * r)
                //    .attr("fill", "none")
                //    .attr("stroke", "red")
                //    .attr("stroke-width", 3);

                // activate tip on visualized svg
                vis.call(tip);

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


                // balls proportional to frequency for minor variants
                if (region.indexOf("minor") != -1) {
                    link.filter(function(d) { return d.target.frequency != "undefined" })
                        .each(plotBallsRadial);

                }
           
                // line connecting the leaves to their labels
                if (leafLabels === true) {
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
                         .attr("text-anchor", function(d) { return d.angle < 180 ? "start" : "end"; })
                         .attr("transform", function(d) {
                             return ("rotate(" + (d.angle - 90) + 
                                     ")translate(" + (r - 170 + 8) + 
                                     ")rotate(" + (d.angle < 180 ? 0 : 180) + ")");
                         })
                         .text(leafLabelFunc)
                         .on("mouseover", moverRadial)
                         .on("mouseout", moutRadial);
                }
           
                // scale bar
                if ((leafLabels === true) & (depth > 1e-6)) {
                    var barLengthData = (30.0 / treeScale).toPrecision(1),
                        barLength = treeScale * barLengthData,
                        barLengthText = String(barLengthData);

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
                       .text(barLengthText)
                       .attr("text-anchor", "middle");
                }

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
                      .attr("cx", t.x)
                      .attr("cy", t.y)
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
                    var r = d.radius,
                        a = (d.angle - 90) / 180 * Math.PI;
                    return {'x': r * Math.cos(a), 'y': r * Math.sin(a)};
                }

                function setProjectRadial(d) {
                    var pr = projectRadial(d);
                    d.x = pr.x;
                    d.y = pr.y;         
                }
                
                function stepRadial(d) {
                  var s = d.source,
                      m = projectRadial({angle: d.target.angle, radius: d.source.radius}),
                      t = d.target,
                      r = d.source.radius,
                      sweep = d.target.angle > d.source.angle ? 1 : 0;
                  return (
                    "M" + s.x + "," + s.y +
                    "A" + r + "," + r + " 0 0," + sweep + " " + m.x + "," + m.y +
                    "L" + t.x + "," + t.y);
                }
                
                function stepAnnoRadial(d, r) {
                    var s = projectRadial({angle: d.angle, radius: d.radius + 5}),
                        t = projectRadial({angle: d.angle, radius: r});
                    return (
                        "M" + s.x + "," + s.y +
                        "L" + t.x + "," + t.y);
                }

                function plotBallsRadial(d, i) {
                    var freq = d.target.frequency,
                        r = 2 + (18 - 5) * ((Math.log(freq) / Math.LN10) + 2.0) / 2.0;
                    vis.append("circle")
                        .datum(d.target)
                        .attr("class", "leaf")
                        .attr("r", r)
                        .attr("cx", d.target.x)
                        .attr("cy", d.target.y)
                        .attr("fill", colorLinkFunc(d, i))
                        .attr("fill-opacity", 0.7);
                
                }

            }

            // RECTANGULAR CHART
            function makeRectangular() {

                var vis = svg.append("g")
                             .attr("transform", "translate(" + (margin.top) + "," + (margin.left) + ")");

                // activate tip on visualized svg
                vis.call(tip);

                // set up d3 cluster
                // note: at present, x and y are swapped to keep consistency with the radial layout
                cluster.size([height, 0.93 * width]);

                // adjust the bar to calibration
                var treeScale = cluster.size()[1] / depth;

                if ((leafLabels === true) & (depth > 1e-6)) {
                    var barLengthData = (30.0 / treeScale).toPrecision(1),
                        barLength = treeScale * barLengthData;

                    // set the bar text, should be the same as the data except in degenerate trees
                    var barLengthText = String(barLengthData);

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
                       .text(barLengthText)
                       .attr("text-anchor", "middle");
                }

                var nodes = cluster.nodes(tree);

                // add scaled depths (in pixels) to the tree
                nodes.map(function(n) {n.y = n.depthScaled * treeScale; });
                
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

                // balls proportional to frequency for minor variants
                if (region.indexOf("minor") != -1) {
                    link.filter(function(d) { return d.target.frequency != "undefined" })
                        .each(plotBallsRect);

                }

                // line connecting the leaves to their labels
                if (leafLabels === true) {
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
                         .text(leafLabelFunc)
                         .on("mouseover", moverRectangular)
                         .on("mouseout", moutRectangular);
                }
           
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
                        "H" + d.target.y
                    );
                }

                function stepAnnoRectangular(d) {
                    // FIXME: the depth is a trick, we should allocate space for the labels and subtract 
                    return (
                        "M" + d.y + "," + d.x +
                        "H" + (10 + depth * treeScale)
                    );
                }


                function plotBallsRect(d, i) {
                    var freq = d.target.frequency,
                        r = 2 + (18 - 5) * ((Math.log(freq) / Math.LN10) + 2.0) / 2.0;
                    vis.append("circle")
                        .attr("class", "leaf")
                        .attr("r", r)
                        .attr("cx", d.target.y)
                        .attr("cy", d.target.x)
                        .attr("fill", colorLinkFunc(d, i));
                
                }

            }

            function leafLabelFunc(d) {
                if (String(region).indexOf('minor') != -1) {
                    var freqLabel = String((d.frequency * 100).toFixed(0)) + "%",
                        timeLabel = String(Math.floor(d.DSI / 30.5)) + " m",
                        label = freqLabel + ", " + timeLabel;
                    
                    return label;
                } else {
                    var name = String(d.name);

                   return name.split('_').join(" ");
                }
            }

        });

        // NOTE: node.depth is a d3 construct to get the integer depth, we use depthScaled
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

        // Set the patient designation of internal nodes to that of its children
        // iff all agree. otherwise, set to undefined
        function setPatientInternal(n) {
            if (n.children) {
                n.children.map(setPatientInternal);
                var patient = n.children.map(function (d) {return d.patient});
                if (patient.every(function (d){return d==patient[0]})){
                    n.patient=patient[0];
                }else{
                    n.patient="undefined";
                }
            }
        }


        // Tooltip function
        function tooltipFunc(d) {
            var n = d.target,
                msg = "";

            if (!(n.children)) {
                var pname =  String(n.name).split('_')[0];
                if (pname[0] == "p")
                    msg = msg + "Patient: " + pname + "</br>";
                else if (pname != "undefined") {
                    if (pname == "LAI-III")
                        msg = msg + "Name: " + pname + " (HXB2)</br>";
                    else
                        msg = msg + "Name: " + pname + "</br>";
                }
            
                if (isNumeric(n.CD4))
                    msg = msg + "CD4+ cell count [cells/ml]: " + n.CD4 + "</br>";

                if (isNumeric(n.VL))
                    msg = msg + "Viral load [virions/ml]: " + n.VL + "</br>";

            }

            if ((typeof(n.subtype) != "undefined") && (n.subtype !== "undefined")) {
                msg = msg + "Subtype: " + n.subtype + "</br>";
            }

            if (isNumeric(n.DSI))
                msg = msg + "Day since infection: " + n.DSI.toFixed(0) + "</br>";

            if (isNumeric(n.frequency))
                msg = msg + "Frequency: " + (100 * n.frequency).toFixed(0) + "%</br>";

            if ((typeof(n.muts) != "undefined") && (n.muts !== "undefined")){
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
            }

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

    chart.leafLabels = function (_) {
        if (!arguments.length) return leafLabels;
        leafLabels = _;
        return chart;
    };

    chart.optimizeSpace = function (_) {
        if (!arguments.length) return optimizeSpace;
        optimizeSpace = _;
        return chart;
    };

    return chart;
}
