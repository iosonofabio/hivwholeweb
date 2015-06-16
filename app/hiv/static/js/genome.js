function emptyGenome(id, keepData) {

 var svg = d3.select('#'+id);
 svg.selectAll("*").remove();

 if ((typeof(keepData) == "undefined") | (!keepData)) svg.datum(null);

}


// FIXME: some genes overlap (p2/p6?)
function updateGenome(id, data) {

    var svg = d3.select('#'+id),
        divWidth = $('#'+id).parent().width();

    var chart = genomeChart().svgWidth(divWidth);

    svg.datum(data)
       .call(chart);

}

/* genome chart closure */
function genomeChart() {

    var svgWidth = 400,
        svgHeight = 290,
        margin = {top: 3, right: 60, bottom: 42, left: 60},
        width = svgWidth - margin.left - margin.right,
        height = svgHeight - margin.top - margin.bottom,
        vpadBlockTop = 5,
        vpadBlock = 5,
        nBlocks = 6,
        fontSize = 5; // px
        heightBlock = ((height - vpadBlockTop) / nBlocks) - vpadBlock;
        // NOTE: each group might take more than one block in height
        featureHierarchy = {
            "fragment": 0,
            "gene": 2,
            "protein": 2,
            "RNA_structure": 5,
            "other": 5,
        },
        drawBorderTop = true,
        drawBorderLeft = true,
        drawBorderRight = true,
        resizeSvg = true,
        x = "undefined",
        /* pre, post callbacks for reusability */
        zoomCallbacks = {};

    function chart(selection) {
        // GENOME CHART FUNCTION
        selection.each(function (data) {

            // get the genome data
            var genome = data.genome,
                features = genome.features;
    
            // tooltip (needs to be initialized)
            var tip = d3.tip()
                 .attr('class', 'd3-tip')
                 .html(function(d) { return d.name + ": " + (+d.location[0][0] + 1) + ", " + d.location[0][1]; });
    
            // outer chart with genometric zoom upon resizing
            var svg = d3.select(this);
            if (resizeSvg)
                svg.attr("preserveAspectRatio", "xMinYMin meet")
                    .attr("viewBox", "0 0 " + svgWidth + " " + svgHeight);
      
            // internal chart, only the bars
            var vis = svg.append("g")
                 .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
            // add tooltip
            vis.call(tip);

            // set a standard x scale if not specified
            if (x == "undefined")
                x = d3.scale.linear()
                     .domain([-50, genome.len + 50])
                     .range([20, width - 20]);
    
            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");
        
            // add axis to the chart
            var xAxisObj = vis.append("g")
                              .attr("class", "d3-axis")
                              .attr("transform", "translate(0," + height + ")")
                              .call(xAxis)
    
            xAxisObj.append("text")
                    .attr("x", width / 2)
                    .attr("y", 40)
                    .style("text-anchor", "middle")
                    .text("Position [bp]");
    
            // add the rest of the rectangle box
            if (drawBorderLeft)
                xAxisObj.append("line")
                 .attr({"class": "bBox",
                        "x1": x.range()[0], "x2": x.range()[0],
                        "y1": 0, "y2": -height});
            if (drawBorderRight)
                xAxisObj.append("line")
                 .attr({"class": "bBox",
                        "x1": x.range()[1], "x2": x.range()[1],
                        "y1": 0, "y2": -height});
            if (drawBorderTop) {
                xAxisObj.append("line")
                 .attr({"class": "bBox",
                        "x1": x.range()[0], "x2": x.range()[1],
                        "y1": -height, "y2": -height});
            }

            // plot the labels for features on the left
            var featureLabels = [
                {
                 "label": "PCR",
                 "y": vpadBlockTop + (vpadBlock + heightBlock) * (featureHierarchy["fragment"] + 1),
                },
                {
                 "label": "RF1",
                 "y": vpadBlockTop + (vpadBlock + heightBlock) * featureHierarchy["gene"] + 0.5 * heightBlock,                
                },
                {
                 "label": "RF2",
                 "y": vpadBlockTop + (vpadBlock + heightBlock) * (featureHierarchy["gene"] + 1) + 0.5 * heightBlock,                
                },
                {
                 "label": "RF3",
                 "y": vpadBlockTop + (vpadBlock + heightBlock) * (featureHierarchy["gene"] + 2) + 0.5 * heightBlock,                
                },
                {
                 "label": "other",
                 "y": vpadBlockTop + (vpadBlock + heightBlock) * featureHierarchy["other"] + 0.5 * heightBlock,                
                },
            ];
            vis.selectAll(".featureLabel")
                .data(featureLabels)
                .enter()
                .append("text")
                .attr("x", -7)
                .attr("y", function(d) { return d.y; })
                .text(function(d) { return d.label; })
                .attr("dy", ".35em")
                .attr("text-anchor", "end");


            // plot the various bars
            plotAllFeatures();

            function plotAllFeatures() {
                plotFeatureGroup("fragment");
                plotFeatureGroup("protein");
                plotFeatureGroup("gene-single");
                plotFeatureGroup("gene-exon1");
                plotFeatureGroup("gene-exon2");
                plotFeatureGroup("RNA_structure");
                plotFeatureGroup("other");
            }

            function clickFeature(d) {
                var self = d3.select(this);
                
                // TODO: use transitions, but they have to be simultaneous which requires some Chinese-boxing

                if (!(self.classed("zoomed"))) {
                    var start = d.location[0][0],
                        end = d.location[0][1],
                        zoomData = {'start': start, 'end': end, 'name': d.name};

                    if (("zoomin" in zoomCallbacks) && ("pre" in zoomCallbacks.zoomin))
                        zoomCallbacks.zoomin.pre(zoomData);

                    // delete boxes fully outside of window
                    vis.selectAll(".featurebox")
                        .filter(function(dOther) { return (dOther.location[0][0] >= end) | (dOther.location[0][1] <= start); })
                        .remove();

                    // in case exon lines are being shown, stop
                    vis.selectAll(".exon-line").remove();

                    // change the x scale
                    var boxWidth = end - start,
                        scalePad = 0.005 * boxWidth;
                    x.domain([-scalePad + start, end + scalePad]);
                    xAxis.scale(x);
                    xAxisObj.call(xAxis);

                    if (("zoomin" in zoomCallbacks) && ("middle" in zoomCallbacks.zoomin))
                        zoomCallbacks.zoomin.middle(zoomData);

                    // move, resize, and recenter boxes
                    function zoomIn(d) {
                        var self = d3.select(this);
                            openLeft = d.location[0][0] < start,
                            openRight = d.location[0][1] > end,
                            startSelf = d3.max([start, d.location[0][0]]),
                            endSelf = d3.min([end, d.location[0][1]]);
                    
                        // move feature group
                        var boxy = +(self.attr("transform").split(",")[1].split(")")[0]);
                        self.attr("transform", "translate(" + x(startSelf) + "," + boxy + ")");

                        // resize rectangle
                        self.select("rect")
                            .attr("width", x(endSelf) - x(startSelf));

                        // leave borders open for features that are not fully in the zoom
                        // NOTE: this is a VERY DIRTY trick due to lacking SVG features
                        if (openLeft & (!openRight))
                            self.select("rect")
                                .attr("stroke-dasharray",
                                      self.select("rect").attr("width") +
                                      ", 0, " +
                                      self.select("rect").attr("height") +
                                      ", 0, " +
                                      self.select("rect").attr("width") +
                                      ", " +
                                      self.select("rect").attr("height")
                                     );
                        if (openRight & (!openLeft))
                            self.select("rect")
                                .attr("stroke-dasharray",
                                      self.select("rect").attr("width") +
                                      ", " +
                                      self.select("rect").attr("height") +
                                      ", " +
                                      self.select("rect").attr("width") +
                                      ", 0, " +
                                      self.select("rect").attr("height")
                                     );
                        if (openRight & openLeft)
                            self.select("rect")
                                .attr("stroke-dasharray",
                                      self.select("rect").attr("width") +
                                      ", " +
                                      self.select("rect").attr("height") +
                                      ", " +
                                      self.select("rect").attr("width") +
                                      ", " +
                                      self.select("rect").attr("height")
                                     );

                        // center text
                        self.select("text")
                            .attr("x", 0.5 * (x(endSelf) - x(startSelf)));

                        self.classed("zoomed", true);

                    }
                    vis.selectAll(".featurebox")
                        .each(zoomIn);

                    if (("zoomin" in zoomCallbacks) && ("post" in zoomCallbacks.zoomin))
                        zoomCallbacks.zoomin.post(zoomData);
    
                } else {
                    if (("zoomout" in zoomCallbacks) && ("pre" in zoomCallbacks.zoomout))
                        zoomCallbacks.zoomout.pre();
    
                    // Restoring could be done more stylish, but it's ok
                    vis.selectAll(".featurebox").remove();
                    
                    // restore the x scale
                    x.domain([-50, genome.len + 50]);
                    xAxis.scale(x);
                    xAxisObj.call(xAxis);

                    plotAllFeatures();

                    if (("zoomout" in zoomCallbacks) && ("post" in zoomCallbacks.zoomout))
                        zoomCallbacks.zoomout.post();
    
                }
            }        

            function moverFeature(d) {
                // change color
                var rect = d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect");
                
                rect.style("fill", "darkred");  
      
                // show tooltip
                tip.show(d);

                // show line connecting the exons
                var exons = [];
                if (d.name.indexOf('exon 1') != -1)
                    exons = ["exon-1", "exon-2"];
                else if (d.name.indexOf('exon 2') != -1)
                    exons = ["exon-2", "exon-1"];
    
                if (exons.length > 0) {
                    var rect2 = d3.select("#" + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "rect");
                    // within zoom views, the other exon might be missing
                    if (rect2.empty())
                        return;

                    var d2 = rect2.data()[0];
                    var box2 = d3.select('#' + d.name.replace(/[ ']/g, "-").replace(exons[0], exons[1]) + "-box");
                    // NOTE: we extract the y coordinate from the affine transformation STRING!
                    var boxy = +(d3.select(this).attr("transform").split(",")[1].split(")")[0]);
                    var box2y = +(box2.attr("transform").split(",")[1].split(")")[0]);
      
                    vis.append('line')
                        .attr("class", "exon-line")
                        .attr("x1", 0.5 * (x(d.location[0][1]) + x(d.location[0][0])))
                        .attr("x2", 0.5 * (x(d2.location[0][1]) + x(d2.location[0][0])))
                        .attr("y1", boxy + 0.5 * heightBlock)
                        .attr("y2", box2y + 0.5 * heightBlock)
                        .style("opacity", 0.7)
                        .style("stroke", "darkred")
                        .style("stroke-width", 2);

                }
            }
    
            function bgColor(d) {
                if (!(data.hasOwnProperty("highlightedRegions")))
                    return "steelblue";
                else if (data.highlightedRegions.indexOf(d.name) == -1)
                    return "steelblue";
                else
                    return "darkorange";
            }

            function moutFeature(d) {
                // change color back
                // if the feature has a native color, use it, otherwise steelblue
                var fearect = d3.select("#" + d.name.replace(/[ ']/g, "-") + "rect");

                fearect.style("fill", bgColor(d));

                // hide tooltip
                tip.hide(d);   
      
                // hide exon line
                if ((d.name.indexOf('exon 1') != -1) | (d.name.indexOf('exon 2') != -1))
                  vis.selectAll(".exon-line").remove();
            }
      
            function plotFeatureGroup(groupname) {
                var featureName = groupname.split("-")[0],
                    group = getFeatureGroup(groupname),
                    yGroup = vpadBlockTop + (vpadBlock + heightBlock) * 
                        featureHierarchy[featureName];

                var fea = vis.selectAll("." + groupname)
                    .data(group)
                    .enter()
                    .append("g")
                    .attr("class", "featurebox " + groupname)
                    .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "-box"; })
                    .attr("transform", function(d) { return "translate(" + x(d.location[0][0]) +
                          "," + (yGroup + heightFeature(d)) + ")"; })
                    .on('mouseover', moverFeature)
                    .on('mouseout', moutFeature)
                    .on("click", function() { d3.event.stopPropagation(); })
                    .on('click.zoom', clickFeature);
    
                var fearect = fea.append("rect")
                    .attr("class", "featurerect")
                    .attr("id", function(d) { return d.name.replace(/[ ']/g, "-") + "rect"; })
                    .attr("x", 0)
                    .attr("y", 0)
                    .attr("width", function(d) { return x(d.location[0][1]) - x(d.location[0][0]); })
                    .attr("height", heightBlock)
                    .style("fill", bgColor);

                // show text only of longer things (that's why we have a tooltip)
                // FIXME: use a criterion that is more zoom-friendly
                fea.append("text")
                    .attr("x", function(d) { return 0.30 * (x(d.location[0][1]) - x(d.location[0][0])); })
                    .attr("y", heightBlock / 2 + fontSize)
                    .text(function(d) {
                        if ((d.location[0][1] - d.location[0][0]) > 350)
                            return d.name;
                        else
                            return "";
                    });
            }

            function getFeatureGroup(groupname) {
                var group = [];
                var gname = "";
                var tmpfea;
                for(i = 0; i < features.length; i++) {
                 gname = features[i].type
                 if (gname == groupname) {
                  group.push(features[i]);
                 } else if (gname + '-single' == groupname) {
                  if (['vpu', 'vpr', 'nef', 'vif'].indexOf(features[i].name) != -1) {
                   group.push(features[i]);
                  }
                 } else if (gname + '-exon1' == groupname) {
                  if (['rev', 'tat'].indexOf(features[i].name) != -1) {
                   tmpfea = {"name": features[i].name + ' exon 1',
                             "type": "gene",
                             "location": [features[i].location[0]]};
                   group.push(tmpfea);
                  }
                 } else if (gname + '-exon2' == groupname) {
                  if (['rev', 'tat'].indexOf(features[i].name) != -1) {
                   tmpfea = {"name": features[i].name + ' exon 2',
                             "type": "gene",
                             "location": [features[i].location[1]]};
                   group.push(tmpfea);
                  }
                 }
                }
                return group;
            }
    
            function heightFeature(feature) {
                if (feature.type == "fragment")
                    return heightFragment(feature);
                else if ((feature.type == "protein") | (feature.type == "gene"))
                    return heightInframe(feature); 
                else
                    return 0;
    
                function heightFragment(fragment) {
                    var fn = +(fragment.name[1]);
                    var height = 0;
                    if ((fn % 2) == 0)
                        height += heightBlock + vpadBlock;
                    return height;
                }
    
                function heightInframe(feature) {
                    var frame;
                    if (['tat exon 2', 'rev exon 2'].indexOf(feature.name) != -1)
                        frame = (+(feature.location[0][1]) - (+genome.framestart)) % 3;
                   else
                     frame = (+(feature.location[0][0]) - (+genome.framestart)) % 3; 
                   return frame * (heightBlock + vpadBlock);
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
        if (!arguments.length) return svgWidth;
        svgWidth = _;
        width = svgWidth - margin.left - margin.right;
        return chart;
    };

    chart.svgHeight = function (_) {
        if (!arguments.length) return svgHeight;
        svgHeight = _;
        height = svgHeight - margin.top - margin.bottom;
        heightBlock = ((height - vpadBlockTop) / nBlocks) - vpadBlock;
        return chart;
    };

    chart.width = function (_) {
        if (!arguments.length) return width;
        width = _;
        svgWidth = width + margin.left + margin.right;
        return chart;
    };

    chart.height = function (_) {
        if (!arguments.length) return height;
        height = _;
        svgHeight = height + margin.top + margin.bottom;
        heightBlock = ((height - vpadBlockTop) / nBlocks) - vpadBlock;
        return chart;
    };

    chart.vpadBlockTop = function (_) {
        if (!arguments.length) return vpadBlockTop;
        vpadBlockTop = _;
        heightBlock = ((height - vpadBlockTop) / nBlocks) - vpadBlock;
        return chart;
    };

    chart.drawBorderTop = function (_) {
        if (!arguments.length) return drawBorderTop;
        drawBorderTop = _;
        return chart;
    }

    chart.drawBorderLeft = function (_) {
        if (!arguments.length) return drawBorderLeft;
        drawBorderLeft = _;
        return chart;
    }

    chart.drawBorderRight = function (_) {
        if (!arguments.length) return drawBorderRight;
        drawBorderRight = _;
        return chart;
    }

    chart.resizeSvg = function (_) {
        if (!arguments.length) return resizeSvg;
        resizeSvg = _;
        return chart;
    }

    chart.x = function (_) {
        if (!arguments.length) return x;
        x = _;
        return chart;
    }

    chart.zoomCallbacks = function (_) {
        if (!arguments.length) return zoomCallbacks;
        zoomCallbacks = _;
        return chart;
    }

    return chart;

}
