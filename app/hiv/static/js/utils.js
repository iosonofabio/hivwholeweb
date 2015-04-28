// support variables and functions for the d3 plots
var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹",
    formatPower = function(d) { 
        var tmp='';
        if (d<0){tmp+='⁻';}
        return tmp+(Math.abs(d) + "").split("").map(function(c) { return superscript[c]; }).join(""); };

