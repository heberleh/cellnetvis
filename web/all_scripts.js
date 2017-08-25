// Author Henry Heberle, 2015
// henry at icmc . usp . br

// This file is part of CellNetVis.

// CellNetVis is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// CellNetVis is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with CellNetVis. If not, see <http://www.gnu.org/licenses/>.

// Those that have employed the tool CellNetVis should mention:
// CellNetVis: a web tool for visualization of biological networks using force-directed layout constrained by cellular components
// Henry Heberle, Marcelo Falsarella Carazzolle, Guilherme P. Telles, Gabriela Vaz Meirelles, Rosane Minghim
// BMC Bioinformatics
// 2017


/* global cell, cytosol_keys,globalnetwork,globalkeepEnabled,unselectedEdge,
 * unselectedEdgeDarker,selectedEdge,globalFile,is_force_running,degree,force,
 * used_organelles,global_cc_data,allcc,w,h,radius,padding,delta,scale,
 * selected_attribute,linkDistance,linkStrength,selected_node,G,bundling,
 * bundling_data,is_force_runing,can_run_forces,showing_bundling,
 * main_components_ids,cytoplasm_id,extracellular_id,
 * plasma_membrane_or_wall_id*/

globalAlpha = 0;
cellids = {}; //organelas
cell = {};
cytosol_keys = {};
globalnetwork = null;
globalkeepEnabled = false;
unselectedEdge = 0.07;
unselectedEdgeDarker = 0.2;
selectedEdge = 0.5;
globalFile = null;
dragging = false;
is_force_running = true;
selected_cc_count = {};
cc_attribute_name = null;

force = null;

used_organelles = {};
global_cc_data = null;
allcc = [];
w = 850,
        h = 850,
        radius = 5,
        padding = 1, // separation between circles
        delta = 2 * radius,
        ddelta = 2 * radius,
        ddeltao = 3 * radius;
scale = d3.scale.category10();

selected_attribute = "Selected CC";
coloring = ['rgb(255,255,204)', 'rgb(255,237,160)', 'rgb(254,217,118)', 'rgb(254,178,76)', 'rgb(253,141,60)', 'rgb(252,78,42)', 'rgb(227,26,28)', 'rgb(189,0,38)', 'rgb(128,0,38)'];

cc_picker = null;

//var nucleus = {xmin: 605, xmax: 997, ymin: h - 27, ymax: h - 300};
//var plasma_membrane = {xmin: 1, xinfmax: 152, yinfmin: h - 15, yinfmax: h - 595, xsupmax: 945, ysupmin: h - 616, ysupmax: h - 700};
//var delta = 0;
linkDistance = 70;
linkStrength = 0.3;
charge = 0;
selected_node = null;
G = new jsnx.Graph();

bundling = null;
bundling_data = null;
is_force_runing = true;
can_run_forces = true;
showing_bundling = false;


transform_svg_separator = ",";
if ((/Edge\/12./i.test(navigator.userAgent)) || (/Edge\/13./i.test(navigator.userAgent))
        || (/Edge\/14./i.test(navigator.userAgent)) || (/Edge\/15./i.test(navigator.userAgent))
        || (/MSIE 10/i.test(navigator.userAgent)) || (/MSIE 11/i.test(navigator.userAgent))
        || (/MSIE 12/i.test(navigator.userAgent)) || (/MSIE 13/i.test(navigator.userAgent))
        || (/MSIE 9/i.test(navigator.userAgent)
                || /rv:11.0/i.test(navigator.userAgent))) {
    transform_svg_separator = " ";
}


/**
 * Hash for faster execution of tick() method that repositions hundreds of nodes
 */
main_components_ids = {"cytoplasm": 0, "extracellular": 1, "cell surface": 1, "plasma_membrane": 2, "cell_wall": 2};
cytoplasm_id = main_components_ids["cytoplasm"];
extracellular_id = main_components_ids["extracellular"];
plasma_membrane_or_wall_id = main_components_ids["plasma_membrane"];

/**
 * Returns the numeric index of a cellular component.
 * @param {string} name The name of a cellular component
 * @returns The numeric index of this component in the hash "main_components_ids" 
 */
function get_component_id(name) {
    if (name in main_components_ids) {
        return main_components_ids[name];
    } else {
        return 3;
    }
}

/**
 * Given a string, change the first char to UpperCase
 */
String.prototype.capitalize = function() {
    return this.charAt(0).toUpperCase() + this.slice(1);
}


/**
 * @description The global force-layout model
 */
force = d3.layout.force()
        .gravity(0.000001)
        .charge(charge)
        .linkDistance(linkDistance)
        .linkStrength(linkStrength)
        .size([w, h]);


/**
 * Resets the force-layout model.
 * @param {list} nodes The nodes the will replace the old ones.
 * @param {list} links The links that will replace the old ones.
 */
function reset_force(nodes, links) {
    if (force !== null) {
        force.stop();
        force.nodes(null);
        force.links(null);
        force.resume();
        force.stop();
        force.nodes(nodes);
        force.links(links);
        force.start();
    }
}

$(document).ready(function () {

    // when document is ready, load the basic cell diagram
    loadDiagram(null, null, false);


    d3.select("#donutchart_div").append("svg").attr("id","donutchart").style("width","300px").style("height","350px");

    // integration with IIS -> checks the ~file~ parameter in the URL
    integrate();

    $('#exact').click(function () {
        search('#search_text');
    });

    $('#search_text')[0].oninput = function () {
        search('#search_text');
    };

    $('#find').click(function () {
        search('#search_text')
    });
    $('#keep').click(keep);
    $('#clear').click(function(){clear_visibility(true)});


    $('#showlabels').click(showAllLabels);
    $('#hidelabels').click(hideAllLabels);
    $('#netmeas').click(compute_network_measures);


    $('#restart_forces').click(restartForces);
    $('#stop_forces').click(stopForces);

    $('#create_edge_bundling').click(function () {
        create_bundling(globalnetwork);
    });

    $('#show_edge_bundling').click(function () {
        draw_edge_bundling(globalnetwork);
    });
    $('#remove_bundling').click(remove_edge_bundling);
    $('#hide_bundling').click(hide_edge_bundling);

    $('#getpic').click(getPicture);

    $('#export_donut').click(exportDonutChart);

    document.getElementById('file').addEventListener('change', handleFileSelect, false);

    $('#alphaCollideSlider').on("change",changeAlpha);
    $('#chargeSlider').on("change",changeCharge);
    $('#strengthSlider').on("change",changeStrength);
    $('#frictionSlider').on("change",changeFriction);

    $('#change_label').on("change",changeLabels);

    $('#create_cerebral').click(create_cerebral_layout);
    d3.select('#cerebralDiv').style("visibility","hidden");

    $('input[type="checkbox"]').on('change', function() {
        $('input[name="' + this.name + '"]').not(this).prop('checked', false);
    });

    document
            .getElementById('cc_question')
            .style.display = 'none';

});


/**
 * Creates a list of numbers from start to count
 * @param {integer} start 
 * @param {integer} count 
 */
function range(start, count) {
    return Array.apply(0, Array(count))
            .map(function (element, index) {
                return index + start;
            });
}

/**
 * If force-layout is running stop it, otherwise start it.
 */
function play_stop_forces() {
    if (is_force_running) {
        is_force_running = false;
        force.stop();
    } else {
        is_force_running = true;
        force.start();
    }
}

/**
 * Update force-layout parameters according to the web interface.
 */
function changeForceParameters() {
    force.linkDistance(d3.select("#distance").node().value);
    force.friction(d3.select("#friction").node().value);
    force.linkStrength(d3.select("#strength").node().value);
    force.start();
}


// integration with IIS
function handleFileIIS(f) {
    $(document).ready(function () {
       // console.log("File: " + f);
       // console.log("handle File IIS");

        //console.log("handle File IIS - reset Global");
        //loadDiagram();
      //  console.log("handle File IIS - load diagram");
        loadDiagram(f, loadNetwork, true);
      //  console.log("handle File IIS - loadNetwork");
    });
}

/**
 * Handles the "file" parameter of the URL, that is, handles the XGMML network given through URL
 * @param {event} evt The event object
 */
function handleFileSelect(evt) {
    $(document).ready(function () {
        try {
            $("#loading-container").css("display", "block");
            //d3.select("#network_group").selectAll("*").remove();
            //d3.select("#chart").selectAll("*").remove();

            //loadDiagram();

            var files = evt.target.files; // FileList object
            // files is a FileList of File objects. List some properties.
            for (var i = 0, f; f = files[i]; i++) {
                loadDiagram(f, loadNetwork, false);
                //console.log("leu arquivo");
                if (i == 1) {
                    alert("Select only one XGMML file.");
                    break;
                }
            }
        } catch (err) {
            console.log(err.message);
        }
    });
}



/**
 * Computes the new position of each node so that they try to avoid collision.
 * @param {numeric} alpha The alpha parameter form 0 to 1
 * @author Adapted from http://bl.ocks.org/GerHobbelt/3116713
 */
function collide(alpha) {
    return function (d) {
        var data = null;

        if (globalkeepEnabled) {
            data = d3.selectAll(".node").filter(function (d) {
                return d.visible;
            }).data();
        }else{
            data = d3.selectAll('.node').data();
        }

        var rb = 2 * radius + padding;
        d3.geom.quadtree(data).visit(function (quad, x1, y1, x2, y2) {
            if (quad.point && (quad.point !== d)) {
                var x = d.x - quad.point.x,
                        y = d.y - quad.point.y,
                        l = Math.sqrt(x * x + y * y);
                if (l < rb) {
                    l = (l - rb) / l * alpha;
                    d.x -= x *= l;
                    d.y -= y *= l;
                    quad.point.x += x;
                    quad.point.y += y;
                }
            }
            return x1 > (d.x + rb) || x2 < (d.x - rb) || y1 > (d.y + rb) || y2 < (d.y - rb);
        });
    };
}



/**
 * Given a node of a network, change its selected cellular component, that is, change its position in the cell diagram.
 * @param {network} network 
 * @param {node} node The selected node.
 * @param {string} selectedCC The cellular component that will be set as Selected CC of this node.
 */
function changeSelectedCC(network, node, selectedCC){

        var old_component = node.attributes["Selected CC"].value.replace(/ /g, "_");

        var current_component = selectedCC.replace(/ /g, "_");

        node.attributes["Selected CC"].value = current_component.replace(/_/g, " ");

        if (current_component == "cell_surface"){
            current_component = "extracellular";
        }else{
            if (!(current_component in cell) || current_component === "" || current_component === null) {
                console.log(current_component+ "is not in cell, placing in cytoplasm instead.");
                current_component = "cytoplasm";
            }
        }

        node.group = node.attributes[d3.select("#color").node().value].value;
        node.component = current_component;
        node.component_id = get_component_id(current_component);
        node.attributes["Changed Selected CC"] =  {"type":"boolean", value: true};

        if (current_component in selected_cc_count) {
            selected_cc_count[current_component] = selected_cc_count[current_component] + 1;
        } else {
            selected_cc_count[current_component] = 1;
        }

        if (old_component in selected_cc_count) {
            selected_cc_count[old_component] = selected_cc_count[old_component] - 1;
        }

        update_valid_ccs_for_diagram(network);
        remove_unused_organelles(); // also shows new used organelles
        updateColor();
        updateDonutChart();
}


/**
 * Given a string that indicates a cellular component that is not default in CellNetVis tries to convert the value to a standard value to works in the system. The function tries to identify substrings, for instance, if a string contains the work "golgi" the function will return "golgi_apparatus".
 * @param {string} value The cellular component that is being verified.
 * @returns A cellular component (id) in the standards of CellNetVis
 */
function convertToCellNetVisCC(value){
    if (value == "unknown" || value == "" || value == null){
        return "not_found";
    }

    if (value in cell){
        return value;
    }

    for(var v in cell){
        if(cell.hasOwnProperty(v)){
            if(v.indexOf(value) >= 0 || value.indexOf(v) >= 0){
                return v;
            }
        }
    }

    if( value.indexOf("cell_surface">=0) ){
       return "extracellular";
    }
    if (value.indexOf("endoplas")>=0){
        return "endoplasmic_reticulum";
    }
    if (value.indexOf("golgi")>=0){
        return "golgi_apparatus";
    }
    if (value.indexOf("mitocho")>=0){
        return "mitochondrion";
    }
    if (value.indexOf("nucle")>=0){
        return "nucleus";
    }
    if (value.indexOf("glycos")>=0){
        return "glycosome";
    }
    if (value.indexOf("glyox")>=0){
        return "glyoxysome";
    }
    if (value.indexOf("endoso")>=0){
        return "endosome";
    }
    if (value.indexOf("plast")>=0){
        return "plastid";
    }
    if (value.indexOf("lyso")>=0){
        return "lysosome";
    }
    if (value.indexOf("perox")>=0){
        return "peroxisome";
    }
    if (value.indexOf("tubule")>=0){
        return "microtubule_organising_centre";
    }
    if (value.indexOf("centro")>=0){
        return "centrosome";
    }
    if (value.indexOf("amylo")>=0){
        return "amyloplast";
    }
    if (value.indexOf("apico")>=0){
        return "apicoplast";
    }
    if (value.indexOf("chloro")>=0){
        return "chloroplast";
    }
    if (value.indexOf("ext")>=0){
        return "extracellular";
    }
    if (value.indexOf("vac")>=0){
        return "vacuole";
    }
    if (value.indexOf("membran")>=0){
        return "plasma_membrane";
    }
    if (value.indexOf("wall")>=0){
        return "cell_wall";
    }
    return "not_found";
}



function updateContentTable(network, node) {
    $("#rows").empty();    

    var html = "<tr><td class=\"text-left\">" + "ID" + "</td>";
    html += "<td class=\"text-left\">" + node.id.replace("n_","") + "</td>";
    html += "</tr>";
    html += "<tr><td class=\"text-left\">" + "Selected CC" + "</td>";
    html += "<td class=\'text-left\'>" +
                "<select id=\'updateCCchooser\'>"
                        +"<option value=\'" + node.attributes["Selected CC"].value + "\'>"
                        +                     node.attributes["Selected CC"].value.replace(/_/g, " ") + "*"
                        + "</option>" ;
    for (var cc in cell){
        if (cc != node.attributes["Selected CC"].value){
            html += "<option value=\'"+cc+"\'>"+cc.replace(/_/g, " ")+"</option>";
        }
    }
    html += "</td>";
    html += "</tr>";



    html += "<tr><td class=\"text-left\">" + "Extracted CCs" + "</td>";
    html += "<td class=\"text-left\">" + node.attributes["CellNetVis CCs"].value.toString().replace(/_/g, " ").replace(/,/g, ", ") + "</td>";
    html += "</tr>";

    $('#rows').append(html);



    $(document).ready(function(){

        $('#updateCCchooser').val(node.attributes["Selected CC"].value);

        $('#updateCCchooser').change(function(){
            changeSelectedCC(network, node, this.value);
        });

    });

    for (name in node.attributes) {
        if (name.toLowerCase() == "selected cc" || name.toLowerCase() == "cellnetvis ccs"){
            continue;
        }
        var html = "<tr><td class=\"text-left\">" + name + "</td>";
        if (name.toLowerCase() === "uniprot id" || name.toLowerCase() === "uniprotid" || name.toLowerCase() === "uniprot_id") {
            html += "<td class=\"text-left\">" + "<a  href='http://www.uniprot.org/uniprot/" + node.attributes[name].value + "' target='_blank' title='Click opens a new tab.'>" + node.attributes[name].value + " (click to search it at UniProt.org)</a> " + "</td>";
        } else {
            html += "<td class=\"text-left\">" + node.attributes[name].value + "</td>";
        }
        html += "</tr>";
        $('#rows').append(html);
    }
}


function update_valid_ccs_for_diagram(network){
    for (cc in used_organelles){
        if (used_organelles.hasOwnProperty(cc)) {
            used_organelles[cc] = false;
        }
    }

    used_organelles["plasma_membrane"] = true;
    used_organelles["extracellular"] = true;
    used_organelles["cytoplasm"] = true;
    used_organelles["cytosol"] = true;
    used_organelles["nucleus"] = true;

    for (var index = 0; index < network.nodes.length; index++){
        if (used_organelles[network.nodes[index].attributes["Selected CC"].value.replace(/ /g, "_")] == false){
            used_organelles[network.nodes[index].attributes["Selected CC"].value.replace(/ /g, "_")] = true;
        }
    }

    restartForces();
}


function deselectNode(d) {

    if (d !== null) {
        d.transition()
                .duration(200)
                .attr("r", radius);
    }
}

function selectNode(d) {
    if (d !== null) {
        d.transition()
                .duration(200)
                .attr("r", radius * 1.3);
    }
}


function restartForces() {
    if (can_run_forces) {
        if (bundling === null) {
            remove_edge_bundling();
        }
        is_force_running = true;
        d3.selectAll(".node").each(function (d) {
            d.fixed = false;
        });
        force.start();
    }
}

function stopForces() {
    force.stop();
    d3.selectAll(".node").each(function (d) {
        d.fixed = true;
    });
    is_force_running = false;
    //force.resume();
}


function constraint_network(network) {
//console.log("iniciou a função constraing_network");

//var nucleus = {xmin: 605, xmax: 997, ymin: h - 27, ymax: h - 300};
//var plasma_membrane = {xmin: 1, xinfmax: 152, yinfmin: h - 15, yinfmax: h - 595, xsupmax: 945, ysupmin: h - 616, ysupmax: h - 700};
    var linkDistance = 45;
    var svg = d3.select("#network_group");

    var link = svg.append("g").classed('edges', true)
            .selectAll("line")
            .data(network.links)
            .enter().append("svg:line")
            .attr("class", "edge")
            .attr("id", function (d) {
                return (d.id);
            });

    d3.selectAll("line").attr('stroke', 'rgba(32,32,32,' + unselectedEdgeDarker + ')');
//    var node = svg.selectAll("circle")
//            .data(network.nodes)
//            .enter().append("svg:circle")
//            .attr("r", r - .75)
//            .attr("class", "node")
//            .style("fill", function(d) {
//                return z(d.group);
//            })
//            .style("stroke", function(d) {
//                return d3.rgb(z(d.group)).darker();
//            })
//            .call(force.drag);




    //========================================== adapted from: http://bl.ocks.org/norrs/2883411    norrs’s block #2883411 June 6, 2012
    var node_drag = d3.behavior.drag()
            .on("dragstart", dragstart)
            .on("drag", dragmove)
            .on("dragend", dragend);

    function dragstart(d, i) {
        if (!showing_bundling && bundling === null) {
            d.fixed = true;
            dragging = true;
        }
    }

    function dragmove(d, i) {
        if (!showing_bundling && bundling === null) {

            d.x += d3.event.dx;
            d.y += d3.event.dy;
            d.px = d.x;
            d.py = d.y;

            tick();
            // d3.select("#"+d.id+"_g").attr("transform_svg_separatorm", function (d) {
            //     return 'translate(' + d.x + transform_svg_separator + d.y + ')';
            // });

            // d3.selectAll("link").filter(function(l) {
            //     return l.source === d;
            // }).attr("x1", d.x).attr("y1", d.y);

            // d3.selectAll("link").filter(function(l) {
            //     return l.target === d;
            // }).attr("x2", d.x).attr("y2", d.y);
        }
        // this is the key to make it work together with updating both px,py,x,y on d !
    }

    function dragend(d, i) {
        if (!showing_bundling && bundling === null) {
            d.fixed = true;
            tick();
            force.resume();
            dragging = false;
        }
    }

//=============================================

    function highlight_possible_ccs(d){
        remove_unused_organelles();

        // reset class 'highglight cc'
        d3.selectAll(".group_component").classed("highlight_cc", false);
               
        // highlight current selected cc of the selected node
        d3.select("#" + convertToCellNetVisCC(d.attributes["Selected CC"].value.toLowerCase().replace(/ /g, "_")) + "_group").style("visibility","visible").classed("highlight_cc", true);

        var localizations = d.attributes["CellNetVis CCs"].value; //this special attribute is a list

        // highlight other CCs
        cell_wall = false;
        for(var loc in localizations){
            var cc = localizations[loc];
            if (cc == "cytosol"){
                cc = "cytoplasm";
            }
            if(cc == "cell_wall"){
                cell_wall=true;
            }
            if (cc != "not_found"){
                d3.select("#" + cc + "_group").style("visibility","visible").classed("highlight_cc", true);
            }
        }
                
        if (is_force_running){
            //stopForces();
            //restartForces();
            d.fixed = false;
        }else{
            if (cell_wall){
                d.fixed = false;
            }
        }
    }

    var gnodes = svg.selectAll('g.gnode')
            .data(network.nodes)
            .enter()
            .append('g')
            .attr("id",function(d){return d.id+"_g"})
            .classed('gnode', true);

    var node = gnodes.append("svg:circle")
            .attr("draggable", true)
            .attr("r", radius - .75)
            .attr("class", "node")
            .attr("id", function (d) {
                return d.id;
            })
            .style("fill", function (d) {
                return scale(d.group);
            })
            .style("stroke", function (d) {
                return d3.rgb(scale(d.group)).darker();
            })
            .on("mouseover", function (d)
            {
                d3.select("#label" + d.id).style("visibility", "visible");
            })
            .on("mouseout", function (d) {
                if (d.showlabel === false) {
                    d3.select("#label" + d.id).style("visibility", "hidden");
                }
            })
            .on("click", function (d) {
                hideAllLabels();

                var thisselected = d.selected;
                d3.selectAll(".node")
                        .filter(function (d) {
                            return d.visible;
                        })
                        .style("stroke-width", "1px")
                        .style("stroke", function (dd) {
                            return d3.rgb(scale(dd.group)).darker();
                        })
                        .style("opacity", 0.2);
                var nds = d3.selectAll(".node").filter(function (d) {
                    return d.visible;
                }).data();
                for (nd in nds) {
                    nds[nd].selected = false;
                }


                if (showing_bundling) {
                    d3.selectAll(".bundling").style('stroke', 'rgb(32,32,32)').style("stroke-opacity", 0.03);
                } else {
                    d3.selectAll(".edge").attr('stroke', 'rgba(32,32,32,' + unselectedEdge + ')');
                }



                if (thisselected === false) {
                    //seleciona nó e vizinhos
                    d.selected = true;

                    highlight_possible_ccs(d);

                    thisnode = d3.select("#" + d.id);
                    thisnode.style("stroke-width", "2px").style('stroke', '#ce6f20');
                    thisnode.style("opacity", 1);
                    d3.select("#label" + d.id).style("visibility", "visible");
                    d.showlabel = true;
                    //network.hash_ids.indexOf($(this).attr("source")

                    //console.log(thisnode.data()[0].attributes);

                    for (var i = 0; i < network.links.length; i++) {
                        if (network.links[i].source.id === d.id) {
                            var node = d3.select("#" + network.links[i].target.id);
                            if (node.data()[0].visible) {
                                node.data()[0].selected = true;
                                node.style("stroke-width", "2px");
                                node.style("opacity", 1);
                                d3.select("#label" + node.data()[0].id).style("visibility", "visible");
                                node.data()[0].showlabel = true;
                                if (showing_bundling) {
                                    d3.selectAll("." + network.links[i].id).style('stroke', 'rgb(32,32,32)').style("stroke-opacity", 0.85);
                                } else {
                                    d3.select("#" + network.links[i].id).attr('stroke', 'rgba(32,32,32,' + selectedEdge + ')');
                                }
                            }

                        } else {
                            if (network.links[i].target.id === d.id) {
                                var node = d3.select("#" + network.links[i].source.id);
                                if (node.data()[0].visible) {
                                    node.data()[0].selected = true;
                                    node.style("stroke-width", "2px");
                                    node.style("opacity", 1);
                                    d3.select("#edge_" + network.links[i].source.id + '_' + network.links[i].target.id).attr('stroke', 'rgba(32,32,32,' + selectedEdge + ')');
                                    d3.select("#label" + node.data()[0].id).style("visibility", "visible");
                                    node.data()[0].showlabel = true;
                                    if (showing_bundling) {
                                        d3.selectAll("." + network.links[i].id).style('stroke', 'rgb(32,32,32)').style("stroke-opacity", 0.85);
                                    } else {
                                        d3.select("#" + network.links[i].id).attr('stroke', 'rgba(32,32,32,' + selectedEdge + ')');
                                    }
                                }
                            }
                        }
                    }
                } else { // if node is already selected


                    if (!globalkeepEnabled) {
                        clear_visibility(false);
                        unselect_nodes();

                        remove_unused_organelles();
                        d3.selectAll(".group_component").classed("highlight_cc", false);

                    } else {
                        d3.selectAll(".node").filter(function (d) {
                            return d.visible;
                        }).style("opacity", 1).style("stroke-width", "1px").style("stroke", function (dd) {
                            return d3.rgb(scale(dd.group)).darker();
                        });
                        if (showing_bundling) {
                            d3.selectAll(".bundling").style('stroke', 'rgb(32,32,32)').style("stroke-opacity", unselectedEdgeDarker);
                        } else {
                            d3.selectAll(".edge").attr('stroke', 'rgba(32,32,32,' + unselectedEdgeDarker + ')');
                        }
                    }

                    if (is_force_running){
                        force.start();
                    }
                    // }else{
                    //     stopForces();
                    // }

                }
//                if (d.showlabel == true) {
//                    d3.select("#label" + d.id).style("visibility", "hidden");
//                    d.showlabel = false;
//                } else {
//                    d3.select("#label" + d.id).style("visibility", "visible");
//                    d.showlabel = true;
//                }
                updateContentTable(network, d);
                //console.log("update table");
            })
            .call(node_drag);



    var labels = gnodes.append("text")
            .text(function (d) {
                return d.label;
            })
            .attr("id", function (d) {
                return "label" + d.id;
            })
            .style("visibility", "hidden")
            .attr('pointer-events', 'none')
            .attr("font-family", "arial")
            .attr("font-weight", "bold")
            .attr("stroke", "white")
            .attr("stroke-width", ".1px")
            .attr("font-size", "11px")
            .attr("class", "node_label");

    force
            .nodes(network.nodes)
            .links(network.links)
            .on("tick", tick)
            .start();

    globalnetwork = network;


    // map each node, computing new x,y values 
    var map = {
        0: function(d, cell, cell_component){
            return mapping_cytoplasm(d.x, d.y, cell, cell_component, Math.sqrt((d.x - cell_component.cx) * (d.x - cell_component.cx) + (d.y - cell_component.cy) * (d.y - cell_component.cy)));
        },

        1: function(d, cell, cell_component){
            return mapping_extracellular(d.x, d.y, cell);
        },

        2: function(d, cell, cell_component){
            return mapping_walls(d.x, d.y, cell[d.component], Math.sqrt((d.x - cell_component.cx) * (d.x - cell_component.cx) + (d.y - cell_component.cy) * (d.y - cell_component.cy)));
        },

        3: function(d, cell, cell_component){
            return mapping_organelles(d.x, d.y, cell_component, Math.sqrt((d.x - cell_component.cx) * (d.x - cell_component.cx) + (d.y - cell_component.cy) * (d.y - cell_component.cy)));
        }
    };


    //  this function is called for each force iteraction
    function tick() {

        if (globalkeepEnabled) {
            d3.selectAll(".node").filter(function (d) {
                return d.visible;
            }).each(function (d) {
                var result;
                var cell_component = cell[d.component];

                result = map[d.component_id](d, cell, cell_component);
                d.x = result.x;
                d.y = result.y;
            });
        } else {
            d3.selectAll(".node").each(function (d) {
                var result;
                var cell_component = cell[d.component];

                result = map[d.component_id](d, cell, cell_component);
                d.x = result.x;
                d.y = result.y;
            });
        }


        if (globalAlpha != 0 && is_force_running && !dragging) {
            d3.selectAll(".node").each(collide(globalAlpha));
        }

        d3.selectAll(".gnode").attr("transform", function (d) {
            if (!isNaN(d.x) && !isNaN(d.y)){
                return 'translate(' + d.x + transform_svg_separator + d.y + ')';
            }else{
                d.x = d.px;
                d.y = d.py;
            }
        });

        d3.selectAll(".edge").attr("x1", function (d) {
                    return d.source.x;
                })
                .attr("y1", function (d) {
                    return d.source.y;
                })
                .attr("x2", function (d) {
                    return d.target.x;
                })
                .attr("y2", function (d) {
                    return d.target.y;
                });
        }
}



function compute_network_measures() {
    if (globalnetwork !== null) {

        $.blockUI({
            message: '<h2>Computing network measures...</h2>',
            onBlock: function(){
                var G = new jsnx.Graph();

                G.addNodesFrom(range(0, globalnetwork.nodes.length));

                degree = [];
                for (var i = 0; i < globalnetwork.nodes.length; i++){
                    degree.push(0);
                }

                for (var i = 0; i < globalnetwork.links.length; i++){
                    G.addEdge(globalnetwork.links[i].source.index, globalnetwork.links[i].target.index);
                    degree[globalnetwork.links[i].source.index] += 1;
                    degree[globalnetwork.links[i].target.index] += 1;
                }

                //var degree = jsnx.degree_centrality(G);
                //var closeness = jsnx.closenessCentrality(G);
                var betweenness = jsnx.betweennessCentrality(G);
                //var eigenvector = jsnx.eigenvectorCentrality(G);
                var clustering = jsnx.clustering(G);
                //betweenness._numberValues[0];

                //adding more atributes
                for (var i = 0; i < globalnetwork.nodes.length; i++) {
                    var node = globalnetwork.nodes[i];
                    node.attributes["[betweenness]"] = {type: "numeric", value: betweenness._numberValues[i]};
                    node.attributes["[clustering]"] = {type: "numeric", value: clustering._numberValues[i]};
                    node.attributes["[degree]"] = {type: "numeric", value: degree[i]};
                }

                $(document).ready(function () {
                    $('#color').append("<option value=\"[betweenness]\">[betweenness]</option>");
                });

                $(document).ready(function () {
                    $('#change_label').append("<option value=\"[betweenness]\">[betweenness]</option>");
                });

                //change_label
                //$(document).ready(function() {
                //    $('#color').append("<option value=\"[eigenvector]\">[eigenvector]</option>");
                //});
                $(document).ready(function () {
                    $('#color').append("<option value=\"[clustering]\">[clustering]</option>");
                });
                $(document).ready(function () {
                    $('#color').append("<option value=\"[degree]\">[degree]</option>");
                });

                $(document).ready(function () {
                    $('#change_label').append("<option value=\"[clustering]\">[clustering]</option>");
                });
                $(document).ready(function () {
                    $('#change_label').append("<option value=\"[degree]\">[degree]</option>");
                });
                $.unblockUI();
            }
        });
    }
}

String.prototype.capitalize = function () {
    return this.replace(/(^|\s)([a-z])/g, function (m, p1, p2) {
        return p1 + p2.toUpperCase();
    });
};

function loadNetwork(f, iis) {
    //console.log("iniciou loadnetwork");
    alreadyHasANetwork = true;
    error = false;
    if ((typeof iis === 'undefined') || (!iis)) {
        var reader = new FileReader();
        reader.onload = function (e) {
            parseXGMML(e.target.result);
            $("#loading-container").css("display", "none");
        };
        if (!error) {
            reader.readAsText(f);
        }
    } else {
        parseXGMML(f);
    }
}



function parseXGMML(xml) {
    $("#cy2").empty();
    d3.select("#cerebralDiv").style("visibility","hidden");


    alreadyHasANetwork = true;
    var network = {};
    network.hash_ids = {};
    var error = false;
    var missingCC = false;

    //console.log("iniciou parser");
    network.nodes = [];
    //Leitura do XGMML
    xmlDoc = $.parseXML(xml),
            $xml = $(xmlDoc),
            $title = $xml.find("title");
    //console.log("leu o xml");
    //console.log($title);

    var att_names = [];

    var att_list_names = [];

    var node_index = 0;
    $(xml).find("node").each(function ()
    {

        //console.log($(this).attr("ID"));
        var node_id =  "n_"+$(this).attr("ID").replace(/\./g, "_").replace(/\//g, "_");

        var node_att_names = [];


        if (node_id in network.hash_ids) {
            alert("XGMML redundant. There is more than one node with ID=" + node_id);
        } else {

            network.hash_ids[node_id] = node_index;
            node_index = node_index + 1;


            var node_label = $(this).attr("label");
            if (node_label == null || node_label == ""){
                //node_label = $(this).attr("name");
                //if (node_label == null || node_label == ""){
                    node_label = node_id.replace("n_","");
                //}
            }

            var node_attributes = {};

            $(this).find("att").each(function () {
                var at_name = $(this).attr("name");

                var low_att = at_name.toLowerCase();
                if (low_att === "selected_cc" || low_att === "selected cc"||
                    low_att === "localization" ) {
                    at_name = "Selected CC";
                }

                if ($(this).attr("type").toLowerCase() !== 'list') {
                    if (at_name === "Selected CC") { //ignores if there are more than one value and overwrite the last
                        var att = {type: $(this).attr("type"), value: $(this).attr("value").toLowerCase()};
                        node_attributes[at_name] = att;
                    } else {
                        if (node_att_names.indexOf(at_name) < 0) {
                            var att = {type: $(this).attr("type"), value: $(this).attr("value")};
                            node_attributes[at_name] = att;
                            node_att_names.push(at_name);
                        } else { //if already exist we are going to make a list of values
                            var l = [].concat(node_attributes[at_name].value).concat($(this).attr("value"));
                            l = $.grep(l, function (n) {
                                return(n);
                            });
                            //console.log(l);
                            var att = {type: $(this).attr("type"), value: l};
                            node_attributes[at_name] = att;
                        }
                    }
                    if (att_names.indexOf(at_name) < 0) {
                        att_names.push(at_name);
                        if (att_list_names.indexOf(at_name) < 0) {
                            $(document).ready(function () {
                                    $('#color').append("<option value=\"" + at_name + "\">" + at_name + "</option>");
                                    $('#change_label').append("<option value=\"" + at_name + "\">" + at_name + "</option>");
                            });
                        }
                    }
                } else {
                    if (at_name !== "Selected CC") {
                        att_list_names.push(at_name);
                    }
                }

            });

            if (!missingCC) {
                //console.log("Armazenamento dos atributos.");
                // console.log(node_attributes);

                if (att_names.indexOf("Selected CC") < 0) {
                    missingCC = true;
                }else{
                    var current_component = node_attributes["Selected CC"].value.replace(/ /g, "_");
                    if (current_component === "cell_surface"){
                        current_component = "extracellular";
                    }else{
                        if (!(current_component in cell) || current_component === "") {
                            current_component = "cytoplasm";
                        }
                    }

                    var node = {
                        id: node_id,
                        name: node_label,
                        label: node_label,
                        showlabel: false,
                        selected: false,
                        sx: 0,
                        sy: 0,
                        visible: true,
                        attributes: node_attributes,
                        group: node_attributes[d3.select("#color").node().value].value,
                        component: current_component,
                        component_id: get_component_id(current_component)};
                    //this component_id is just a index of the main type of CC... only to
                    //compute the TICK first step-> decide if it is in cytoplasm,
                    //cell wall or plasma mebrane or organelles

                    network.nodes.push(node);

                    //console.log("Nó adicionado à rede.");

                    // read Selected CC without replacing unknown or not defined CC by Cytoplasm
                    var component = node_attributes["Selected CC"].value.replace(/ /g, "_").toLowerCase();
                    if (component in selected_cc_count) {
                        selected_cc_count[component] = selected_cc_count[component] + 1;
                    } else {
                        selected_cc_count[component] = 1;
                    }
                    used_organelles[component] = true;
                }
            }
            if (missingCC){
                node_attributes["Changed Selected CC"] = {"type":"boolean", "value":false};

                var node = {
                    id: node_id,
                    name: node_label,
                    label: node_label,
                    showlabel: false,
                    selected: false,
                    sx: 0,
                    sy: 0,
                    visible: true,
                    attributes: node_attributes,
                    group: null,
                    component: null,
                    component_id: null};
                //this component_id is just a index of the main type of CC... only to
                //compute the TICK first step-> decide if it is in cytoplasm,
                //cell wall or plasma mebrane or organelles
                network.nodes.push(node);
                //console.log("Nó adicionado à rede.");
            }
        }
        resetForces();
    });


    // get list of UNIPROT/etc ID
    var query_ids = [];
    var auxnode = network.nodes[0];
    var attnames = Object.keys(auxnode.attributes);
    var query_type = ""; // Ensembl, UniProt, Entrez or InnateDB
    var attname_query = "";
    if (missingCC){
        for (var i = 0; i < attnames.length; i++){
            attname_query = attnames[i];
            if (attnames[i].toLowerCase().indexOf("uniprot") !== -1){
                query_type = "UniProt";
                break;
            }
            if (attnames[i].toLowerCase().indexOf("ensembl") !== -1){
                query_type = "Ensembl";
                break;
            }
            if (attnames[i].toLowerCase().indexOf("entrez") !== -1){
                query_type = "Entrez";
                break;
            }
            if (attnames[i].toLowerCase().indexOf("innatedb") !== -1){
                query_type = "InnateDB";
                break;
            }
        }
    }


    var showPopupCC = function(event) {
        $('#selectCC').empty();

        $.each(att_names, function (index, value) {
            $('#selectCC').append($('<option/>', {
                value: value,
                text : value
            }));
        });


        d3.select("#cc_question").style("display","block");
        // document
        //     .getElementById('cc_question')
        //     .style.display = 'block';
    };


    $('#submit_cc').click(function() {
        //console.log("picker...");

        cc_picker = "Innate DB";

        // if (document.getElementById("check_InnateDB").checked){
        //     cc_picker = "Innate DB";
        // }else{

        //     if(document.getElementById("check_FirstCC").checked){
        //         cc_picker = "First CC";
        //         cc_attribute_name = $('#selectCC :selected').text();
        //     }else{
        //         if(document.getElementById("check_AllInCytosol").checked){
        //             cc_picker = "All in Cytosol";
        //         }else{
        //             if(document.getElementById("check_AllInExtracellular").checked){
        //                 cc_picker = "All in Extracellular";
        //             }else{
        //                 console.log("error. cc_picker == null");
        //             }
        //         }
        //     }
        // }
        // console.log("cc_picker: "+cc_picker);

        // document
        //     .getElementById('cc_question')
        //     .style.display = 'none';

        d3.select("#cc_question").style("display","none");
        complete_network();
    });

    if (missingCC){

        $('#color').append("<option value=\"" + "Changed Selected CC" + "\">" + "Changed Selected CC" + "</option>");
        $('#change_label').append("<option value=\"" + "Changed Selected CC" + "\">" + "Changed Selected CC" + "</option>");

        document
        .getElementById('popup-click')
        .addEventListener('click', showPopupCC);

         document
        .getElementById('popup-click').click();


    }else{
        //console.log("ignoring the dialog.")
        complete_network();
    }

    function create_selected_all_in(component){

        for (var index = 0; index < network.nodes.length; index++){

            var selectedCC = {'type':"string", "value": component};
            network.nodes[index].attributes["Selected CC"] = selectedCC;

            var current_component = selectedCC.value;

            network.nodes[index].group = network.nodes[index].attributes[d3.select("#color").node().value].value;
            network.nodes[index].component = current_component;
            network.nodes[index].component_id = get_component_id(current_component);

            if (current_component in selected_cc_count) {
                selected_cc_count[current_component] = selected_cc_count[current_component] + 1;
            } else {
                selected_cc_count[current_component] = 1;
            }
            used_organelles[current_component] = true;
        }
        $('#color').append("<option value=\"Selected CC\">Selected CC</option>");
        $('#change_label').append("<option value=\"Selected CC\">Selected CC</option>");
    }


    function create_selected_cc_using_first_cc(){

        var aux = network.nodes[0].attributes[cc_attribute_name].value;
        var comma = aux.split(",");
        var dotcomma = aux.split(";");
        var separator = ";";
        if (comma.length > dotcomma.length){
            separator = ",";
        }

        for (var index = 0; index < network.nodes.length; index++){

            var selectedCC = {'type':"string", "value":"cytoplasm"}; // DEFAULT VALUE

            //console.log(cc_attribute_name);
            //console.log(network.nodes[index].attributes[cc_attribute_name].value);

            var node_ccs = network.nodes[index].attributes[cc_attribute_name].value.split(separator);

            var new_node_ccs = [];

            for (var i = 0; i < node_ccs.length; i++){
                var value = node_ccs[i].replace(/ /g, "_");
                if ( value in cell || value == "cytosol" || value == "cell_surface"){
                    new_node_ccs.push(node_ccs[i]);
                }else{
                   new_node_ccs.push(convertToCellNetVisCC(value));
                }
            }

            if (new_node_ccs.length > 0){
                selectedCC = {'type':"string", "value": new_node_ccs[0]};
            }

            var ccs_string = "";
            if(new_node_ccs.length > 1){
                for (var c = 0; c < new_node_ccs.length; c++){
                    ccs_string +=  new_node_ccs[c]+ ", "
                }
            }
            ccs_string+=new_node_ccs[new_node_ccs.length-1];

            network.nodes[index].attributes["CellNetVis CCs"] = {"type":"string","value": ccs_string};

            //console.log(selectedCC.value);
            network.nodes[index].attributes["Selected CC"] = selectedCC;

            var current_component = selectedCC.value.replace(/ /g, "_");

            if (current_component == "cell_surface"){
                current_component = "extracellular";
            }else{
                if (!(current_component in cell) || current_component === "" || current_component === null) {
                    //console.log("not in cell");
                    //console.log(current_component);
                    current_component = "cytoplasm";
                }
            }

            network.nodes[index].group = network.nodes[index].attributes[d3.select("#color").node().value].value;
            network.nodes[index].component = current_component;
            network.nodes[index].component_id = get_component_id(current_component);

            current_component = network.nodes[index].attributes["Selected CC"].value.replace(/ /g, "_");

            if (current_component in selected_cc_count) {
                selected_cc_count[current_component] = selected_cc_count[current_component] + 1;
            } else {
                selected_cc_count[current_component] = 1;
            }
            used_organelles[current_component] = true;

        }
        $('#color').append("<option value=\"Selected CC\">Selected CC</option>");
        $('#change_label').append("<option value=\"Selected CC\">Selected CC</option>");
    }

    function create_selected_cc_using_innatedb(){

        console.log("Quering Localization from InnateDB.");

        var genes = [];
        //console.log("att to query");
        //console.log(attname_query);

        if (attname_query.toLowerCase() === "innatedb reference"){
            for (var i = 0; i < network.nodes.length; i++){
                genes.push(network.nodes[i].attributes[attname_query].value.split("?id=")[1].toString());
            }
        }else{
            for (var i = 0; i < network.nodes.length; i++){
                genes.push(network.nodes[i].attributes[attname_query].value.toString());
            }
        }

        var query = 'xref='+query_type+'&ids=' + JSON.stringify(genes).replace('\[', '').replace('\]', '').replace(/\"/g, '');

        console.log(query);

        $.ajax({
            dataType: "json",
            async: false,
            type: "GET",
            url: "http://www.innatedb.com/cerebralLocalizationWS.do",
            data: query,
            success: function(data){
                for (var index = 0; index < genes.length; index++){

                    network.nodes[index].attributes["Selected CC"] = {'type':"string", "value":data[genes[index]].toLowerCase()};
    // console.log(index + "," +data[genes[index]]+ ", "+genes[index]+"  -  "+network.nodes[index].attributes["Localization"].value);

                    var current_component = network.nodes[index].attributes["Selected CC"].value.replace(/ /g, "_");

                    if (current_component == "cell_surface"){
                        current_component = "extracellular";
                    }else{
                        if (!(current_component in cell) || current_component === "" || current_component === null) {
                            //console.log("not in cell");
                            //console.log(current_component);
                            current_component = "cytoplasm";
                        }
                    }

                    network.nodes[index].group = network.nodes[index].attributes[d3.select("#color").node().value].value;
                    network.nodes[index].component = current_component;
                    network.nodes[index].component_id = get_component_id(current_component);

                    current_component = network.nodes[index].attributes["Selected CC"].value.replace(/ /g, "_");

                    if (current_component in selected_cc_count) {
                        selected_cc_count[current_component] = selected_cc_count[current_component] + 1;
                    } else {
                        selected_cc_count[current_component] = 1;
                    }
                    used_organelles[current_component] = true;
                }
            },
            error: function(jqXHR, textStatus, errorThrown)
            {
                alert("Error in innateDB web service: http://www.innatedb.com/cerebralLocalizationWS.do. It is not possible to query Selected CC from this server.");
                console.log(errorThrown);
            }
        });

        $('#color').append("<option value=\"Selected CC\">Selected CC</option>");
        $('#change_label').append("<option value=\"Selected CC\">Selected CC</option>");
    }


    function complete_network(){
        //console.log("completing network...");

        console.log("Missing cc: " +missingCC);

        if (missingCC){
            // // 1. use the first value from attribute chosed by user
            // if(cc_picker == "First CC"){
            //     create_selected_cc_using_first_cc();
            // }else{
            //     // 4. query from InnateDB
            //     if(cc_picker == "InnateDB" && query_type != ""){
            //         create_selected_cc_using_innatedb();
            //     }else{
            //         if (cc_picker == "All in Cytosol"){
            //             create_selected_all_in("cytoplasm");
            //         }else{
            //             if (cc_picker == "All in Extracellular"){
            //                 create_selected_all_in("extracellular");
            //             }
            //         }
            //     }
            // }
            if(query_type != ""){
                create_selected_cc_using_innatedb();
            }else{
                alert("We could not query the Localization attribute from InnateDB since your network does not contain any of the attributes: UniProt, Ensembl, Entrez and InnateDB. Please consider reading the Help page for more information.");
                resetGlobalVariables();
                return;
            }
            console.log("Selected CC attribute created.")
        }



        var hash_names = {};
        for (name in network.nodes[0].attributes){
            if (network.nodes[0].attributes.hasOwnProperty(name)) {
                hash_names[name] = 1;
            }
        }

        var cc_cellnetvis_name = "not_found";
        var possible_names = {"cellnetvis_ccs":1, "cellular_component":2, "cellular_component_(go)":3, "go_localizations":4};

        for(name in hash_names){
            if(hash_names.hasOwnProperty(name)){
                var aux = name.toLowerCase().replace(/ /g, "_");
                if (aux in possible_names){
                    cc_cellnetvis_name = name;
                    break;
                }else{
                    if(name.indexOf("components") >= 0 || name.indexOf("localizations") >= 0 ){
                        cc_cellnetvis_name = name;
                        break;
                    }
                }
            }
        }

        if (cc_cellnetvis_name == "not_found"){
            for (var index = 0; index < network.nodes.length; index++){
                network.nodes[index].attributes["CellNetVis CCs"] = {"type":"list","value": ["not_found"]};
            }
        }else{

            var cc_attribute_name = cc_cellnetvis_name;
            //console.log(cc_attribute_name);
            //console.log(network.nodes[0].attributes);
            var aux = network.nodes[0].attributes[cc_attribute_name].value;
            var comma = aux.split(",");
            var dotcomma = aux.split(";");
            var separator = ";";
            if (comma.length > dotcomma.length){
                separator = ",";
            }

            for (var index = 0; index < network.nodes.length; index++){

                var node_ccs = network.nodes[index].attributes[cc_attribute_name].value.split(separator);

                var new_node_ccs = [];

                for (var i = 0; i < node_ccs.length; i++){
                    var value = node_ccs[i].toLowerCase().replace(/ /g, "_");
                    if (value == "cytosol" || value == "cell_surface"){
                        new_node_ccs.push(node_ccs[i]);
                    }else{
                        var new_cc = convertToCellNetVisCC(value);
                        if (new_cc != "not_found" && new_node_ccs.indexOf(new_cc)<0){
                            new_node_ccs.push(new_cc);
                        }
                    }
                }

                network.nodes[index].attributes["CellNetVis CCs"] = {"type":"list","value": new_node_ccs};
            }
        }

        if (!error) {
            used_organelles["plasma_membrane"] = true;
            used_organelles["extracellular"] = true;
            used_organelles["cytoplasm"] = true;
            used_organelles["cytosol"] = true;
            used_organelles["nucleus"] = true;

            remove_unused_organelles();


            var not_cytosol = {"cytoplasm": "", "plasma_membrane": "", "extracellular": "", "cytosol": "", "cell_wall": "", "cell_surface":""};
            for (key in cell) {
                if (!(key in not_cytosol)) {
                    cytosol_keys[key] = "";
                }
            }

            network.links = [];

            $(xml).find("edge").each(function ()
            {
                var source_id = "n_"+$(this).attr("source").replace(/\./g, "_").replace(/\//g, "_");
                var target_id = "n_"+$(this).attr("target").replace(/\./g, "_").replace(/\//g, "_");

                if (source_id !== target_id) { //remove self-edge
                    var sourceIdx = network.hash_ids[source_id];
                    var targetIdx = network.hash_ids[target_id];
                    var edge = {
                        label: $(this).attr("label"),
                        id: "edge_" + source_id + "_" + target_id,
                        source: sourceIdx,
                        target: targetIdx,
                        value: 1
                    };
                    //console.log(edge);

                    network.links.push(edge);


                }
            });

            constraint_network(network);


            createDonutChart();

            //console.log("donut criado.");
            updateNumbers(network);

            d3.select('#color').property('value', 'Selected CC');
            d3.select('#change_label').property('value', 'Selected CC');
            updateColor();


            if (network.nodes.length > 1000 || network.links.length > 10000){
                $("#alphaCollide").slider('disable');
            }else{
                $("#alphaCollide").slider('enable');
            }
        }
        network = null;
        resetForces();
    }
}


function updateNodesGroup(attribute) {
    //console.log(attribute);
    var possible_values = getAttributeValues(attribute);
    //console.log(possible_values);

    var is_numeric = isNumericVector(possible_values);
    //console.log(is_numeric);

    if (is_numeric) {
        var max = Math.max.apply(Math, possible_values);
        var min = Math.min.apply(Math, possible_values);
        if (min === max) {
            max = max + 1;
        }
        scale = d3.scale.quantize().domain([min, max]).range(coloring);
    } else {
        if (possible_values.length < 11) {
            if (attribute.toLowerCase() === "selected_cc" || attribute.toLowerCase() === "selected cc") {
                scale = d3.scale.category20();//.domain(allcc);
            } else {
                scale = d3.scale.category10().domain(possible_values);
            }
        } else {
            if (attribute.toLowerCase() === "selected_cc" || attribute.toLowerCase() === "selected cc") {
                scale = d3.scale.category20();//.domain(allcc);
            } else {
                scale = d3.scale.category20().domain(possible_values);
            }

        }
    }

//atualiza a variável group de cada nó
    d3.selectAll('.node').attr('group', function (d) {
        d.group = d.attributes[attribute].value;
        return d.group;
    });
    d3.selectAll('.node').style("fill", function (d) {
        return scale(d.group);
    }).style("stroke", function (d) {
        return d3.rgb(scale(d.group)).darker();
    });
}

function isNumericVector(vector) {
    for (var i = 0; i < vector.length; i++) {
        if ($.isNumeric(vector[i]) === false) {
            return  false;
        }
    }
    return true;
}

function getAttributeValues(attribute) {
    var list = d3.selectAll('.node').data();
    var values = [];
    //console.log("attribute value: " + attribute);
    for (var i = 0; i < list.length; i++) {
        if (values.indexOf(list[i].attributes[attribute].value) < 0) {
            values.push(list[i].attributes[attribute].value);
        }
    }
//console.log("get values attributes");
//console.log(values);
    return values;
}



function updateColor() {
    var attribute = d3.select("#color").node().value;
    updateNodesGroup(attribute);
//    updateD3Scale(scale);
}

function showAllLabels() {
//d3.selectAll(".node_label").style("visibility", "visible");
    var list = d3.selectAll('.node').filter(function (d) {
        return d.visible;
    }).data();
    for (var i = 0; i < list.length; i++) {
        list[i].showlabel = true;
        d3.select("#label" + list[i].id).style("visibility", "visible");
    }
}

function hideAllLabels() {
    d3.selectAll(".node_label").style("visibility", "hidden");
    var list = d3.selectAll('.node').data();
    for (var i = 0; i < list.length; i++) {
        list[i].showlabel = false;
    }
}

function setLabels(attr) {
//rascunho da função
    var nds = d3.selectAll(".node").data();
    for (nd in nds) {
        nds[nd].label = nds[nd].attributes[attr].value; //verificar se é isso mesmo
    }
}

function updateDiagramLabels() {
    var nds = d3.selectAll(".node").data();
    for (nd in nds) {
        d3.select("#label" + nds[nd].id).text(nds[nd].label);
    }
}

function remove_unused_organelles() {
    //used_organelles
    for (var key in used_organelles) {
        if (used_organelles.hasOwnProperty(key)) {
            if (used_organelles[key] === false) {
                d3.select("#" + key + "_group").classed("removed", true).style("visibility","hidden");
            }else{
                d3.select("#" + key + "_group").classed("removed", false).style("visibility","visible");
            }
        }
    }
}

function updateNumbers(network) {
    d3.select("#Nnodes").text(network.nodes.length + " nodes");
    d3.select("#Nedges").text(network.links.length + " links");
}



function loadDiagram(file, func, iis) {
    d3.xml("circular_cell_structure.svg", "image/svg+xml", function (xml) {
        var importedNode = document.importNode(xml.documentElement, true);
        $(document).ready(function () {
            resetGlobalVariables();
            d3.select("#chart").selectAll("*").remove();
            //console.log("removeu diagrama");
            d3.select("#chart").node().appendChild(importedNode);
            //console.log("adicionou novo diagrama");
            var cellt = d3.selectAll("#cell_group")[0][0].getAttribute("transform");
           // console.log("cell transform:")
            //console.log(cellt);

            cellt = cellt.replace("translate(", "");
            cellt = cellt.replace(")", "");
            cellt = cellt.split(transform_svg_separator);
            var cell_t = [parseFloat(cellt[0]), parseFloat(cellt[1])];

            function updateCXCY() {
                var group_components = d3.selectAll(".group_component")[0];
                for (var i = 0; i < group_components.length; i++) { //para cada grupo
                    var group_component = group_components[i];
                    var t = group_component.getAttribute("transform");  //seleciona a partição celular deste grupo

                    t = t.replace("translate(", "");
                    t = t.replace(")", "");

                    t = t.split(transform_svg_separator);

                    var component_name = group_component.id.replace("_group", "");
                    var component = d3.select("#" + component_name)[0][0];

                    var cx = parseFloat(component.getAttribute("cx")) + parseFloat(t[0]) + cell_t[0];
                    var cy = parseFloat(component.getAttribute("cy")) + parseFloat(t[1]) + cell_t[1];
                    var rmax = parseFloat(component.getAttribute("r"));
                    var id = component.getAttribute("id");
                    cellids[id] = i;

                    if (used_organelles[id] !== true) {
                        used_organelles[id] = false;
                    }

                    cell[id] = {cx: cx, cy: cy, rmax: rmax, rmin: 0, id: id};


                }

                //set rmin for plasma_membrane = rmax from cytoplasm
                cell["plasma_membrane"].rmin = cell["cytoplasm"].rmax + 0.6 * delta;
                cell["cell_wall"].rmin = cell["plasma_membrane"].rmax + 0.6 * delta;
            }

            //d3.select("#cell_group").attr("transform", "translate(" + 0 + "," + 0 + ")");
            var svg = d3.select("#diagram").append("svg:g").attr("id", "network_group").attr("transform", "translate(" + 0 + transform_svg_separator + 0 + ")");

            function move() {
                if (bundling === null) {
                    var circle = d3.select("#" + this.id.replace("_group", "")).node();
                    var data = cell[this.id.replace("_group", "")];
                    //var t = d3.transform(d3.select(this));
                    var dx = d3.event.x - circle.getAttribute("cx");
                    var dy = d3.event.y - circle.getAttribute("cy");

                    d3.select(this).attr("transform", "translate(" + dx + transform_svg_separator + dy + ")");

                    updateCXCY();
                    force.start();
                }
            }
            //console.log("Definiu função move");

            var not_move = {"cytoplasm_group": "", "plasma_membrane_group": "", "cell_wall_group": "", "extracellular_group": ""};
            var group_components = d3.selectAll(".group_component")[0];
            //console.log(group_components);
            for (var i = 0; i < group_components.length; i++) {
                var group_component = group_components[i];
                if (!(group_component.id in not_move)) {
                    d3.select("#" + group_component.id).call(d3.behavior.drag()
                            .on("dragstart", function () {
                                d3.event.sourceEvent.stopPropagation();
                            })
                            .on("drag", move))
                            .on("dragend", function () {
                                force.stop();
                                force.start();
                            });
                }
            }
            updateCXCY();
        });
        if (func !== null && file !== null) {
            func(file, iis);
        }

        $('#extracellular').click(unselect_nodes);
        $('#cytoplasm').click(unselect_nodes);
        $('#plasma_membrane').click(unselect_nodes);
        $('#cell_wall').click(unselect_nodes);
    });

}






function mapping_walls(x, y, cell_component, R) {
//NÃO SE APLICA AO CYTOPLASM   ou EXTRACELLULAR
    if (R > cell_component.rmax - ddelta) { //se estiver fora do circulo, pela borda mais externa, seta o raio para essa borda
        return {x: (cell_component.rmax - ddelta) / R * (x - cell_component.cx) + cell_component.cx, y: (cell_component.rmax - ddelta) / R * (y - cell_component.cy) + cell_component.cy};
    } else {
        if (R < cell_component.rmin + ddelta) {
            return {x: (cell_component.rmin + ddelta) / R * (x - cell_component.cx) + cell_component.cx, y: (cell_component.rmin + ddelta) / R * (y - cell_component.cy) + cell_component.cy};
        } else {
            return {x: x, y: y};
        }
    }
}

//retorna o correto valor de x,y após aplicar as restrições
function mapping_organelles(x, y, cell_component, R) {
//NÃO SE APLICA AO CYTOPLASM   ou EXTRACELLULAR

    if (R > cell_component.rmax - ddeltao) { //se estiver fora do circulo, pela borda mais externa, seta o raio para essa borda
        return {x: (cell_component.rmax - ddeltao) / R * (x - cell_component.cx) + cell_component.cx, y: (cell_component.rmax - ddeltao) / R * (y - cell_component.cy) + cell_component.cy};
    } else {
        return {x: x, y: y};
    }
}


function mapping_cytoplasm(x, y, cell, cell_component, R) {
    //se estiver fora do cytoplasm
    if (R > cell_component.rmax - delta) { //se estiver fora do circulo, pela borda mais externa, seta o raio para essa borda
        return {x: (cell_component.rmax - ddelta) / R * (x - cell_component.cx) + cell_component.cx, y: (cell_component.rmax - ddelta) / R * (y - cell_component.cy) + cell_component.cy};
    } else {
    //se é pra estar no citosol mas está em uma organela, tira dessa organela:

        for (key in cytosol_keys) {
            cell_component = cell[key];
            if (used_organelles[key] === true) {
                R = Math.sqrt((x - cell_component.cx) * (x - cell_component.cx) + (y - cell_component.cy) * (y - cell_component.cy));
                if (R < (cell_component.rmax + ddelta)) {
                    return {x: (cell_component.rmax + ddelta) / R * (x - cell_component.cx) + cell_component.cx, y: (cell_component.rmax + ddelta) / R * (y - cell_component.cy) + cell_component.cy};
                }
            }
        }
    }
    // se não, se não tiver de transladar o vértice...
    return {x: x, y: y};
}


function mapping_extracellular(x, y, cell) {
    var border = null;
    if (used_organelles["cell_wall"] === true) {
        border = cell["cell_wall"];
    } else {
        border = cell["plasma_membrane"];
    }

    var R = Math.sqrt((x - border.cx) * (x - border.cx) + (y - border.cy) * (y - border.cy));

    if (R < (border.rmax + ddelta)) { //se estiver dentro do circulo, pela borda mais externa, seta o raio para essa borda
        return {x: (border.rmax + ddelta) / R * (x - border.cx) + border.cx, y: (border.rmax + ddelta) / R * (y - border.cy) + border.cy};

    } else {
        //se estiver fora do retângulo (bordas máximas)
        var xt = 0;
        var yt = 0;
        if (x < (ddelta)) {
            xt = ddelta;
        } else {
            if (x > (w - ddelta)) {
                xt = w - ddelta;
            } else {
                xt = x;
            }
        }
        if (y < (ddelta)) {
            yt = ddelta;
        } else {
            if (y > (h - ddelta)) {
                yt = h - ddelta;
            } else {
                yt = y;
            }
        }
        return {x: xt, y: yt};
    }
}



function updateDonutChart(){
    $("#chart2 svg").empty();
    createDonutChart();
}



function createDonutChart() {
    var dataset = [];
    var label = "";
    for (var key in selected_cc_count) {
        if (selected_cc_count.hasOwnProperty(key) && selected_cc_count[key] > 0) {
            if (key === "") {
                dataset.push({"label": "*undefined", "value": selected_cc_count[key]});
            } else {
                if (key.indexOf("endoplasmic_reticulum") > -1) {
                    label = "ER";
                } else {
                    if (key.indexOf("microtubule_organising_centre") > -1) {
                        label = "MTOC";
                    } else {
                        label = key.replace(/_/g, " ").toLowerCase().capitalize();
                    }
                }
                dataset.push({"label": label, "value": selected_cc_count[key]});
            }
        }
    }

    global_cc_data = dataset;
    allcc = [];
    global_cc_data.forEach(function (d) {
        allcc.push(d.label);
    });

    nv.addGraph(function () {
        var chart = nv.models.pieChart()
                .x(function (d) {
                    return d.label;
                })
                .y(function (d) {
                    return d.value;
                })
                .showLabels(true)     //Display pie labels
                .labelThreshold(.05)  //Configure the minimum slice size for labels to show up
                .labelType("percent") //Configure what type of data to show in the label. Can be "key", "value" or "percent"
                .donut(true)          //Turn on Donut mode. Makes pie chart look tasty!
                .donutRatio(0.35)     //Configure how big you want the donut hole size to be.
                ;

        d3.select("#chart2 svg")
                .datum(dataset)
                .transition().duration(350)
                .call(chart);

        //var t = d3.select("#chart2 svg").select("g").attr("transform");
        //t = t.replace("translate(", "");
        //t = t.replace(")", "");
        //t = t.split(",");
        //var trans = "translate("+(parseInt(t[0])-65)+","+t[1]+")";
        //d3.select("#chart2 svg").select("g").attr("transform", trans);

        var total = 0;
        dataset.forEach(function (d) {
            total = total + d.value;
        });
        var tp = function (key, y, e, graph) {
            return '<h3 style="background-color: '
                    + e.color + '">' + key + '</h3>'
                    + '<p>' + parseInt(y) + '</p>'
                    + '<p>' + (y * 100 / total).toFixed(2) + '%</p>';
        };

        chart.tooltipContent(tp);
        return chart;
    });
}


function create_bundling(network) {
    if (is_force_running) {
        alert("Can't compute edge bundling while Force algorithm is running. Please click on [Stop] button from Force section. After that, try computing the bundles again. After computing the bundling, you need to click on [Show]. To [Start] the force algorithm again you need to click on [Destroy] edge bundling.");
        return null;
    }

    alert("Edge bundling will be computed and may take some seconds or minutes. When it ends, a message will pop-up.");
    var fbundling = null;

    d3.selectAll(".node").each(function (d) {
        can_run_forces = false;
        d.draggable = false;
    });

    var nnodes = {};
    var nodes = d3.selectAll(".node").each(function (d) {
        nnodes[d.index] = {x: d.x, y: d.y};
    });

    var edges = network.links;
    var eedges = [];
    for (var i = 0; i < edges.length; i++) {
        if (edges[i].source.visible && edges[i].target.visible)
            eedges.push({source: edges[i].source.index, target: edges[i].target.index, id: edges[i].id});
    }

    $.blockUI({
        //message: '<h2>Computing edge bundling...</h2>',
        onBlock: function(){
            //Run the FDEB algorithm using default values on the data
            fbundling = d3.ForceEdgeBundling().nodes(nnodes).edges(eedges);
            if(fbundling != null) {$.unblockUI()};
        },
        onUnblock: function(){
            bundling_data = fbundling;
            bundling = bundling_data();
            $("#hide_bundling").prop("disabled", false);
            $("#show_edge_bundling").prop("disabled", false);
            $("#remove_bundling").prop("disabled", false);
            $("#restart_forces").prop("disabled", true);
            alert("Edge bundling computed. Click [Show] on [Edge Bundling] menu to see the result.");
        }
    });
}


function draw_edge_bundling(network) {
    if (bundling === null) {
        $("#show_edge_bundling").prop("disabled", true);
        return;
    }

    //set normal edges to hidden
    d3.selectAll(".edge").attr("visibility", "hidden");

    //create bundling edges with specific class: bundling
    var d3line = d3.svg.line()
            .x(function (d) {
                return d.x;
            })
            .y(function (d) {
                return d.y;
            })
            .interpolate("linear");

    var svg = d3.select(".edges");
//    bundling.forEach(function(edge_subpoint_data) {
//        // for each of the arrays in the results
//        // draw a line between the subdivions points for that edge
//
//        svg.append("path")
//                .attr("d", d3line(edge_subpoint_data))
//                .style("stroke-width", 1)
//                .style("stroke", "#ff2222")
//                .style("fill", "none")
//                .style('stroke-opacity', 0.15)
//                .attr("class", "bundling " + edge_id); //use opacity as blending
//    });

    d3.selectAll(".bundling").remove();
    var edges = bundling_data.edges();
    for (var i = 0; i < bundling.length; i++) {
        svg.append("path")
                .attr("d", d3line(bundling[i]))
                .style("stroke-width", 1)
                .style('stroke', 'rgb(32,32,32)')
                .style("stroke-opacity", unselectedEdgeDarker)
                .style("fill", "none")
                .attr("class", "bundling " + edges[i].id); //use opacity as blending
    }

    d3.selectAll(".edge").each(function (d) {
        d3.select("." + d.id).attr("visibility", d.source.visible && d.target.visible ? "visible" : "hidden");
    });

    showing_bundling = true;

}


function remove_edge_bundling() {
    if (bundling === null) {
        return;
    }

    //set normal edges to visible and hide bundling
    hide_edge_bundling();

    //remove all lines bundling class
    d3.selectAll(".bundling").remove();

    //set button show bundling to disabled
    $("#hide_bundling").prop("disabled", true);
    $("#show_edge_bundling").prop("disabled", true);
    $("#remove_bundling").prop("disabled", true);
    $("#restart_forces").prop("disabled", false);

    can_run_forces = true;

    d3.selectAll(".node").each(function (d) {
        d.draggable = true;
    });

    bundling = null;


}

function hide_edge_bundling() {
    if (bundling === null) {
        $("#hide_bundling").prop("disabled", true);
        return;
    }
    //set edge bundling to hidden
    d3.selectAll(".bundling").attr("visibility", "hidden");

    //set normal edges to visible
    d3.selectAll(".edge").attr("visibility", function (d) {
        return d.source.visible && d.target.visible ? "visible" : "hidden";
    });

    showing_bundling = false;
}


// IIS Integration
function getParams() {

    var params = {},
            pairs = document.URL.split('?')
            .pop()
            .split('&');

    for (var i = 0, p; i < pairs.length; i++) {
        p = pairs[i].split('=');
        params[ p[0] ] = p[1];
    }
    //console.log("Parâmetros: na URL.");
    //console.log(params);
    return params;
}




function integrate() {
    params = getParams();

    var file = "";
    if ("file" in params) {
        file = unescape(params["file"]);
        $("#loading-container").css("display", "block");
        $.ajax({
            url: file,
            type: 'GET',
            cache: false,
            dataType: 'text',
            timeout: 100000,
            error: function () {
                alert('Error loading XGMML document');
            },
            success: function (xml) {
                //$("div#debug").text(xml);
                handleFileIIS(xml);
            },
            error: function (xml) {
                alert("Error loading file.");
            },
            complete: function () {
                $("#loading-container").css("display", "none");
            }
        });
    }
}




function resetGlobalVariables() {
    cell = {};
    used_organelles = {};
    global_cc_data = null;
    allcc = [];
    selected_attribute = "node.fillColor";
    selected_node = null;

    if (globalnetwork !== null) {
        if (globalnetwork.nodes !== null)
            globalnetwork.nodes = null;
        if (globalnetwork.links !== null)
            globalnetwork.links = null;
        if (globalnetwork.hash_ids !== null)
            globalnetwork.hash_ids = null;
    }
    globalnetwork = null;
    globalkeepEnabled = false;
    cellids = {}; //organelas
    cytosol_keys = {};
    globalFile = null;
    is_force_running = true;
    $('#color').empty();

    if (force !== null) {
        force.stop();
        force.nodes(null);
        force.links(null);
        force = null;
    }


    linkDistance = 45;
    force = d3.layout.force()
            .gravity(0.000001)
            .charge(charge)
            .linkDistance(linkDistance)
            .linkStrength(linkStrength)
            .size([w, h]);
    resetForces();

    $("#hide_bundling").prop("disabled", true);
    $("#show_edge_bundling").prop("disabled", true);
    $("#remove_bundling").prop("disabled", true);
    $("#restart_forces").prop("disabled", false);

    is_force_running = true;
    bundling = null;
    showing_bundling = false;
    bundling_data = null;
    can_run_forces = true;

    selected_cc_count = {};

    $('#cy2').empty();
}






function getPicture() {
    var serializer = new XMLSerializer();
    var svg = document.getElementById("diagram");
    if (svg !== null) {
        var chooser = document.getElementById("figureformat");
        var choosen_type = chooser.options[chooser.selectedIndex].value;
        var serialized = serializer.serializeToString(svg);
        var new_blob = new Blob([serialized], {type: "image/svg+xml;charset=" + document.characterSet});

        var prefix_document_name = document.getElementById("svg-name").value;
        var document_name = null;
        if (prefix_document_name === null || prefix_document_name === "") {
            document_name = "diagram" + choosen_type;
        } else {
            document_name = prefix_document_name + choosen_type;
        }

        if (choosen_type === ".svg") {
            //alert("Saving diagram as "+document_name+". It may take some seconds.");
            saveAs(new_blob, document_name);
            return false;
        } else {
            if (choosen_type === ".png") {
                var image = new Image();
                image.src = 'data:image/svg+xml;base64,' + window.btoa(serialized);
                image.onload = function () {
                    var scale = 5;
                    var canvas = document.createElement('canvas');
                    canvas.width = image.width * scale;
                    canvas.height = image.height * scale;

                    var context = canvas.getContext('2d');
                    context.scale(scale, scale);
                    context.drawImage(image, 0, 0);
                    // set to draw behind current content
                    context.globalCompositeOperation = "destination-over";

                    // set background color
                    context.fillStyle = '#FFFFFF'; // <- background color

                    // draw background / rect on entire canvas
                    context.fillRect(0, 0, canvas.width, canvas.height);

                    var a = document.createElement('a');
                    a.download = document_name;
                    a.href = canvas.toDataURL('image/png');
                    document.body.appendChild(a);
                    a.click();
                };
                return false;
            } else {
                var file = "";
                if (choosen_type === '.txt') { //exportar listas de elementos das intersecções em vez de diagrama propriamente dito
                    for (var i = 0; i < labelsDiagram.length; i++) {
                        var id = labelsDiagram[i].toUpperCase();
                        var set = intersectionsSet[id];

                        if (set !== null) {
                            var textName = "[";
                            textName = textName + document.getElementById("name" + id[0]).value;
                            for (var j = 1; j < id.length; j++) {
                                textName = textName + "] and [" + document.getElementById("name" + id[j]).value;
                            }
                            textName = textName + "]";

                            file = file + textName + ": ";
                            file = file + set.toString();
                            file = file + "\n";
                            //console.log(file);
                        }

                    }
                    saveAs(new Blob([file], {type: "text/plain;charset=" + document.characterSet}), document_name);
                    return false;
                }
            }
        }
    }
}



function exportDonutChart() {

    var serializer = new XMLSerializer();

    var fixer = document.querySelector("#donutchart");
    fixer.removeAttribute("xmlns");
    fixer.removeAttribute("xmlns:xlink");

    prefix = "http://www.w3.org/2000/xmlns/";
    if (!fixer.hasAttributeNS(prefix, "xmlns")) {
        fixer.setAttributeNS(prefix, "xmlns", "http://www.w3.org/2000/svg");
    }

    if (!fixer.hasAttributeNS(prefix, "xmlns:xlink")) {
        fixer.setAttributeNS(prefix, "xmlns:xlink", "http://www.w3.org/1999/xlink");
    }

    //var svg = document.getElementById("donutchart");
    var svg = fixer;
    var svg_aux = d3.select("#donutchart");
    var width = d3.select("#donutchart").style("width");
    var height = d3.select("#donutchart").style("height");
    svg_aux.attr('width', width);
    svg_aux.attr('height', height);
    svg_aux.style("font-size","12px").style("font-family","Arial");


    if (svg !== null) {

        var chooser = document.getElementById("donut_figureformat");
        var choosen_type = chooser.options[chooser.selectedIndex].value;

        d3.select("#donutchart").attr("width",300);
        d3.select("#donutchart").attr("height",350);
        if (choosen_type === ".png") {
            d3.select("#donutchart").style("width",null).style("height",null);
        }

        var serialized = serializer.serializeToString(svg);
        var new_blob = new Blob([serialized], {type: "image/svg+xml;charset=" + document.characterSet});
        var url = window.URL.createObjectURL(new_blob);

        var prefix_document_name = document.getElementById("donut-name").value;
        var document_name = null;
        if (prefix_document_name === null || prefix_document_name === "") {
            document_name = "cc-distribution" + choosen_type;
        } else {
            document_name = prefix_document_name + choosen_type;
        }

        if (choosen_type === ".svg") {
            //alert("Saving diagram as "+document_name+". It may take some seconds.");
            saveAs(new_blob, document_name);
            return false;
        } else {
            if (choosen_type === ".png") {
                var image = new Image();
                image.src = 'data:image/svg+xml;base64,' + window.btoa(serialized);
                image.onload = function () {
                    var scale = 5;
                    var canvas = document.createElement('canvas');
                    canvas.width = image.width * scale;
                    canvas.height = image.height * scale;

                    var context = canvas.getContext('2d');
                    context.scale(scale, scale);
                    context.drawImage(image, 0, 0);
                    // set to draw behind current content
                    context.globalCompositeOperation = "destination-over";

                    // set background color
                    context.fillStyle = '#FFFFFF'; // <- background color

                    // draw background / rect on entire canvas
                    context.fillRect(0, 0, canvas.width, canvas.height);

                    var a = document.createElement('a');
                    a.download = document_name;
                    a.href = canvas.toDataURL('image/png');
                    document.body.appendChild(a);
                    a.click();
                    //d3.select("#donutchart").style("width",width).style("height",height);
                };
                return false;
            } else { // if .csv
                return getPercents(global_cc_data, document_name);
            }
        }
    }
}




function getPercents(dataset, document_name) {
    if (dataset !== null) {

        var total = 0;
        dataset.forEach(function (d) {
            total = total + d.value;
        });

        str = "cellular component,count,percent\n";
        dataset.forEach(function (d) {
            str = str + d.label + "," + d.value + "," + (d.value / total).toFixed(4) + "\n";
        });
        saveAs(new Blob([str], {type: "text/plain;charset=" + document.characterSet}), document_name);
    }
    return false;
}


function search(inputID) {
    if (globalnetwork !== null) {
        clear_visibility(false);

        var exact = d3.select("#exact").node().checked;
        //ID do campo onde o texto de busca se encontra
        var str = d3.select(inputID).node().value.toUpperCase();


        var nds = d3.selectAll(".node").data();

//    hideAllLabels();

        if (showing_bundling) {
            d3.selectAll(".bundling").style('stroke', 'rgb(32,32,32)').style("stroke-opacity", unselectedEdge);
        } else {
            d3.selectAll(".edge").attr('stroke', 'rgba(32,32,32,' + unselectedEdge + ')');
        }

        d3.selectAll(".node")
////            .filter(function(d) {
////                return d.visible;
////            })
//            .style("stroke-width", "1px")
//            .style("stroke", function(dd) {
//                return d3.rgb(scale(dd.group)).darker();
//            })
                .style("opacity", 0.2);

//    for (nd in nds) {
//        nds[nd].selected = false;
//    }

        var splited_str = str.split("\n");
        var selected_nodes_ids = [];
        for (nd in nds) {
            for (var i in splited_str) {
                if (splited_str[i] === null) {
                    break;
                } else {
                    if (exact) {
                        if (nds[nd].label.toUpperCase().valueOf() === splited_str[i].toUpperCase().valueOf()) {
                            nds[nd].selected = true;
                            nds[nd].showlabel = true;
                            thisnode = d3.select("#" + nds[nd].id);
                            thisnode.style("stroke-width", "2px").style("opacity", 1);
                            d3.select("#label" + nds[nd].id).style("visibility", "visible");
                            selected_nodes_ids.push(nds[nd].id);
                            break;
                        }
                    } else {
                        if (nds[nd].label.toUpperCase().search(splited_str[i]) > -1) {
                            nds[nd].selected = true;
                            nds[nd].showlabel = true;
                            thisnode = d3.select("#" + nds[nd].id);
                            thisnode.style("stroke-width", "2px").style("opacity", 1);
                            d3.select("#label" + nds[nd].id).style("visibility", "visible");
                            selected_nodes_ids.push(nds[nd].id);
                            break;
                        }
                    }
                }
            }
        }
        selected_nodes_dic = {};
        for (var k = 0; k < selected_nodes_ids.length; k++) {
            selected_nodes_dic[selected_nodes_ids[k]] = "";
        }

        for (var i in globalnetwork.links) {
            if (globalnetwork.links[i].source.id in selected_nodes_dic && globalnetwork.links[i].target.id in selected_nodes_dic) {
                if (showing_bundling) {
                    d3.select("." + globalnetwork.links[i].id).style('stroke', 'rgb(32,32,32)').style("stroke-opacity", selectedEdge);
                } else {
                    d3.select("#" + globalnetwork.links[i].id).attr('stroke', 'rgba(32,32,32,' + selectedEdge + ')');
                }
            }
        }
    }
}


function updateGraphVisibility() {

    d3.selectAll(".node").style("visibility", function (d) {
        return d.visible ? "visible" : "hidden";
    });

    if (showing_bundling) {
        d3.selectAll(".edge").each(function (d) {
            d3.select("." + d.id).attr("visibility", d.source.visible && d.target.visible ? "visible" : "hidden");
        });
    } else {
        d3.selectAll(".edge").attr("visibility", function (d) {
            return d.source.visible && d.target.visible ? "visible" : "hidden";
        });
    }
}

//tornar os vértices não selecionados invisíveis (e suas arestas)
function keep() {
    globalkeepEnabled = true;
    var nds = d3.selectAll(".node").data();

    for (nd in nds) {
        if (nds[nd].selected) {
            nds[nd].visible = true;
        } else {
            nds[nd].visible = false;
        }
    }

    updateGraphVisibility();
    hideAllLabels();
    showAllLabels();
    force_only_selected();
}

function clear_visibility(force_start) {
    if (globalnetwork !== null) {
        globalkeepEnabled = false;

        var nds = d3.selectAll(".node").data();

        for (nd in nds) {
            nds[nd].visible = true;
            nds[nd].selected = false;
        }

        d3.selectAll(".node").style("opacity", 1).style("stroke-width", "1px").style("stroke", function (dd) {
            return d3.rgb(scale(dd.group)).darker();
        });

        if (showing_bundling) {
            d3.selectAll(".bundling").style('stroke', 'rgb(32,32,32)').style("stroke-opacity", 0.03);
        } else {
            d3.selectAll(".edge").attr('stroke', 'rgba(32,32,32,' + unselectedEdgeDarker + ')');
        }

        if (globalnetwork.nodes.length > 1000 || globalnetwork.links.length > 10000){
            $("#alphaCollide").slider('disable');
        }else{
            $("#alphaCollide").slider('enable');
        }

        if (force_start){
            force
                .nodes(globalnetwork.nodes)
                .links(globalnetwork.links)
                .start();
            restartForces();
        }

        // remove highlight of components
        d3.selectAll(".group_component").classed("highlight_cc", false);

        hideAllLabels();
        updateGraphVisibility();
    }
}


function unselect_nodes() {
    if (globalnetwork !== null) {

        var nds = d3.selectAll(".node").data();

        for (nd in nds) {
            if(nds[nd].visible == true){
                nds[nd].selected = false;
            }           
        }

        d3.selectAll(".node").style("opacity", 1).style("stroke-width", "1px").style("stroke", function (dd) {
            return d3.rgb(scale(dd.group)).darker();
        });

        if (showing_bundling) {
            d3.selectAll(".bundling").style('stroke', 'rgb(32,32,32)').style("stroke-opacity", 0.03);
        } else {
            d3.selectAll(".edge").attr('stroke', 'rgba(32,32,32,' + unselectedEdgeDarker + ')');
        }

        if (globalnetwork.nodes.length > 1000 || globalnetwork.links.length > 10000){
            $("#alphaCollide").slider('disable');
        }else{
            $("#alphaCollide").slider('enable');
        }

        // if (force_start){
        //     force
        //         .nodes(globalnetwork.nodes)
        //         .links(globalnetwork.links)
        //         .start();
        // }

        // remove highlight of components
        d3.selectAll(".group_component").classed("highlight_cc", false);

        hideAllLabels();
        updateGraphVisibility();
    }
}

function force_only_selected() {
    if (globalnetwork !== null) {
        newnodes = [];
        for (var i in globalnetwork.nodes) {
            if (globalnetwork.nodes[i].selected) {
                newnodes.push(globalnetwork.nodes[i]);
            }
        }

        var newlinks = [];
        for (var j in globalnetwork.links) {
            var edge = globalnetwork.links[j];
            if (edge.source.selected && edge.target.selected) {
                newlinks.push(edge);
            }
        }

        if (globalnetwork.nodes > 1500 && globalnetwork.nodes < 2000){
            if (newnodes > 100){
                $("#alphaCollide").slider('disable');
            }else{
                $("#alphaCollide").slider('enable');
            }
        }else{
            if (globalnetwork.nodes > 1000){
                if (newnodes > 200){
                    $("#alphaCollide").slider('disable');
                }else{
                     $("#alphaCollide").slider('enable');
                }
            }else{
                $("#alphaCollide").slider('enable');
            }
        }

        reset_force(newnodes, newlinks);
    }
}

function changeAlpha(){
    globalAlpha = $("#alphaCollide").slider("getValue");
    force.start();
}


function changeCharge(){
    force.charge($("#charge").slider("getValue"));
    force.start();
}

function changeStrength(){
    var tstrength = $("#strength").slider("getValue");
    force.linkStrength(tstrength);
    force.start();
}

function changeFriction(){
    force.friction($("#friction").slider("getValue"));
    force.start();
}


function resetForces(){

    if (globalnetwork!==null && globalnetwork.nodes.length > 300){
        $("#alphaCollide").slider('setValue',0);
        d3.select("#alphaLabel").text("0.0");
        changeAlpha();
    }else{
        $("#alphaCollide").slider('setValue',0.6);
        d3.select("#alphaLabel").text("0.6");
        changeAlpha();
    }

    if (globalnetwork!==null && globalnetwork.nodes.length > 900){
        $("#charge").slider('setValue', -7);
        d3.select("#chargeLabel").text("-7");
        changeCharge();

        $("#friction").slider('setValue',0.9);
        d3.select("#frictionLabel").text("0.9");
        changeFriction();
    }else{
        $("#charge").slider('setValue',-63);
        d3.select("#chargeLabel").text("-63");
        changeCharge();

        $("#friction").slider('setValue',0.7);
        d3.select("#frictionLabel").text("0.7");
        changeFriction();
    }

    $("#strength").slider('setValue',0.3);
    d3.select("#strengthLabel").text("0.3");
    changeStrength();



}

// function removeLinks(){
//     force.links([]);
//     d3.selectAll(".edge").style("visibility", "hidden");
//     force.start();
// }

// function restoreLinks(){
//     force.links(globalnetwork.links);
//     d3.selectAll(".edge").style("visibility", "visible");
//     force.start();
// }





function create_cerebral_layout(){

    if (globalnetwork == null || globalnetwork.links.length == 0){
        alert("There is no network loaded. Please, load a network.");
        return;
    }

    var nodes = globalnetwork.nodes;
    var links = globalnetwork.links;

    // creating elements vector
    var elements = [];


    //{"data": {"id":"Q15349","name":"RPS6KA2","color":"#d1d1e0","localization":"Cytosol"},"group":"nodes"},
    for (var i = 0; i < nodes.length; i++){
        var node = {};

        node.data = {};
        node.data.id = nodes[i].id;
        node.data.name = nodes[i].label;
        node.data.color = "#d1d1e0";
        node.data.localization = nodes[i].attributes["Selected CC"].value.capitalize();

        node.group = "nodes";

        elements.push(node);
    }

    // {"data": {"id":"44","name":"interaction","source":"P00533","target":"P04792","idgroup":"44"},"group":"edges"},
    for (var i = 0; i < links.length; i++){
        var link = {};

        link.data = {};
        link.data.id = links[i].id;
        link.data.name = "interaction";
        link.data.source = links[i].source.id;
        link.data.target = links[i].target.id;
        link.data.idgroup = links[i].id;

        link.group = "edges";

        elements.push(link);
    }

    //console.log(elements);

    // fiding all layers
    var layers = {};
    for (var i = 0; i < nodes.length; i++){
        var cc = nodes[i].attributes["Selected CC"].value;
        if (cc !== null && cc !== ""){
            layers[cc.capitalize()]=0;
        }
    }

    var partial_layers = Object.keys(layers);
    var partial_layers_aux = {};

    for(var i=0; i< partial_layers.length;i++){partial_layers_aux[ partial_layers[i]]="";}


    var final_layers = [];
    if("Extracellular" in partial_layers_aux){
        final_layers.push("Extracellular");
    }
    if("Plasma Membrane" in partial_layers_aux){
        final_layers.push("Plasma Membrane");
    }
    if("Cytoplasm" in partial_layers_aux){
        final_layers.push("Cytoplasm");
    }
    if("Cytosol" in partial_layers_aux){
        final_layers.push("Cytosol");
    }
    if("Nucleus" in partial_layers_aux){
        final_layers.push("Nucleus");
    }
    for(var i = 0; i < partial_layers.length; i++){
        if (final_layers.indexOf(partial_layers[i])==-1){
            final_layers.push(partial_layers[i]);
        }
    }
    //console.log(final_layers);



    options.layers = final_layers;

    //console.log(options.layers);
    // creating layout
    $('#cy2').empty();
    $('#cy2').cytoscape({
                layout: options,
                showOverlay: false,
                zoom: 1,
                style: cerebral_style,
                elements: elements,
                ready: function() {
                    cy = this;
                    cerebral_ready(cy);
                    //console.log("ready");
                }
    });

    d3.select('#cerebralDiv').style("visibility","visible");
}

function changeLabels(){
    //console.log("changeLabels");

    var attribute = d3.select("#change_label").node().value;
    changeNodesLabels(attribute);
}

function changeNodesLabels(attribute){
    //console.log("changeNodesLabels");

    d3.selectAll('.node')//.selectAll(".node_label")
        .each(function (d) {
                d.label = d.attributes[attribute].value;
                d3.select("#label"+d.id).text(d.label);
    });
}

