import json
import logging

#################################################

""" Utilities for HTStream submodule """

#################################################

# Logger Initialization
log = logging.getLogger(__name__)


###################################
# convert json
def resolve(pairs):
    resolved_dict = {}
    index_dict = {}

    # iterates through json key value pairs, resolves key conflits
    for k, v in pairs:
        if k in index_dict.keys() and "hts_" in k:
            resolved_dict[k + "_" + str(index_dict[k])] = v
            index_dict[k] += 1

        elif "hts_" in k:
            resolved_dict[k + "_1"] = v
            index_dict[k] = 2

        else:
            resolved_dict[k] = v

    return resolved_dict


###################################
# Json and stats parsing functions
def parse_json(name, f):
    app_dict = {}
    apps = json.loads(f)

    # Will fail if old format is usef
    try:
        # Allows for multiple instances of app, just adds number suffix
        for a in apps:
            i = 1
            app_name = a["Program_details"]["program"] + "_" + str(i)

            if app_name in app_dict.keys():
                while app_name in app_dict.keys():
                    i += 1
                    app_name = a["Program_details"]["program"] + "_" + str(i)

            app_dict[app_name] = a

    except:
        # Used to parse older json files. Will likely be removed in future.
        app_dict = json.loads(f, object_pairs_hook=resolve)
        log.warning("Sample " + name + " uses old json format. Please update to a newer version of HTStream.")
        raise

    return app_dict


###################################
# prints keys in a pretty way
def key_print(dictionary):
    string = ""

    for key in dictionary.keys():
        string += key + ", "

    string = string[:-2] + "."

    return string


###################################
# Checks if read lengths are uniform
def uniform(json, read):
    midpoint = 0

    # Check if read lengths are uniform across all samples
    for key in json.keys():
        temp = json[key][read][0]["shape"][-1] * 2

        if midpoint == 0:
            midpoint = temp

        elif midpoint == temp:
            midpoint = midpoint

        else:
            midpoint = -1
            break

    return midpoint


###################################
# Multiplot html formatter
def multi_plot_html(header, samples, btn_1, btn_2, id_1, id_2, graph_1, graph_2, exempt=True):
    # section header
    wrapper_html = header

    # Buttons
    wrapper_html += '<div class="btn-group hc_switch_group {}">\n'.format("htstream_exempt")
    wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this, \'{t}\')" id="{i}_btn">{b}</button>\n'.format(
        i=id_1, t=id_2, b=btn_1
    )
    wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this, \'{t}\')" id="{i}_btn">{b}</button>\n'.format(
        i=id_2, t=id_1, b=btn_2
    )
    wrapper_html += "</div>\n"

    # this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
    # line graph div
    wrapper_html += '<div id="{b}" class="htstream_fadein">'.format(b=id_1)
    wrapper_html += graph_1 + "</div>"

    # get unique id and create unique id for dropdown
    id_suffix = id_1.split("_")[-1]
    id_3 = "htstream_stats_dropdown_" + id_suffix

    # split up graph html
    graph_2 = graph_2.split("\n")
    tmp = graph_2[0].split(">")

    # # add dropdown html
    tmp[0] += """><div class="btn-group">
                <button type="button" class="btn btn-default dropdown-toggle" id="{i}" data-toggle="dropdown">{s} <span class="caret"></span></button>
                    <ul class="dropdown-menu scrollable-menu" role="menu">""".format(i=id_3, s=samples[0])

    # populate dropdown buttons for each sample
    for s in samples:
        tmp[0] += """<li style="border-bottom: 1px solid #bdbcbc; margin-botton:8px;">
                        <a><button class="hts-btn" id="{i}" onclick="hts_btn_click(this)">{s}</button></a>
                     </li>""".format(i=s + "_" + id_suffix, s=s)

    # close all elements and rejoin
    tmp[0] += """</ul>
            </div"""

    # hide buttons
    tmp[1] += ' style="display:none;"'

    # rejoin all html
    graph_2[0] = ">".join(tmp)
    graph_2 = "\n".join(graph_2)

    # this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
    # line graph div
    wrapper_html += '<div id="{b}" class="htstream_fadein" style="display:none;"><br>'.format(b=id_2)
    wrapper_html += graph_2 + "</div>"

    return wrapper_html
