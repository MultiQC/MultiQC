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
    repeated_apps = []

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
                    repeated_apps.append(a["Program_details"]["program"])

            app_dict[app_name] = a

    except:
        # Used to parse older json files. Will likely be removed in future.
        app_dict = json.loads(f, object_pairs_hook=resolve)
        log.warning("Sample " + name + " uses old json format. Please update to a newer version of HTStream.")
        raise

    return app_dict, repeated_apps


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
