#!/usr/bin/env python

""" MultiQC code to export data to MegaQC / flat JSON files """


import gzip
import io
import json
import os

import requests

from . import config

log = config.logger


# Custom encoder to handle lambda functions
class MQCJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if callable(obj):
            try:
                return obj(1)
            except Exception:
                return None
        return json.JSONEncoder.default(self, obj)


def multiqc_dump_json(report):
    exported_data = dict()
    export_vars = {
        "report": [
            "data_sources",
            "general_stats_data",
            "general_stats_headers",
            "multiqc_command",
            "plot_data",
            "saved_raw_data",
        ],
        "config": [
            "analysis_dir",
            "creation_date",
            "git_hash",
            "intro_text",
            "report_comment",
            "report_header_info",
            "script_path",
            "short_version",
            "subtitle",
            "title",
            "version",
        ],
    }
    for s in export_vars:
        for k in export_vars[s]:
            try:
                d = None
                if s == "config":
                    d = {f"{s}_{k}": getattr(config, k)}
                elif s == "report":
                    d = {f"{s}_{k}": getattr(report, k)}
                if d:
                    json.dumps(d, cls=MQCJSONEncoder, ensure_ascii=False)  # Test that exporting to JSON works
                    exported_data.update(d)
            except (TypeError, KeyError, AttributeError):
                log.warning(f"Couldn't export data key '{s}.{k}'")
        # Get the absolute paths of analysis directories
        exported_data["config_analysis_dir_abs"] = list()
        for d in exported_data.get("config_analysis_dir", []):
            try:
                exported_data["config_analysis_dir_abs"].append(os.path.abspath(d))
            except Exception:
                pass
    return exported_data


def multiqc_api_post(exported_data):
    headers = {"Content-Type": "application/json", "content-encoding": "gzip"}
    if config.megaqc_access_token is not None:
        headers["access_token"] = config.megaqc_access_token
    post_data = json.dumps({"data": exported_data}, cls=MQCJSONEncoder, ensure_ascii=False, indent=2)
    post_data = post_data.encode("utf-8", "ignore")

    # Gzip the JSON for massively decreased filesize
    sio_obj = io.BytesIO()
    gzfh = gzip.GzipFile(fileobj=sio_obj, mode="w")
    gzfh.write(post_data)
    gzfh.close()
    request_body = sio_obj.getvalue()

    log.debug("Sending data to MegaQC")
    log.debug(f"MegaQC URL: {config.megaqc_url}")
    try:
        r = requests.post(config.megaqc_url, headers=headers, data=request_body, timeout=config.megaqc_timeout)
    except (requests.exceptions.ConnectTimeout, requests.exceptions.ReadTimeout) as e:
        log.error(f"Timed out when sending data: {e}")
    except requests.exceptions.ConnectionError:
        log.error(f"Couldn't connect to MegaQC URL {config.megaqc_url}")
    except Exception as e:
        log.error(f"Error sending data: {e}")
    else:
        try:
            api_r = json.loads(r.text)
        except Exception:
            log.error(f"Error: JSON response could not be parsed (status code: {r.status_code})")
            return None
        if r.status_code == 200:
            if api_r["success"]:
                log.info(f"{api_r['message']}")
            else:
                log.error(f"Error - {api_r['message']}")
        else:
            if r.status_code == 403:
                if config.megaqc_access_token is not None:
                    log.error("Error 403: Authentication error, megaqc_access_token not recognised")
                else:
                    log.error("Error 403: Authentication error, megaqc_access_token is required")
            else:
                log.debug(f"MegaQC API status code was {r.status_code}")
                log.error(f"Error - {api_r.get('message', 'Unknown problem')}")
