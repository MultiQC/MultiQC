""" MultiQC code to export data to MegaQC / flat JSON files """

import gzip
import io
import json

import requests

from . import config
from .util_functions import dump_json

log = config.logger


def multiqc_api_post(exported_data):
    headers = {"Content-Type": "application/json", "content-encoding": "gzip"}
    if config.megaqc_access_token is not None:
        headers["access_token"] = config.megaqc_access_token
    post_data = dump_json({"data": exported_data}, ensure_ascii=False, indent=2)
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
