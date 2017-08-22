#!/usr/bin/env python

""" MultiQC code to export data to MegaQC / flat JSON files """

from __future__ import print_function
import json
import requests

from multiqc import config
log = config.logger

# Custom encoder to handle lambda functions
class MQCJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if callable(obj):
            try:
                return obj(1)
            except:
                return None
        return json.JSONEncoder.default(self, obj)

def multiqc_dump_json(report):
    exported_data = dict()
    export_vars = {
        'report': [
            'data_sources',
            'general_stats_data',
            'general_stats_headers',
            'multiqc_command',
            'plot_data',
            'saved_raw_data',
        ],
        'config': [
            'analysis_dir',
            'creation_date',
            'git_hash',
            'intro_text',
            'report_comment',
            'report_header_info',
            'script_path',
            'short_version',
            'subtitle',
            'title',
            'version',
        ]
    }
    for s in export_vars:
        for k in export_vars[s]:
            try:
                if s == 'config':
                    d = {'{}_{}'.format(s, k): getattr(config, k)}
                elif s == 'report':
                    d = {'{}_{}'.format(s, k): getattr(report, k)}
                json.dumps(d, cls=MQCJSONEncoder, ensure_ascii=False) # Test that exporting to JSON works
                exported_data.update(d)
            except (TypeError, KeyError, AttributeError):
                log.warn("Couldn't export data key '{}.{}'".format(s, k))
    return exported_data


def multiqc_api_post(exported_data):
    headers = { 'Content-Type': 'application/json' }
    if config.megaqc_access_token is not None:
        headers['access_token'] = config.megaqc_access_token
    post_data = json.dumps({'data': exported_data}, cls=MQCJSONEncoder, ensure_ascii=False)
    log.info("Sending data to MegaQC")
    try:
        r = requests.post(config.megaqc_url, headers=headers, data=post_data, timeout=config.megaqc_timeout)
    except (requests.exceptions.ConnectTimeout, requests.exceptions.ReadTimeout) as e:
        log.warn("Timed out when sending data: {}".format(e))
    except requests.exceptions.ConnectionError:
        log.warn("Couldn't connect to MegaQC URL {}".format(config.megaqc_url))
    except Exception as e:
        log.warn("Error sending data: {}".format(e))
    else:
        try:
            api_r = json.loads(r.text)
        except KeyError:
            api_r = {'message': 'JSON response could not be parsed: {}'.format(r.text)}
        if r.status_code == 200:
            log.info('{}'.format(api_r['message']))
        else:
            if r.status_code == 403:
                if config.megaqc_access_token is not None:
                    log.warn('Error 403: Authentication error, megaqc_access_token not recognised')
                else:
                    log.warn('Error 403: Authentication error, megaqc_access_token is required')
            else:
                log.warn('Error {}: {}'.format(r.status_code, api_r['message']))


