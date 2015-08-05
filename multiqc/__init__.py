#!/usr/bin/env python

""" MultiQC modules """

class BaseMultiqcModule(object):

    def __init__(self):
        pass

    def dict_to_table(self, d, colheaders=None, tclasses=['table'], sort_rows=False):
        """ Takes a 2D dictionary and creates a HTML table.
        1st set of keys = row headers
        2nd set of keys = column headers
        """
        t = '<table class="{}">'.format(' '.join(tclasses))
        thead = False
        for rh, cols in d.iteritems():
            if not thead:
                t += '<tr><td></td>'
                for c in cols.keys():
                    if colheaders is not None and c in colheaders:
                        t += '<th>{}</th>'.format(colheaders[c])
                    else:
                        t += '<th>{}</th>'.format(c)
                t += '</tr>'
                thead = True

            t += '<tr><th>{}</th><td>{}</td></tr>'.format(rh, '</td><td>'.join(cols.values()))
        t += '</table>'
        return t
