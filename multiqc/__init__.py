#!/usr/bin/env python

""" MultiQC modules """

class BaseMultiqcModule(object):

    def __init__(self):
        pass

    def dict_to_table(self, d, colheaders=None, table_attrs='class="table"', header_attrs={}, cell_attrs={}):
        """ Takes a 2D dictionary and creates a HTML table.
        1st set of keys = row headers
        2nd set of keys = column headers
        """
        t = '<table {}>'.format(table_attrs)
        thead = False
        for rh, cols in d.iteritems():
            if not thead:
                t += '<thead><tr><th></th>'
                for k in cols.keys():
                    try:
                        a = header_attrs[k]
                    except KeyError:
                        a = ''
                    if colheaders is not None and k in colheaders:
                        t += '<th {}>{}</th>'.format(a, colheaders[k])
                    else:
                        t += '<th {}>{}</th>'.format(a, k)
                t += '</tr></thead><tbody>'
                thead = True

            t += '<tr><th>{}</th>'.format(rh)
            for k, v in cols.iteritems():
                try:
                    a = cell_attrs[k]
                except KeyError:
                    a = ''
                t += '<td {}>{}</td>'.format(a, v)
            t += '</tr>'
        t += '</tbody></table>'
        return t
