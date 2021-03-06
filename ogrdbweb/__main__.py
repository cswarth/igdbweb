#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Convert IGDB json files into static website.

Usage:
    python -m ogrdbweb

'''
from __future__ import print_function

import os
import sys
import argparse
from flask_frozen import Freezer

from ogrdbweb import app


# pass commandline arguments to flask app
# http://flask.pocoo.org/snippets/133/
#
# Start this app with,
# 	$ python -m ogrdbweb 
   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ogrdbbweb', description=__doc__)

    parser.add_argument('-H', '--host', default="0.0.0.0",
                        help='Hostname of the Flask app [%(default)s]')
    parser.add_argument('-P', '--port', default=5000,
                        help='Port for the Flask app  [%(default)d]')
    parser.add_argument('-d', '--debug', default=False, action="store_true",
                        help='turn on debugging output')

    parser.add_argument('-f', '--freeze', default=False, action="store_true",
                        help="""freeze pages""")
    parser.add_argument('-t', '--template', default='template.jinja',
                        help="""Jinja2 Tempate file[default: %(default)s]""")
    parser.add_argument('-c', '--content', default="content",
                        help="""Directory where content can be found: [default "%(default)s"]""")
    parser.add_argument('-n', '--dryrun', action='store_true',
                        help="""do everythign short of writing to the filesystem.""")
    parser.add_argument('-o', '--output', default='./build',
                        help='directory where static web files will be placed.  Equivalent to `FREEZER_DESTINATION` Frozen-Flask configuration parameter. [default %(default)s]')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output')

    global options
    options = parser.parse_args()
    app.config['OPTIONS'] = options
    app.config['FREEZER_DESTINATION'] = options.output

    if options.freeze:
        freezer = Freezer(app)
        print(freezer.root)
        freezer.freeze()
    else:
        app.run(
            debug=options.debug,
            host=options.host,
            port=int(options.port)
        )
