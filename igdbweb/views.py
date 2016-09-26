#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import os
import sys

import gene
import datetime
import getpass
import hashlib
from cStringIO import StringIO
from jinja2 import Environment, FileSystemLoader
from flask import render_template, abort, Response, redirect, url_for, request, g, jsonify
from flask_breadcrumbs import register_breadcrumb

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from igdbweb import app

@app.route('/index')
@app.route("/")
@register_breadcrumb(app, '.', 'Index')
def index():
    genes = app.config['GENES']
    print("generating index page for {}".format(genes))
    renderdict = {
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'IgDB',
        'genes': genes
        }

    return render_template('index.html', **renderdict)

def view_gene_dlc():
    id = request.view_args['id']
    genes = app.config['GENES']
    gene = genes[id]
    return [{'text': id, 'url': url_for('gene_page', id=id)}]

@app.route("/gene/<id>.html")
@register_breadcrumb(app, '.gene', '<ignored>',
                     dynamic_list_constructor=view_gene_dlc)
def gene_page(id=None):
    print("generating gene page for {}".format(id))
    genes = app.config['GENES']
    gene = genes[id]
    
    renderdict = {
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'IgDB',
        }
    
    renderdict.update(gene.__dict__)
    return render_template('gene.html', **renderdict)
            

@app.route("/fasta/<id>.fa")
def gene_fasta(id=None):
    print("generating fasta page for {}".format(id))
    genes = app.config['GENES']
    gene = genes[id]

    def to_fasta(seq):
        fp = StringIO()
        record = SeqRecord(Seq(seq, IUPAC.extended_dna),
                                   id=id, description="")
        SeqIO.write(record, fp, 'fasta')
        return fp.getvalue()

        
    fasta = to_fasta(gene.sequence)
    return Response(fasta, mimetype="application/octet-stream")



@app.route("/phylip/<id>.phy")
def gene_phylip(id=None):
    print("generating phylip page for {}".format(id))
    genes = app.config['GENES']
    gene = genes[id]

    def to_phylip(seq):
        fp = StringIO()
        seq = seq.replace('.', '-')
        record = SeqRecord(Seq(seq, IUPAC.extended_dna),
                                   id=id)
        SeqIO.write(record, fp, 'phylip-relaxed')
        return fp.getvalue()

    phylip = to_phylip(gene.sequence)

    return Response(phylip, mimetype="application/octet-stream")


            

