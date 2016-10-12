#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Writing of D3-format tree to a stream.
"""

import re
import warnings
import json
from dendropy.utility import error
from dendropy.utility import textprocessing
from dendropy.dataio import jsonwriter

##############################################################################
## D3Writer

class D3Writer(jsonwriter.JsonWriter):
    """
    Formatter for D3 data -- a special case of the JsonWriter
    """

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------
        X suppress_leaf_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels will not be rendered for leaves.
            Default is |False|: render leaf taxon labels. See notes below for
            details.
        X suppress_leaf_node_labels : boolean, default: |True|
            If |False|, then node labels (if available) will be printed for
            leaves. Defaults to |True|: do not render leaf node labels. See
            notes below for details.
        suppress_internal_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels will not be printed for internal
            nodes. Default is |False|: print taxon labels for internal nodes.
            See notes below for details.
        suppress_internal_node_labels : boolean, default: |False|
            If |True|, then node labels will not be printed for internal nodes.
            Default is |False|: print node labels for internal nodes. See notes
            below for details.
        X suppress_edge_lengths : boolean, default: |False|
            If |True|, will not write edge lengths. Default is |False|: edge
            lengths will be written.
        store_tree_weights : boolean, default: |False|
            If |True|, tree weights are written. Default is |False|: tree
            weights will not be written.
        X suppress_annotations : boolean, default: |True|
            If |False|, metadata annotations will be written out as special
            comments. Defaults to |True|: metadata annotations will be ignored.
        suppress_item_comments : boolean, default: |True|
            If |False|, any additional comments associated with trees, nodes,
            edges, etc. will be written. Default is |True|: comments will be
            ignored.
        node_label_element_separator : string, default: ' '
            If both ``suppress_leaf_taxon_labels`` and
            ``suppress_leaf_node_labels`` are |False|, then this will be the
            string used to join them. Defaults to ' ' (space).
        node_label_compose_fn : function object or |None|, default: |None|
            If not |None|, should be a function that takes a |Node|
            object as an argument and returns the string to be used to
            represent the node in the tree statement. The return value from
            this function is used unconditionally to print a node
            representation in a tree statement, by-passing the default
            labelling function, ignoring ``suppress_leaf_taxon_labels``,
            ``suppress_leaf_node_labels=True``, ``suppress_internal_taxon_labels``,
            ``suppress_internal_node_labels``, etc. Defaults to |None|.
        edge_label_compose_fn : function object or |None|, default: |None|
            If not |None|, should be a function that takes an Edge object as
            an argument, and returns the string to be used to represent the
            edge length in the tree statement.
        real_value_format_specifier : string, default: ''
            Format specification for real/float values. Will be applied to edge
            lengths (if ``edge_label_compose_fn`` is not given) as well as
            annotations. The format specifier should be given in Python's
            string format specification mini-language. E.g. ".8f", ".4E",
            "8.4f".
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """
        ioservice.DataWriter.__init__(self)
        ioservice.DataWriter.__init__(self, **kwargs)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        pass

    def _write_tree_list(self, stream, tree_list):
        """
        Writes a |TreeList| in D3 schema to ``stream``.
        """
        pass

    def _write_tree(self, stream, tree):
        # Initialize a tree's metadata
        pass
