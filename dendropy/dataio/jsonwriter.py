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
from dendropy.dataio import ioservice

##############################################################################
## JsonWriter

class JsonWriter(ioservice.DataWriter):
    """
    Formatter for JSON data.
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
        X store_tree_weights : boolean, default: |False|
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
        self.suppress_leaf_taxon_labels = kwargs.pop("suppress_leaf_taxon_labels", False)
        self.suppress_leaf_node_labels = kwargs.pop("suppress_leaf_node_labels", True)
        self.suppress_internal_taxon_labels = kwargs.pop("suppress_internal_taxon_labels", False)
        self.suppress_internal_node_labels = kwargs.pop("suppress_internal_node_labels", False)
        self.suppress_edge_lengths = kwargs.pop("suppress_edge_lengths", False)
        self.preserve_spaces = kwargs.pop("preserve_spaces", False)
        self.store_tree_weights = kwargs.pop("store_tree_weights", False)
        self.suppress_annotations = kwargs.pop("suppress_annotations", True)
        self.suppress_item_comments = kwargs.pop("suppress_item_comments", True)
        self.node_label_element_separator = kwargs.pop("node_label_element_separator", ' ')
        self.node_label_compose_fn = kwargs.pop("node_label_compose_fn", None)
        self.edge_label_compose_fn = kwargs.pop("edge_label_compose_fn", None)
        self._real_value_format_specifier = ""
        self._real_value_formatter = None
        self.real_value_format_specifier = kwargs.pop("real_value_format_specifier", self._real_value_format_specifier)
        if self.edge_label_compose_fn is None:
            self.edge_label_compose_fn = self._format_edge_length
        self.check_for_unused_keyword_arguments(kwargs)

    def _get_real_value_format_specifier(self):
        return self._real_value_format_specifier

    def _set_real_value_format_specifier(self, f):
        if f is None:
            f = ""
        self._real_value_format_specifier = f
        s = "{:" + self._real_value_format_specifier + "}"
        self._real_value_formatter = s.format
    real_value_format_specifier = property(_get_real_value_format_specifier, _set_real_value_format_specifier)

    def _format_edge_length(self, edge):
        """
        Note: instance method to allow overriding.
        """
        return self._real_value_formatter(edge.length)

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):
        for tree_list in tree_lists:
            if (self.attached_taxon_namespace is not None
                    and tree_list.taxon_namespace is not self.attached_taxon_namespace):
                continue
            self._write_tree_list(stream, tree_list)

    def _write_tree_list(self, stream, tree_list):
        """
        Writes a |TreeList| in D3 schema to ``stream``.
        """
        # Write out a tree list. If only one tree is in the list, write out tree only.
        if len(tree_list) > 1:
            metadata = dict(name=tree_list.label)
            # Add treelist annotations
            if self.suppress_annotations is False and tree_list.has_annotations:
                metadata.update(annotations=tree_list.annotations.values_as_dict())
            # Start adding tree.
            trees = []
            for tree in tree_list:
                tree_metadata = {}
                self._get_tree_metadata(tree, tree_metadata)
                trees.append(tree_metadata)
            stream.write(json.dumps(metadata))
            stream.write("\n")
        else:
            self._write_tree(stream, tree_list[0])

    def _write_tree(self, stream, tree):
        # Initialize a tree's metadata
        metadata = {}
        self._get_tree_metadata(tree, metadata)
        tree_string = json.dumps(metadata)
        stream.write(tree_string)

    def _get_tree_metadata(self, tree, out):
        """Add tree metadata to out.
        """
        out.update(name=tree.label)
        # Add annotations to tree.
        if self.suppress_annotations is False and tree.has_annotations:
            out.update(annotations=tree.annotations.values_as_dict())
        # Add tree weight.
        if self.store_tree_weights is True and tree.weight is not None:
            out.update(weight=tree.weight)
        # Add node metadata to tree:
        tree_metadata = {}
        for node in tree.nodes():
            # Iterate over deepest nodes.
            if node.parent_node is None:
                self._get_node_metadata(node, tree_metadata)
        out.update(data=tree_metadata)

    def _get_node_metadata(self, node, out):
        """Add node metadata to out.
        """
        # Add annotations
        if self.suppress_annotations is False and node.has_annotations:
            out.update(annotations=node.annotations.values_as_dict())
        # Add edge lengths to children nodes
        if self.suppress_edge_lengths is False:
            if node.edge_length is not None:
                out.update(length=node.edge_length)
        # Handle nodes and leafs differently
        if node.is_leaf():
            # Add labels
            if self.suppress_leaf_node_labels is False:
                out.update(name=node.label)
            if self.suppress_leaf_taxon_labels is False and node.taxon is not None:
                self._get_taxon_metadata(node.taxon, out)
        else:
            # Add labels to node
            if self.suppress_internal_node_labels is False:
                out.update(name=node.label)
            # Add taxon info to node
            if self.suppress_internal_taxon_labels is False and node.taxon is not None:
                self._get_taxon_metadata(node.taxon, out)
            # Get children data
            self._get_children_metadata(node, out)

    def _get_taxon_metadata(self, taxon, out):
        """Adds new attribute to node out. Key=`taxon`, value=taxon metadata
        """
        metadata = dict(name=taxon.label)
        # Add annotations
        if self.suppress_annotations is False and taxon.has_annotations:
            metadata.update(annotations=taxon.annotations.values_as_dict())
        out.update(taxon=metadata)

    def _get_children_metadata(self, parent, out):
        """Construct a `children` list. Get childrend metadata. Add to out.
        """
        children = []
        nodes = parent.child_nodes()
        for child in nodes:
            metadata = {}
            # Get node data for child.
            self._get_node_metadata(child, metadata)
            # Add parent label to each child
            if self.suppress_internal_node_labels is False:
                metadata.update(parent=parent.label)
            children.append(metadata)
        out.update(children=children)
