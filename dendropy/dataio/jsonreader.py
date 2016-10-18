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
Parsing of JSON-format tree from a stream.
"""

import re
import warnings
import json
from dendropy.utility import error
from dendropy.utility import deprecate
from dendropy.utility.textprocessing import StringIO
from dendropy.dataio import ioservice
from dendropy.datamodel.treemodel import Node
from dendropy.datamodel.taxonmodel import Taxon
##############################################################################
## JsonReader

class JsonReader(ioservice.DataReader):
    """
    Parser for JSON-formatted data.
    """


    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------
        rooting : string, {['default-unrooted'], 'default-rooted', 'force-unrooted', 'force-rooted'}
            Specifies how trees in the data source should be intepreted with
            respect to their rooting:

                'default-unrooted' [default]:
                    All trees are interpreted as unrooted unless a '[&R]'
                    comment token explicitly specifies them as rooted.
                'default-rooted'
                    All trees are interpreted as rooted unless a '[&U]'
                    comment token explicitly specifies them as unrooted.
                'force-unrooted'
                    All trees are unconditionally interpreted as unrooted.
                'force-rooted'
                    All trees are unconditionally interpreted as rooted.

        edge_length_type : type, default: ``float``
            Specifies the type of the edge lengths (``int`` or ``float``). Tokens
            interpreted as branch lengths will be cast to this type.
            Defaults to ``float``.
        suppress_edge_lengths : boolean, default: |False|
            If |True|, edge length values will not be processed. If |False|,
            edge length values will be processed.
        extract_comment_metadata : boolean, default: |True|
            If |True| (default), any comments that begin with '&' or '&&' will
            be parsed and stored as part of the annotation set of the
            corresponding object (accessible through the ``annotations``
            attribute of the object). This requires that the comment
            contents conform to a particular format (NHX or BEAST: 'field =
            value'). If |False|, then the comments will not be parsed,
            but will be instead stored directly as elements of the ``comments``
            list attribute of the associated object.
        store_tree_weights : boolean, default: |False|
            If |True|, process the tree weight (e.g. "[&W 1/2]") comment
            associated with each tree, if any. Defaults to |False|.
        encode_splits : boolean, default: |False|
            If |True|, split hash bitmasks will be calculated and attached to
            the edges.
        finish_node_fn : function object, default: |None|
            If specified, this function will be applied to each node after
            it has been constructed.
        case_sensitive_taxon_labels : boolean, default: |False|
            If |True|, then taxon labels are case sensitive (e.g., "P.regius"
            and "P.REGIUS" wil be treated as different operation taxonomic
            unit concepts). Otherwise, taxon label intepretation will be made
            without regard for case.
        preserve_underscores : boolean, default: |False|
            If |True|, unquoted underscores in labels will *not* converted to
            spaces. Defaults to |False|: all underscores not protected by
            quotes will be converted to spaces.
        suppress_internal_node_taxa : boolean, default: |True|
            If |False|, internal node labels will be instantantiated into
            |Taxon| objects. If |True|, internal node labels
            will *not* be instantantiated as strings.
        suppress_leaf_node_taxa : boolean, default: |False|
            If |False|, leaf (external) node labels will be instantantiated
            into |Taxon| objects. If |True|, leaff (external) node
            labels will *not* be instantantiated as strings.
        terminating_semicolon_required : boolean, default: |True|
            If |True| [default], then a tree statement that does not end in a
            semi-colon is an error. If |False|, then no error will be raised.
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.
        """
        # base
        ioservice.DataReader.__init__(self)
        self.edge_length_type = kwargs.pop("edge_length_type", float)
        self.suppress_edge_lengths = kwargs.pop("suppress_edge_lengths", False)
        self.extract_comment_metadata = kwargs.pop('extract_comment_metadata', True)
        self.store_tree_weights = kwargs.pop("store_tree_weights", False)
        #self.default_tree_weight = kwargs.pop("default_tree_weight", self.__class__._default_tree_weight)
        self.finish_node_fn = kwargs.pop("finish_node_fn", None)
        self.case_sensitive_taxon_labels = kwargs.pop('case_sensitive_taxon_labels', False)
        self.preserve_unquoted_underscores = kwargs.pop('preserve_underscores', False)
        self.suppress_internal_node_taxa = kwargs.pop("suppress_internal_node_taxa", True)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_external_node_taxa", False) # legacy (will be deprecated)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_leaf_node_taxa", self.suppress_leaf_node_taxa)
        self.terminating_semicolon_required = kwargs.pop("terminating_semicolon_required", True)
        self.check_for_unused_keyword_arguments(kwargs)

        # per-tree book-keeping
        self._tree_statement_complete = None
        self._parenthesis_nesting_level = None
        self._seen_taxa = None

    def tree_iter(self,
            stream,
            tree_factory):
        """
        Iterator that yields trees in NEWICK-formatted source.

        Parameters
        ----------
        stream : file or file-like object
            A file or file-like object opened for reading.
        taxon_namespace : |TaxonNamespace|
            Operational taxonomic unit namespace to use for taxon management.
        tree_factory : function object
            A function that returns a new |Tree| object when called
            without arguments.

        Returns
        -------
        iter : :py`collections.Iterator` [|Tree|]
            An iterator yielding |Tree| objects constructed based on
            data in ``stream``.
        """
        treelist_metadata = json.loads(stream.getvalue())
        if type(treelist_metadata) is not list:
            treelist_metadata = [treelist_metadata]
        #while True:
        for tree_metadata in treelist_metadata:
            tree = self._parse_tree_metadata(
                    tree_metadata=tree_metadata,
                    tree_factory=tree_factory)
            yield tree
            if tree is None:
                raise StopIteration

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        taxon_namespace = taxon_namespace_factory(label=None)
        tree_list = tree_list_factory(label=None, taxon_namespace=taxon_namespace)
        tree_factory = tree_list.new_tree
        for tree in self.tree_iter(stream=stream,
                tree_factory=tree_factory):
            pass
        product = self.Product(
                taxon_namespaces=None,
                tree_lists=[tree_list],
                char_matrices=None)
        return product


    def _parse_tree_metadata(self,
            tree_metadata,
            tree_factory):
        """
        Parses a single tree statement from a token stream and constructs a
        corresponding Tree object.
        """
        # Get the data of the seed
        seed_node = tree_metadata.get("data")
        # Parse node
        node = self._parse_node_metadata(seed_node)
        # Add node to tree
        tree = tree_factory(seed_node=node, label=tree_metadata.get("name"))
        # Add annotations.
        if tree_metadata.get("annotations") is not None:
            self._parse_annotations(tree, tree_metadata.get("annotations"))
        return tree

    def _parse_node_metadata(self, node_metadata):
        """Parse node json and construct a node object"""
        # Add taxon.
        kwargs = dict(label=node_metadata.get("name"),
            edge_length=node_metadata.get("length"))
        # Parse taxon
        if node_metadata.get("taxon") is not None:
            taxon_metadata = node_metadata.get("taxon")
            kwargs.update(taxon=Taxon(label=taxon_metadata.get("name")))
        # Add node.
        node = Node(**kwargs)
        if node_metadata.get("children") is not None:
            self._parse_children_metadata(node, node_metadata.get("children"))
        # Add annotations.
        if node_metadata.get("annotations") is not None:
            self._parse_annotations(node, node_metadata.get("annotations"))
        return node

    def _parse_children_metadata(self, parent_node, children_metadata):
        """Parse multiple children json and construct children nodes underneath
        parent node.
        """
        for child in children_metadata:
            node = self._parse_node_metadata(child)
            parent_node.add_child(node)

    def _parse_annotations(self, obj, annotations):
        """Parse annotations json and construct annotations object."""
        for key, value in annotations.items():
            obj.annotations.add_new(key, value)
