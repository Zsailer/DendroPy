#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests basic Tree structure and iteration.
"""

import unittest
import dendropy

###############################################################################
## Test Tree
# Following tree:
#
#                  a
#                 / \
#                /   \
#               /     \
#              /       \
#             /         \
#            /           \
#           /             c
#          b             / \
#         / \           /   \
#        /   e         /     f
#       /   / \       /     / \
#      /   /   \     g     /   h
#     /   /     \   / \   /   / \
#    i    j     k  l  m  n   o   p
#
#  Can be specified as:
#
#      a -> b -> i;
#      b -> e -> j;
#      e -> k;
#      a -> c;
#      c -> g;
#      c -> f;
#      g -> l;
#      g -> m;
#      f -> n;
#      f -> h -> o;
#      h -> p;
class TestTreeStructure(unittest.TestCase):
    dot_str = "a -> b -> i; b -> e -> j; e -> k; a -> c; c -> g; c -> f; g -> l; g -> m; f -> n; f -> h -> o; h -> p;"
    newick_unweighted_edges_str = "((i, (j, k)e)b, ((l, m)g, (n, (o, p)h)f)c)a;"
    newick_weighted_edges_str = "((i:1, (j:2, k:3)e:4)b:5, ((l:6, m:7)g:8, (n:9, (o:10, p:11)h:12)f:13)c:14)a:15;"
    preorder_sequence = ["a", "b", "i", "e", "j", "k", "c", "g", "l", "m", "f", "n", "h", "o", "p"]
    postorder_sequence = ["i", "j", "k", "e", "b", "l", "m", "g", "n", "o", "p", "h", "f", "c", "a"]
    leaf_sequence = ["i", "j", "k", "l", "m", "n", "o", "p"]
    levelorder_sequence = ["a", "b", "c", "i", "e", "g", "f", "j", "k", "l", "m", "n", "h", "o", "p"]
    internal_levelorder_sequence = ["a", "bc", "egf", "h"]
    inorder_sequence = ["i", "b", "j", "e", "k", "a", "l", "g", "m", "c", "n", "f", "o", "h", "p"]
    ageorder_sequence = ["i", "j", "k", "l", "m", "n", "o", "p", "e", "g", "h", "b", "f", "c", "a"]
    node_children = {
            "a" : ["b", "c"],
            "b" : ["i", "e"],
            "c" : ["g", "f"],
            "e" : ["j", "k"],
            "f" : ["n", "h"],
            "g" : ["l", "m"],
            "h" : ["o", "p"],
            "i" : [],
            "j" : [],
            "k" : [],
            "l" : [],
            "m" : [],
            "n" : [],
            "o" : [],
            "p" : [],
            }
    node_siblings = {
            "a": [],
            "b": ["c"],
            "c": [],
            "e": [],
            "f": [],
            "g": ["f"],
            "h": [],
            "i": ["e"],
            "j": ["k"],
            "k": [],
            "l": ["m"],
            "m": [],
            "n": ["h"],
            "o": ["p"],
            "p": [],
            }
    node_edge_lengths = {
            "a": 15.0,
            "b": 33.0,
            "c": 14.0,
            "e": 14.0,
            "f": 13.0,
            "g": 30.0,
            "h": 12.0,
            "i": 17.0,
            "j":  3.0,
            "k":  3.0,
            "l":  6.0,
            "m":  6.0,
            "n": 23.0,
            "o": 11.0,
            "p": 11.0,
            }
    node_ages = {
            "a": 50.0,
            "b": 17.0,
            "c": 36.0,
            "e":  3.0,
            "f": 23.0,
            "g":  6.0,
            "h": 11.0,
            "i":  0.0,
            "j":  0.0,
            "k":  0.0,
            "l":  0.0,
            "m":  0.0,
            "n":  0.0,
            "o":  0.0,
            "p":  0.0,
            }
    node_ancestors = {
            "a": [],
            "b": ["a"],
            "c": ["a"],
            "e": ["b", "a"],
            "f": ["c", "a"],
            "g": ["c", "a"],
            "h": ["f", "c", "a"],
            "i": ["b", "a"],
            "j": ["e", "b", "a"],
            "k": ["e", "b", "a"],
            "l": ["g", "c", "a"],
            "m": ["g", "c", "a"],
            "n": ["f", "c", "a"],
            "o": ["h", "f", "c", "a"],
            "p": ["h", "f", "c", "a"],
            }

    # def get_tree(self):
    #     tree = dendropy.Tree()
    #     def add_child_node(parent, label, edge_length):
    #         nd = tree.node_factory()
    #         nd.label = label
    #         nd.edge.length = edge_length
    #         parent.add_child(nd)
    #         return nd
    #     a = tree.seed_node
    #     a.label = "a"
    #     a.edge.length = 15.0
    #     b = add_child_node(a, label="b", edge_length=xxx["b"])
    #     c = add_child_node(a, label="c", edge_length=xxx["c"])
    #     e = add_child_node(b, label="i", edge_length=xxx["i"])
    #     i = add_child_node(b, label="c", edge_length=xxx["c"])
    #     j = add_child_node(e, label="j", edge_length=xxx["j"])
    #     k = add_child_node(e, label="k", edge_length=xxx["k"])
    #     f = add_child_node(c, label="f", edge_length=xxx["f"])
    #     g = add_child_node(c, label="g", edge_length=xxx["g"])
    #     l = add_child_node(g, label="l", edge_length=xxx["l"])
    #     m = add_child_node(g, label="m", edge_length=xxx["m"])
    #     h = add_child_node(f, label="h", edge_length=xxx["h"])
    #     n = add_child_node(f, label="n", edge_length=xxx["n"])
    #     o = add_child_node(h, label="o", edge_length=xxx["o"])
    #     p = add_child_node(h, label="p", edge_length=xxx["p"])
    #     return tree

    def get_tree(self):
        tree = dendropy.Tree()
        a = tree.seed_node
        a.label = "a"
        a.edge.length = 15.0
        b = a.new_child(label="b", edge_length=self.node_edge_lengths["b"])
        assert b.label == "b"
        assert b.edge.length == self.node_edge_lengths[b.label]
        assert b.parent_node is a
        assert b.edge.tail_node is a
        assert b in a._child_nodes
        c = a.new_child(label="c", edge_length=self.node_edge_lengths["c"])
        assert c.label == "c"
        assert c.edge.length == self.node_edge_lengths[c.label]
        assert c.parent_node is a
        assert c.edge.tail_node is a
        assert c in a._child_nodes
        i = b.new_child(label="i", edge_length=self.node_edge_lengths["i"])
        assert i.label == "i"
        assert i.edge.length == self.node_edge_lengths[i.label]
        assert i.parent_node is b
        assert i.edge.tail_node is b
        assert i in b._child_nodes
        e = b.new_child(label="e", edge_length=self.node_edge_lengths["e"])
        assert e.label == "e"
        assert e.edge.length == self.node_edge_lengths[e.label]
        assert e.parent_node is b
        assert e.edge.tail_node is b
        assert e in b._child_nodes
        j = e.new_child(label="j", edge_length=self.node_edge_lengths["j"])
        assert j.label == "j"
        assert j.edge.length == self.node_edge_lengths[j.label]
        assert j.parent_node is e
        assert j.edge.tail_node is e
        assert j in e._child_nodes
        k = e.new_child(label="k", edge_length=self.node_edge_lengths["k"])
        assert k.label == "k"
        assert k.edge.length == self.node_edge_lengths[k.label]
        assert k.parent_node is e
        assert k.edge.tail_node is e
        assert k in e._child_nodes
        g = c.new_child(label="g", edge_length=self.node_edge_lengths["g"])
        assert g.label == "g"
        assert g.edge.length == self.node_edge_lengths[g.label]
        assert g.parent_node is c
        assert g.edge.tail_node is c
        assert g in c._child_nodes
        f = c.new_child(label="f", edge_length=self.node_edge_lengths["f"])
        assert f.label == "f"
        assert f.edge.length == self.node_edge_lengths[f.label]
        assert f.parent_node is c
        assert f.edge.tail_node is c
        assert f in c._child_nodes
        l = g.new_child(label="l", edge_length=self.node_edge_lengths["l"])
        assert l.label == "l"
        assert l.edge.length == self.node_edge_lengths[l.label]
        assert l.parent_node is g
        assert l.edge.tail_node is g
        assert l in g._child_nodes
        m = g.new_child(label="m", edge_length=self.node_edge_lengths["m"])
        assert m.label == "m"
        assert m.edge.length == self.node_edge_lengths[m.label]
        assert m.parent_node is g
        assert m.edge.tail_node is g
        assert m in g._child_nodes
        n = f.new_child(label="n", edge_length=self.node_edge_lengths["n"])
        assert n.label == "n"
        assert n.edge.length == self.node_edge_lengths[n.label]
        assert n.parent_node is f
        assert n.edge.tail_node is f
        assert n in f._child_nodes
        h = f.new_child(label="h", edge_length=self.node_edge_lengths["h"])
        assert h.label == "h"
        assert h.edge.length == self.node_edge_lengths[h.label]
        assert h.parent_node is f
        assert h.edge.tail_node is f
        assert h in f._child_nodes
        o = h.new_child(label="o", edge_length=self.node_edge_lengths["o"])
        assert o.label == "o"
        assert o.edge.length == self.node_edge_lengths[o.label]
        assert o.parent_node is h
        assert o.edge.tail_node is h
        assert o in h._child_nodes
        p = h.new_child(label="p", edge_length=self.node_edge_lengths["p"])
        assert p.label == "p"
        assert p.edge.length == self.node_edge_lengths[p.label]
        assert p.parent_node is h
        assert p.edge.tail_node is h
        assert p in h._child_nodes
        tree._debug_check_tree()
        leaf_nodes = set([i, j, k, l, m, n, o, p])
        internal_nodes = set([b, c, e, f, g, h])
        all_nodes = leaf_nodes | internal_nodes | set([a])
        return tree, all_nodes, leaf_nodes, internal_nodes

    ###########################################################################
    ## Node and Edge Collection Access

    def test_get_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.nodes()
        self.assertEqual(len(nodes), len(anodes))
        self.assertEqual(set(nodes), anodes)
        obs_labels = [nd.label for nd in nodes]

    def test_get_nodes_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.nodes(filter_fn = lambda x : x.edge.length > 13)
        exp_nodes = set([nd for nd in anodes if nd.edge.length > 13])
        for nd in nodes:
            self.assertTrue(nd.edge.length > 13)
        self.assertEqual(len(nodes), len(exp_nodes))
        self.assertEqual(set(nodes), exp_nodes)

    def test_get_leaf_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.leaf_nodes()
        self.assertEqual(len(nodes), len(lnodes))
        self.assertEqual(set(nodes), lnodes)

    def test_get_internal_nodes_with_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.internal_nodes()
        nlnodes = inodes | set([tree.seed_node])
        self.assertEqual(len(nodes), len(nlnodes))
        self.assertEqual(set(nodes), nlnodes)

    def test_get_internal_nodes_no_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = tree.internal_nodes(True)
        self.assertEqual(len(nodes), len(inodes))
        self.assertEqual(set(nodes), inodes)

    def test_get_edges(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.edges()
        eset = set([nd.edge for nd in anodes])
        self.assertEqual(len(edges), len(eset))
        self.assertEqual(set(edges), eset)

    def test_get_edges_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.edges(filter_fn=lambda x : x.length > 13)
        exp_edges = set([nd.edge for nd in anodes if nd.edge.length > 13])
        for edge in edges:
            self.assertTrue(edge.length > 13)
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_leaf_edges(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.leaf_edges()
        exp_edges = set([nd.edge for nd in lnodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_internal_edges_with_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.internal_edges()
        nlnodes = inodes | set([tree.seed_node])
        exp_edges = set([nd.edge for nd in nlnodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_internal_edges_no_root(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        edges = tree.internal_edges(True)
        exp_edges = set([nd.edge for nd in inodes])
        self.assertEqual(len(edges), len(exp_edges))
        self.assertEqual(set(edges), exp_edges)

    def test_get_child_nodes(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for node in tree:
            child_labels = [ch.label for ch in node.child_nodes()]
            expected_children = self.node_children[node.label]
            self.assertEqual(len(child_labels), len(expected_children))
            self.assertEqual(set(child_labels), set(expected_children))

    ###########################################################################
    ## (Taxon-free) Node Finders

    def test_find_node(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node(lambda x: x.label == "c")
        self.assertEqual(node.label, "c")

    def test_find_node_nonexisting(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node(lambda x: x.label == "zzz")
        self.assertIs(node, None)

    def test_find_node_with_label(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node_with_label("c")
        self.assertEqual(node.label, "c")

    def test_find_node_with_label_nonexisting(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        node = tree.find_node_with_label("zzz")
        self.assertIs(node, None)

    ###########################################################################
    ## Iterators

    ### Default Iterator ###

    def test_default_iteration(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.preorder_sequence)

    ### Preorder Node Iterator ###

    def test_preorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.preorder_sequence)

    def test_preorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Preorder Internal Node Iterator ###

    def test_preorder_internal_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_internal_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x]]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_internal_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.preorder_internal_node_iter(exclude_seed_node=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_node_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.preorder_internal_node_iter(exclude_seed_node=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    ### Postorder Node Iterator ###

    def test_postorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.postorder_sequence)

    def test_postorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Postorder Internal Node Iterator ###

    def test_postorder_internal_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_internal_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x]]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_internal_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.postorder_internal_node_iter(exclude_seed_node=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_node_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.postorder_internal_node_iter(exclude_seed_node=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    ### Level-Order Node Iterator ###

    def test_levelorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.levelorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.levelorder_sequence)

    def test_levelorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.levelorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.levelorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### In-Order Node Iterator ###

    def test_inorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.inorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.inorder_sequence)

    def test_inorder_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.inorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.inorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Leaf Node Iterator ###

    def test_leaf_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.leaf_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.leaf_sequence)

    def test_leaf_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.leaf_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.leaf_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Age-Order Node Iterator ###

    def test_node_ages(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        tree.calc_node_ages()
        nodes = [nd for nd in tree.ageorder_node_iter()]
        for nd in nodes:
            self.assertEqual(nd.age, self.node_ages[nd.label])

    def test_ageorder_node_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.ageorder_sequence)

    def test_ageorder_node_iter_unfiltered_no_leaves(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter(include_leaves=False)]
        visited_labels = [nd.label for nd in nodes]
        expected = [label for label in self.ageorder_sequence if self.node_children[label]]
        self.assertEqual(visited_labels, expected)

    def test_ageorder_node_iter_unfiltered_reversed(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [nd for nd in tree.ageorder_node_iter(descending=True)]
        visited_labels = [nd.label for nd in nodes]
        nda = [ (self.node_ages[x], x) for x in self.preorder_sequence ]
        nda.sort(key=lambda x: x[0], reverse=True)
        exp = [x[1] for x in nda]
        self.assertEqual(visited_labels, exp)

    def test_leaf_node_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.edge.length > 13
        nodes = [nd for nd in tree.ageorder_node_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.ageorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Preorder Edge Iterator ###

    def test_preorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.preorder_sequence)

    def test_preorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Preorder Internal Edge Iterator ###

    def test_preorder_internal_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x]]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(exclude_seed_edge=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    def test_preorder_internal_edge_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.preorder_internal_edge_iter(exclude_seed_edge=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.preorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    ### Postorder Edge Iterator ###

    def test_postorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.postorder_sequence)

    def test_postorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Postorder Internal Edge Iterator ###

    def test_postorder_internal_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x]]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13)]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_without_root_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                self.node_children[x] and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    def test_postorder_internal_edge_iter_without_root_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.postorder_internal_edge_iter(exclude_seed_edge=True, filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.postorder_sequence if
                (self.node_children[x] and self.node_edge_lengths[x] > 13) and x != "a"]
        self.assertEqual(visited_labels, exp_labels)

    ### Level-Order Edge Iterator ###

    def test_levelorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.levelorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.levelorder_sequence)

    def test_levelorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.levelorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.levelorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### In-Order Edge Iterator ###

    def test_inorder_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.inorder_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.inorder_sequence)

    def test_inorder_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.inorder_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.inorder_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ### Leaf Edge Iterator ###

    def test_leaf_edge_iter_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        nodes = [edge.head_node for edge in tree.leaf_edge_iter()]
        visited_labels = [nd.label for nd in nodes]
        self.assertEqual(visited_labels, self.leaf_sequence)

    def test_leaf_edge_iter_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        f = lambda x: x.length > 13
        nodes = [edge.head_node for edge in tree.leaf_edge_iter(filter_fn=f)]
        visited_labels = [nd.label for nd in nodes]
        exp_labels = [x for x in self.leaf_sequence if self.node_edge_lengths[x] > 13]
        self.assertEqual(visited_labels, exp_labels)

    ###########################################################################
    ## Special Iterators

    def test_child_iterator_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            expected_children = self.node_children[nd.label]
            children = [ch.label for ch in nd.child_iter()]
            self.assertEqual(children, expected_children)

    def test_child_iterator_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_children = [label for label in self.node_children[nd.label] if self.node_edge_lengths[label] > 13]
            children = [ch.label for ch in nd.child_iter(filter_fn=filter_fn)]
            self.assertEqual(children, expected_children)

    def test_ancestor_iterator_exclusive_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=False)]
            expected_ancestors = self.node_ancestors[nd.label]
            self.assertEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_exclusive_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_ancestors = self.node_ancestors[nd.label]
            expected_ancestors = [nda for nda in expected_ancestors if self.node_edge_lengths[nda] > 13]
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=False, filter_fn=filter_fn)]
            self.assertEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_inclusive_unfiltered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        for nd in anodes:
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=True)]
            expected_ancestors = [nd.label] + self.node_ancestors[nd.label]
            self.assertEqual(ancestors, expected_ancestors)

    def test_ancestor_iterator_inclusive_filtered(self):
        tree, anodes, lnodes, inodes = self.get_tree()
        filter_fn = lambda x: x.edge.length > 13
        for nd in anodes:
            expected_ancestors = [nd.label] + self.node_ancestors[nd.label]
            expected_ancestors = [nda for nda in expected_ancestors if self.node_edge_lengths[nda] > 13]
            ancestors = [ch.label for ch in nd.ancestor_iter(inclusive=True, filter_fn=filter_fn)]
            self.assertEqual(ancestors, expected_ancestors)


if __name__ == "__main__":
    unittest.main()
