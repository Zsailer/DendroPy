#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Wrappers for interacting with NCBI databases.
"""

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

import urllib
import dendropy
import re

GB_FASTA_DEFLINE_PATTERN = re.compile(r'^gi\|(\d+)\|gb\|([\w\d]+).(\d+)\|(.*)$')

def compose_taxon_label_from_gb_defline(gb_defline,
        num_desc_components=3,
        separator='_',
        gbnum_in_front=True):
    """
    If `gb_defline` matches a GenBank FASTA-format defline structure, then this returns a
    label:

        <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

    So, for example, given the following FASTA label:

        gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

    the corresponding taxon 3-component (default) label will be:

        EU105975_Homo_sapiens_Ache

    If `gb_defline` does *not* match a GenBank FASTA-format defline structure, then the string
    is returned unchanged.
    """
    m = GB_FASTA_DEFLINE_PATTERN.match(gb_defline)
    if m is not None:
        groups = m.groups()
        desc_parts = [s.strip() for s in groups[-1].split() if s]
        if gbnum_in_front:
            label_parts = [groups[1]] + desc_parts[:num_desc_components]
        else:
            label_parts = desc_parts[:num_desc_components] + [groups[1]]
        return separator.join(label_parts)
    else:
        return gb_defline

def relabel_taxa_from_defline(taxon_set,
        num_desc_components=3,
        separator='_',
        gbnum_in_front=True):
    """
    Examines the labels of each `Taxon` object in `taxon_set`, and if
    conforming to a GenBank pattern, translates the labels to a standard
    format:

        <GB-ACCESSION-ID><SEPARATOR><DESC_COMPONENT(1)><SEPARATOR><DESC_COMPONENT(2)>...<DESC_COMPONENT(n)>

    So, for example, given the following FASTA label:

        gi|158931046|gb|EU105975.1| Homo sapiens Ache non-coding region T1584 genomic sequence

    the corresponding taxon 3-component (default) label will be:

        EU105975_Homo_sapiens_Ache

    """
    for taxon in taxon_set:
        taxon.label = compose_taxon_label_from_gb_defline(
                gb_defline=taxon.label,
                num_desc_components=num_desc_components,
                separator=separator,
                gbnum_in_front=gbnum_in_front)
    return taxon_set

class Entrez(object):
    """
    Wraps up all interactions with Entrez.
    Example usage::

        >>> from dendropy.interop import ncbi
        >>> e = ncbi.Entrez(generate_labels=True,
        ... label_gbnum_in_front=False,
        ... sort_taxa_by_label=True)
        >>> d = e.fetch_nucleotide_accession_range(105474, 106046, prefix="EU")

    """

    BASE_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    DATABASES = [
        'pubmed',
        'protein',
        'nucleotide',
        'nuccore',
        'nucgss',
        'nucest',
        'structure',
        'genome',
        'biosystems',
        'books',
        'cancerchromosomes',
        'cdd',
        'gap',
        'dbvar',
        'domains',
        'epigenomics',
        'gene',
        'genomeprj',
        'gensat',
        'geo',
        'gds',
        'homologene',
        'journals',
        'mesh',
        'ncbisearch',
        'nlmcatalog',
        'omia',
        'omim',
        'pepdome',
        'pmc',
        'popset',
        'probe',
        'proteinclusters',
        'pcassay',
        'pccompound',
        'pcsubstance',
        'seqannot',
        'snp',
        'sra',
        'taxonomy',
        'toolkit',
        'toolkitall',
        'unigene',
        'unists',
        'linkoutpubmed',
        'linkoutseq',
        'linkoutother',
        ]

    def __init__(self,
            generate_labels=False,
            label_num_desc_components=3,
            label_separator='_',
            label_gbnum_in_front=True,
            sort_taxa_by_label=False):
        """
        Instantiates a broker that queries NCBI and returns data.  If
        ``generate_labels`` is ``True``, then appropriate labels for sequences
        will be automatically composed for each sequence based on the GenBank
        FASTA defline. ``label_num_desc_components`` specifies the number of
        components from the defline to use. ``label_gbnum_in_front`` specifies
        whether the GenBank accession number should form the beginning
        (``True``) or tail (``False``) end of the label. ``sort_taxa_by_label``
        specifies whether the sequences should be sorted by their final label
        values.
        """
        self.generate_labels = generate_labels
        self.label_num_desc_components = label_num_desc_components
        self.label_separator = label_separator
        self.label_gbnum_in_front = label_gbnum_in_front
        self.sort_taxa_by_label = sort_taxa_by_label

    def _fetch(self, db, ids, rettype):
        """
        Raw fetch. Returns file-like object opened for reading on string
        returned by query.
        """
        if isinstance(ids, str):
            id_list = ids
        else:
            id_list = ",".join([str(i) for i in set(ids)])
        params = {'db': db,
                'id': id_list,
                'rettype': rettype,
                'retmode': 'text'}
        query_url = Entrez.BASE_URL + "/efetch.fcgi?" + urllib.urlencode(params)
        query = urllib.urlopen(query_url)
        return query

    def fetch_nucleotide_accession_ids(self, ids, prefix=None, **kwargs):
        """
        Returns a DnaCharacterMatrix object populated with sequences from the
        Entrez nucleotide database with accession numbers given by `ids` (a
        list of numbers). If `prefix` is given, it is pre-pended to all
        values given in the id list. Any other keyword arguments given are
        passed to the constructor of ``DnaCharacterMatrix``.
        """
        if prefix is not None:
            ids = ["%s%s" % (prefix,i) for i in ids]
        results = self._fetch(db='nucleotide', ids=ids, rettype='fasta')
        d = dendropy.DnaCharacterMatrix.get_from_stream(results, 'fasta', **kwargs)
        if self.generate_labels:
            relabel_taxa_from_defline(d.taxon_set,
                    num_desc_components=self.label_num_desc_components,
                    separator=self.label_separator,
                    gbnum_in_front=self.label_gbnum_in_front)
        if self.sort_taxa_by_label:
            d.taxon_set.sort(key=lambda x: x.label)
        return d

    def fetch_nucleotide_accession_range(self, start, stop, prefix=None, **kwargs):
        """
        Returns a DnaCharacterMatrix object populated with sequences from the
        Entrez nucleotide database with accession numbers between ``start``
        and, up to, but *not* including, ``end`` (i.e., behavior is identical
        to Python's built-in ``range``). If `prefix` is given, then it
        is pre-pended to the ids. Any other keyword arguments given are passed
        to thee constructor of ``DnaCharacterMatrix``.
        """
        ids = range(start, stop)
        return self.fetch_nucleotide_accession_ids(ids=ids, prefix=prefix, **kwargs)

