import os
import re

# for checking : http://www.julesberman.info/cgi-bin/postback.cgi

class Taxonomy(object):
    """This class should ease the retrieval and manipulation of Taxons

    There are many resources to retrieve information about a Taxon.
    For instance, from BioServices, one can use UniProt, Ensembl, or
    EUtils. This is convenient to retrieve a Taxon (see :meth:`fetch_by_name`
    and :meth:`fetch_by_id` that rely on Ensembl). However, you can
    also download a flat file from EBI ftp server, which
    stores a set or records (1.3M at the time of the implementation).

    Note that the Ensembl database does not seem to be as up to date
    as the flat files but entries contain more information.

    for instance taxon 2 is in the flat file but not available through
    the :meth:`fetch_by_id`, which uses ensembl.

    So, you may access to a taxon in 2 different ways getting differnt
    dictionary. However, 3 keys are common (id, parent, scientific_name)

    ::

        >>> t = taxonomy.Taxonomy()
        >>> t.fetch_by_id(9606)   # Get a dictionary from Ensembl
        >>> t.records[9606] # or just try with the get
        >>> t[9606]
        >>> t.get_lineage(9606)

    """
    def __init__(self, verbose=True, online=True):
        """.. rubric:: constructor

        :param offline: if you do not have internet, the connction to Ensembl
            may hang for a while and fail. If so, set **offline** to True

        """
        if online:
            from bioservices import Ensembl
            self.ensembl = Ensembl(verbose=False)
        self.records = {} # empty to start with.
        self.verbose = verbose

    def _load_flat_file(self, overwrite=False):
        """Loads entire flat file from EBI

        Do not overwrite the file by default.
        """
        import ftplib
        output_filename = 'taxonomy.dat'
        self.name = output_filename
        from sequana import sequana_config_path
        self.filename = sequana_config_path + os.sep + self.name
        if os.path.exists(self.filename) and overwrite is False:
            return

        url = 'ftp.ebi.ac.uk' # /pub/databases/taxonomy/'
        self.ftp = ftplib.FTP(url)
        self.ftp.login()
        self.ftp.cwd('pub')
        self.ftp.cwd('databases')
        self.ftp.cwd('taxonomy')

        print('Downloading and saving in %s' % self.filename)
        self.ftp.retrbinary('RETR taxonomy.dat',
                open(self.filename, 'wb').write)

    def load_records(self, overwrite=False):
        """Load a flat file and store records in :attr:`records`

        """
        self._load_flat_file(overwrite=overwrite)
        self.records = {}

        # TODO: check if it exists otherwise, load it ?
        if os.path.exists(self.filename) is False:
            self.load()

        with open(self.filename) as f:
            data = f.read().strip()

        data = data.split("//\n") # the sep is //\n
        self._child_match = re.compile('ID\s+\:\s*(\d+)\s*')
        self._parent_match = re.compile('PARENT ID\s+\:\s*(\d+)\s*')
        self._rank_match = re.compile('RANK\s+\:\s*([^\n]+)\s*')
        self._name_match = re.compile('SCIENTIFIC NAME\s+\:\s*([^\n]+)\s*')

        from easydev import Progress
        pb = Progress(len(data))

        if self.verbose:
            print('Loading all taxon records.')
        for i, record in enumerate(data[0:]):
            # try/except increase comput. time by 5%
            try:
                dico = self._interpret_record(record)
                identifier = int(dico['id'])
                self.records[identifier] = dico
            except Exception as err:
                print(err)
                print('Could not parse the following record '  + \
                      'Please fill bug report on http://github.com/biokit')
                print(record)
            if self.verbose:
                pb.animate(i+1)
        if self.verbose:
            print()

    def _interpret_record(self, record):
        data = {'raw': record}
        # in principle, we should check for the existence of a match
        # but this takes time. All entries must have an ID so no
        # need to check for it. Same for parent and scientific name.
        # Does not save that much time though.
        m = self._child_match.search(record)
        if m: data['id'] = m.group(1)
        m = self._parent_match.search(record)
        if m: data['parent'] = m.group(1)
        m = self._name_match.search(record)
        if m: data['scientific_name'] = m.group(1)
        m = self._rank_match.search(record)
        if m: data['rank'] = m.group(1)

        return data

    def fetch_by_id(self, taxon):
        """Search for a taxon by identifier

        :return; a dictionary.

        ::

            >>> ret = s.search_by_id('10090')
            >>> ret['name']
            'Mus Musculus'

        """
        res = self.ensembl.get_taxonomy_by_id(taxon)
        return res

    def fetch_by_name(self, name):
        """Search a taxon by its name.

        :param str name: name of an organism. SQL cards possible e.g.,
            _ and % characters.
        :return: a list of possible matches. Each item being a dictionary.

        ::

            >>> ret = s.search_by_name('Mus Musculus')
            >>> ret[0]['id']
            10090

        """
        res = self.ensembl.get_taxonomy_by_name(name)
        return res

    def on_web(self, taxon):
        """Open UniProt page for a given taxon"""
        # Should work for python2 and 3
        import webbrowser
        try:
            from urllib.request import urlopen
            from urllib.error import HTTPError, URLError
        except:
            from urllib2 import urlopen, HTTPError, URLError
        try:
            urlopen('http://www.uniprot.org/taxonomy/%s' % taxon)
            webbrowser.open("http://www.uniprot.org/taxonomy/%s" % taxon)
        except HTTPError as err:
            print("Invalid taxon")
        except URLError as err:
            print(err.args)

    def get_lineage(self, taxon):
        """Get lineage of a taxon

        :param int taxon: a known taxon
        :return: list containing the lineage

        """
        # important to reinit the second argument to []
        taxon = int(taxon)
        lineage = self._gen_lineage_and_rank(taxon, [])
        lineage = [x[0] for x in lineage]
        return lineage

    def _gen_lineage_and_rank(self, taxon, lineage_rank=[]):
        # recursively filling the lineage argument
        if len(self.records) == 0:
            self.load_records()

        try:
            record = self.records[taxon]
        except:
            return [('unknown', 'no rank')]

        parent = int(record['parent'])
        if parent not in [0]:
            lineage_rank.append((record['scientific_name'], record['rank']))
            return self._gen_lineage_and_rank(parent, lineage_rank)
        else:
            lineage_rank.reverse()
            return lineage_rank

    def get_lineage_and_rank(self, taxon):
        """Get lineage and rank of a taxon

        :param int taxon:
        :return: a list of tuples. Each tuple is a pair of taxon name/rank
            The list is the lineage for to the input taxon.

        """
        taxon = int(taxon)
        lineage = self._gen_lineage_and_rank(taxon, [])
        return lineage

    def get_children(self, taxon):
        if len(self.records) == 0:
            self.load_records()
        taxon = str(taxon)
        children = [self.records[k] for k in self.records.keys()
                if self.records[k]['parent'] == taxon]
        children = [child['id'] for child in children]
        return children

    def get_family_tree(self, taxon):
        """root is taxon and we return the corresponding tree"""
        # should limit the tree size
        # uniprot flat files has no record about children, so we would
        # need to reconstruct the tree
        tree = {}
        children = self.get_children(taxon)
        if len(children) == 0:
            return tree
        else:
            return [self.get_family_tree(child) for child in children]

    def __getitem__(self, iden):
        if len(self.records) == 0:
            self.load_records()
        return self.records[iden]


class Taxon(object):
    """Some codecs"""
    def __init__(self, taxid):
        self.taxid = taxid


    def to_genbank(self, retmax=10000):
        """Draft: from a TaxID, uses EUtils to retrieve
        the GenBank identifiers

        :Inspiration: https://gist.github.com/fjossinet/5673672
        """
        from bioservices import EUtils
        e = EUtils()
        idlist = e.ESearch(db='nucleotide',
                term='txid%s[Organism:exp]'%self.taxid,
                restart=0, retmax=retmax)['idlist']
        results = e.ESummary(db='nucleotide', id=idlist, retmax=retmax)
        return results


class Lineage(object):
    def __init__(self, lineage):
        self.lineage = lineage[:]

    def __str__(self):
        txt = ""
        for i,this in enumerate(self.lineage):
            N = i
            txt += N*" " + "%s\n" % this
        return txt
