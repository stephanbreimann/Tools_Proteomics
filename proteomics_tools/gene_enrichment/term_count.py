"""
This is a script for implementing a Term Count class for semantic analysis
"""
from collections import Counter


# II Main Functions
class TermCounts:
    """TermCounts counts the term counts for each term
    In: a) dict_go: go (complete go library)
        b) dict_org_go: organism specfic version of go (annotations)
    Q:  a) https://github.com/DessimozLab/go-handbook/blob/master/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb"""

    def __init__(self, dict_go=None, dict_org_go=None):
        # Backup
        self._dict_go = dict_go
        # Initialise the counters
        self._counts = Counter()
        self._aspect_counts = Counter()
        # Fill the counters
        self._count_terms(dict_go, dict_org_go)

    def _count_terms(self, dict_go, dict_org_go):
        """Fills in the counts and overall aspect counts."""
        for x in dict_org_go:
            try:
                # Extract term information
                go_id = dict_org_go[x]['GO_ID']
                namespace = dict_go[go_id].namespace
                self._counts[go_id] += 1
                rec = dict_go[go_id]
                parents = rec.get_all_parents()
                for p in parents:
                    self._counts[p] += 1
                self._aspect_counts[namespace] += 1
            except KeyError:    # GO id not in given go db
                pass

    def get_count(self, go_id=None):
        """Returns the count of that GO term observed in the annotations"""
        return self._counts[go_id]

    def _get_total_count(self, namespace):
        """Gets the total count that's been precomputed"""
        return self._aspect_counts[namespace]

    def get_term_freq(self, go_id=None):
        """Returns the frequency at which a particular GO term has been observed in the annotations"""
        try:
            namespace = self._dict_go[go_id].namespace  # e.g. biological process
            freq = float(self.get_count(go_id)) / float(self._get_total_count(namespace))
        except ZeroDivisionError:
            freq = 0
        return freq
