"""
This is a script for the base classes of the gene enrichment module adapted from
    https://github.com/DessimozLab/go-handbook/blob/master/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb
"""
import os
import pandas as pd
from matplotlib import image as mpimg, pyplot as plt

import proteomics_tools.gene_enrichment._utils as ut


# II Main Functions
class BaseGO(ut.Utils):
    """Base class for querying GO terms
    Ref: https://github.com/DessimozLab/go-handbook/blob/master/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb"""
    def __init__(self, dict_go=None):
        ut.Utils.__init__(self)
        self.dict_go = dict_go
        # Namespace
        self.dict_go_ns = {"BP": "biological_process",
                           "CC": "cellular_component",
                           "MF": "molecular_function"}
        self.dict_evidence_codes = {"experiment": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"],
                                    "high_throughput_experiment": ["HTP", "HDA", "HMP", "HGI", "HEP"],
                                    "phylogenetic": ["IBA", "IBD", "IKR", "IRD"],
                                    "computation": ["ISS", "ISO", "ISA", "ISM", "IGC", "RCA"],
                                    "author_statement": ["TAS", "NAS"],
                                    "curator_statement": ["IC", "ND"],
                                    "electronic": ["IEA"]}

    # 1. Querying GO
    def _go_terms_to_list(self, go_cat=None, set_go_terms=None, relation="parent"):
        """Convert go_term into nested list that can be used for input of pd dataframe"""
        def rep(term, cat_str):
            return term.replace("level-", "").replace("depth-", "").replace(" [{}]".format(cat_str), "")
        if relation not in ["parent", "child"]:
            raise ValueError("'relation' must be parent or child")
        list_relations = [[rep(x, self.dict_go_ns[go_cat]) for x in str(self.dict_go[term]).split("\t")] +
                          [relation] for term in set_go_terms]
        return list_relations

    def count_go_term(self, query_term="growth"):
        """Count occurrence of go term"""
        count = 0
        for go_term in self.dict_go.values():
            if query_term in go_term.name:
                count += 1
        return count

    # 2. Parent
    def common_parents(self, go_ids=None):
        """ Find the common ancestors in the GO tree of the list of terms in the input."""
        # Find candidates from first
        rec = self.dict_go[go_ids[0]]
        candidates = rec.get_all_parents()
        candidates.update({go_ids[0]})

        # Find intersection with second to nth term
        for go_id in go_ids[1:]:
            rec = self.dict_go[go_id]
            parents = rec.get_all_parents()
            parents.update({go_id})
            # Find the intersection with the candidates, and update.
            candidates.intersection_update(parents)
        return candidates

    def nearest_common_parent(self, go_ids=None):
        """This function gets the nearest/deepest common ancestor using the above function.
        Only returns single most specific - assumes unique exists. In WordNet called least common subsumer (lcs)"""
        # Take the element at maximum depth
        ncp = max(self.common_parents(go_ids=go_ids), key=lambda t: self.dict_go[t].depth)
        return ncp

    # 3. GO hierarchy
    def get_parents(self, go_id=None):
        """Get all go parents"""
        rec = self.dict_go[go_id]
        set_parents = rec.get_all_parents()
        return set_parents

    def get_children(self, go_id=None):
        """Get all go children"""
        rec = self.dict_go[go_id]
        set_parents = rec.get_all_children()
        return set_parents

    def get_go_hierarchy(self, go_id=None, go_cat="BP"):
        """Find parents for term1 from go_term_set"""
        rec = self.dict_go[go_id]
        set_parents = rec.get_all_parents()
        set_children = rec.get_all_children()
        list_go_parents = self._go_terms_to_list(go_cat=go_cat,
                                                 set_go_terms=set_parents,
                                                 relation="parent")
        list_go_children = self._go_terms_to_list(go_cat=go_cat,
                                                  set_go_terms=set_children,
                                                  relation="child")
        list_go_terms = list_go_parents + list_go_children
        columns = ["GO_term", "level", "depth", self.dict_go_ns[go_cat], "relation"]
        df_hierarchy = pd.DataFrame(list_go_terms, columns=columns)
        df_hierarchy.sort_values(by="depth", inplace=True)
        df_hierarchy.reset_index(inplace=True, drop=True)
        return df_hierarchy

    def plot_lineage(self, go_id=None, save=False, show=True):
        """Plot GO hierarchy (lineage of parents and children) for given go_id"""
        go_term = self.dict_go[go_id]
        fig_lineage = ut.FOLDER_DATA + go_id + '-lineage.png'
        self.dict_go.draw_lineage([go_term], lineage_img=fig_lineage)
        img = mpimg.imread(fig_lineage)
        imgplot = plt.imshow(img)
        plt.axis("off")
        title = "GO hierarchy for {}".format(go_id)
        plt.title(title, fontsize=10)
        plt.tight_layout()
        if not save:
            os.remove(fig_lineage)
        if show:
            plt.show()
        else:
            return plt.gcf()

    # 4. Annotations
    @staticmethod
    def _check_evidence_output(out):
        """Check if evidence output in list"""
        list_out = ["evidence_code", "evidence_category"]
        if out not in list_out:
            raise ValueError("{} not in {}".format(out, list_out))

    @staticmethod
    def _plot_evidence_code(title=None):
        """Pie plot evidence code"""
        plt.ylabel("")
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def evidence_code(self, df_id_entry=None, out="evidence_category", show=True, organism="human"):
        """Get overview of evidence """
        self._check_evidence_output(out)
        counts = df_id_entry["Evidence"].value_counts()  # http://geneontology.org/docs/guide-go-evidence-codes/
        list_evidence = []
        for code, count in counts.iteritems():
            for cat in self.dict_evidence_codes:
                if code in self.dict_evidence_codes[cat]:
                    list_evidence.append([cat, code, count])
        df_evidence = pd.DataFrame(data=list_evidence, columns=["category", "code", "count"]).sort_values(by="category")
        df_cat = df_evidence.groupby("category").sum()
        if out == "evidence_code":
            if show:
                counts.plot.pie(autopct='%1.1f%%', labeldistance=1.1, explode=[0.05] * len(counts))
                plt.legend(bbox_to_anchor=(1.2, 0.75), fontsize=8)
                title = "Count of evidence code for {}".format(organism)
                self._plot_evidence_code(title=title)
            return counts
        else:
            if show:
                df_cat.plot(kind="pie", y="count", legend=False, autopct='%1.1f%%')
                title = "Count of evidence category for {}".format(organism)
                self._plot_evidence_code(title=title)
            return df_cat
