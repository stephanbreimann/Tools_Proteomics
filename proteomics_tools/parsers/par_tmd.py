"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np


# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions
def reverse_string(string):
    """Reverse string"""
    return string[::-1]


def retrieve_string_starting_at_end(seq, start=None, end=None):
    """Reverse_string_start_end"""
    reversed_seq = reverse_string(seq)
    reversed_seq_part = reversed_seq[start:end]
    seq = reverse_string(reversed_seq_part)
    return seq


def standardize_col(df, reset_index=True):
    """Standardize columns uniprot_id and gen_name"""
    if reset_index is True:
        df.reset_index(drop=True, inplace=True)
    df.columns = [x.lower() for x in list(df.columns)]
    dict_name_list = {'entry': ['tmd_id', 'id', 'uniprot_id'],
                      'name': ['gen_name', 'entry_name', "entry name"],
                      'mean': ['total_mean'],
                      'sequence': ['seq'],
                      'start': ['start_region', 'start_tmd'],
                      'stop': ['stop_region', 'stop_tmd'],
                      'subcellular location': ["subcellular location [cc]"]}
    for key in dict_name_list:
        for name in dict_name_list[key]:
            if name in list(df.columns):
                df.rename(columns={name: key}, inplace=True)
    return df


def get_dict_part_seq(tmd=None, jmd_n=None, jmd_c=None, ext_len=4):
    """Get dictionary for part to sequence
    In: a) tmd: sequence of tmd
        b) jmd_n: sequence of jmd_n
        c) jmd_c: sequence of jmd_c
        d) ext_len: length of extending part (starting from C and N terminal part of TMD)"""
    tmd_n = tmd[0:round(len(tmd) / 2)]
    tmd_c = tmd[round(len(tmd) / 2):]
    ext_n = retrieve_string_starting_at_end(jmd_n, start=0, end=ext_len)  # helix_stop motif for TMDs
    ext_c = jmd_c[0:ext_len]  # anchor for TMDs
    tmd_e = ext_n + tmd + ext_c
    part_seq_dict = {'tmd': tmd, 'tmd_e': tmd_e,
                     'tmd_n': tmd_n, 'tmd_c': tmd_c,
                     'jmd_n': jmd_n, 'jmd_c': jmd_c,
                     'tmd_jmd': jmd_n + tmd + jmd_c,
                     'jmd': jmd_n + jmd_c,
                     'jmd_n_tmd_n': jmd_n + tmd_n, 'tmd_c_jmd_c': tmd_c + jmd_c,
                     'anchor': ext_c, 'helix_stop': ext_n,
                     'tmd_c_anchor': tmd_c + ext_c, 'helix_stop_tmd_n': ext_n + tmd_n}

    return part_seq_dict


class TmdIn:
    """Interface class for TMD classes.
    In: a) entry: id of tmd
        b1) tmd_start: start position within the sequence of the tmd
        b2) tmd_stop: stop position within the sequence of the tmd
        b3) sequence: sequence of tmd
        c) jmd_len: length of jmd"""
    def __init__(self, entry, tmd_start, tmd_stop, seq, jmd_len, ext_len):
        self._entry = entry
        self._tmd_start = tmd_start
        self._tmd_stop = tmd_stop
        self._seq = seq
        self._jmd_len = jmd_len
        self._jmd_n_length = self._tmd_start - 1
        self._jmd_c_length = len(self._seq) - self._tmd_stop
        if jmd_len <= 4:
            raise ValueError("Minimum length for 'jmd_len' is 5")
        if ext_len >= jmd_len:
            raise ValueError("'ext_len' must be smaller than jmd_len")
        self._ext_len = ext_len


class ParserBase:
    """Base class for all parsers"""
    def __init__(self, jmd_len=20):
        self.jmd_len = jmd_len

    def _get_tmd_jmd_by_split(self, seq):
        """Split sequence in jmd_n, tmd and jmd_c"""
        if len(seq)/3 < self.jmd_len:
            n = int(len(seq)/3) - 1
        else:
            n = self.jmd_len
        jmd_n = seq[0:n]
        tmd = seq[n:-n]
        jmd_c = seq[-n:]
        return jmd_n, tmd, jmd_c


# II Main Functions
# 1. TMD class
class Tmd(TmdIn):
    """Class of Transmembrane Domains (tmd) with their juxta & transmembrane domain (JMD, TMD).
    Part to sequence can be constructed for a specific tmd."""
    def __init__(self, **kwargs):
        TmdIn.__init__(self, **kwargs)
        self.dict_part_seq = self._dict_part_seq()

    def _tmd(self):
        """Get sequence of original tmd (regarding uniprot annotation)"""
        start = self._tmd_start - 1
        stop = self._tmd_stop
        return self._seq[start:stop]

    def _jmd_n(self):
        """Get N-terminal JMD"""
        if self._jmd_n_length >= self._jmd_len:
            start = self._tmd_start - (self._jmd_len + 1)
            stop = self._tmd_start - 1
            jmd_n = self._seq[start:stop]
        else:
            start = 0
            stop = self._tmd_start - 1
            part = self._seq[start:stop]
            jmd_n = 'X' * (self._jmd_len - len(part)) + part  # Add X for missing AA in JMD
        return jmd_n

    def _jmd_c(self):
        """Get C-terminal juxtamembraneregion"""
        if self._jmd_c_length >= self._jmd_len:
            start = self._tmd_stop
            stop = self._tmd_stop + self._jmd_len
            jmd_c = self._seq[start:stop]
        else:
            start = self._tmd_stop
            stop = self._tmd_stop + self._jmd_c_length
            part = self._seq[start:stop]
            jmd_c = part + 'X' * (self._jmd_len - len(part))  # Add X for missing AA in JMD
        return jmd_c

    def _dict_part_seq(self):
        """Save dictionary of parts to their sequence when initiating object"""
        dict_part_seq = get_dict_part_seq(tmd=self._tmd(), jmd_n=self._jmd_n(),
                                          jmd_c=self._jmd_c(), ext_len=self._ext_len)
        return dict_part_seq


# 2. Adder class
class Adder:
    """Class for converting tmd information from tmd dataframe into unified form. Used in TmdParser class.
    In: a) df: data frame with tmd information (e.g. sequecne, start, stop)
        b) jmd_len: length of jmd
        c) art_mut: internal boolean just set True if artificial or mutated tmds are added
    !:  a) Output of functions modified data frame"""

    def __init__(self, df=None, jmd_len=12, art_mut=False, ext_len=4):
        self.df = standardize_col(df)
        self.col_names = [col.lower() for col in list(df)]
        if jmd_len is None:
            raise ValueError("'jmd_len' should not be None")
        self.jmd_len = jmd_len
        self._ext_len = ext_len
        self._art_mut = art_mut

    def _col_split(self, col_name='transmembrane'):
        """Split col_split to get start and end position of tmd
        In: a) col_name: name of column which contains information of stop and start (normally, not to change)"""
        pd.options.mode.chained_assignment = None # disable warning of copy of slice
        if 'start' in self.col_names and 'stop' in self.col_names:
            pass
        # Get start and stop from transmembrane column
        elif col_name in self.col_names:
            self.df['start'] = [int(i[-2]) for i in self.df[col_name].str.split('[\s+\-\,\;]', n=2)]
            self.df['stop'] = [int(i[-1]) for i in self.df[col_name].str.split('[\s+\-\,\;]', n=2)]
        # Get start and stop using jmd_len (if no transmembrane column but jmd_c and jmd_n are given)
        elif 'jmd_n' in self.col_names and 'jmd_c' in self.col_names:
            self.df['start'] = [self.jmd_len + 1] * len(self.df)
            self.df['stop'] = [self.jmd_len + len(tmd) for tmd in self.df['tmd']]
        # Get start and stop using jmd_len (if artificial sequence is given)
        elif self._art_mut is True:
            self.df['start'] = [self.jmd_len + 1] * len(self.df)
            self.df['stop'] = [len(seq) - self.jmd_len for seq in self.df['sequence']]
        else:
            raise Exception("'{}' is wrong col_name for split (default is 'transmembrane')".format(col_name))
        return self.df

    def _entry(self):
        """Set id name to entry"""
        if 'entry' not in self.col_names:
            self.df.insert(loc=1, column='entry', value= np.nan * len(self.df))
            print("No 'entry' in column names")
        return self.df

    def _add_name_class(self, class_name='NaN', replace=False):
        """Add names & class to df
        In: a) class_name: name of class to insert
            b) replace: boolean to choice if class name should be replaced"""

        def add_names(df=None, col_names=None):
            """Add name to df"""
            if 'name' not in col_names:
                df.insert(loc=1, column='name', value=['NaN'] * len(df))
                print("No 'name' in column names")
            else:
                try:
                    with pd.option_context('mode.chained_assignment', None):
                        #df['name'] = [name.replace('_HUMAN', '') for name in list(df['name'])]
                        df['name'] = [name for name in list(df['name'])]
                except:
                    pass
            return df

        def add_class(df=None, col_names=None, class_name='NaN', replace=True):
            """Add class to df"""
            if 'class' not in col_names:
                df.insert(loc=2, column='class', value=[class_name] * len(df))
            if class_name is not None and replace is True:
                df['class'] = [class_name] * len(df)
            return df

        self.df = add_names(df=self.df, col_names=self.col_names)
        self.df = add_class(df=self.df, col_names=self.col_names,
                            class_name=class_name,
                            replace=replace)
        return self.df

    def _adjust_jmd(self):
        """Adjust length of jmd and insert 'X' if too short"""
        jmd_n_list = []
        for jmd_n in self.df['jmd_n']:
            dif_len = abs(self.jmd_len - len(jmd_n))
            if len(jmd_n) < self.jmd_len:
                jmd_n_list.append('X' * dif_len + jmd_n)
            elif len(jmd_n) > self.jmd_len:
                jmd_n_list.append(jmd_n[dif_len:])
            else:
                jmd_n_list.append(jmd_n)
        self.df['jmd_n'] = jmd_n_list
        jmd_c_list = []
        for jmd_c in self.df['jmd_c']:
            dif_len = abs(self.jmd_len - len(jmd_c))
            if len(jmd_c) < self.jmd_len:
                jmd_c_list.append(jmd_c + 'X' * dif_len)
            elif len(jmd_c) > self.jmd_len:
                jmd_c_list.append(jmd_c[0:self.jmd_len])
            else:
                jmd_c_list.append(jmd_c)
        self.df['jmd_c'] = jmd_c_list

    def _add_tmd_jmd(self):
        """Get jmd_n, tmd, and jmd_n using TMD object"""
        jmd_n_tmd_jmd_c_array = []
        for i, row in self.df.iterrows():
            dict_part_seq = Tmd(entry=row['entry'], seq=row['sequence'],
                                tmd_start=row['start'], tmd_stop=row['stop'],
                                ext_len=self._ext_len,
                                jmd_len=self.jmd_len).dict_part_seq
            jmd_n_tmd_jmd_c_array.append([dict_part_seq[part] for part in ['jmd_n', 'tmd', 'jmd_c']])
        df = pd.DataFrame(jmd_n_tmd_jmd_c_array, columns=['jmd_n', 'tmd', 'jmd_c'])
        self.df = pd.concat([self.df, df], axis=1, sort=False)

    def _add_seq(self):
        """Get sequence joining jmd_n + tmd + jmd_c"""
        self.df['sequence'] = self.df['jmd_n'] + self.df['tmd'] + self.df['jmd_c']

    def unified(self, class_name='NaN', replace=False, col_split='transmembrane'):
        """Main function for uniforming data of tmd dataframe
        In: a) class_name: class name of modified data
            b) replace: boolean for changing class_name
            c) col_name: name of column which contains information of stop and start (normally, not to change)"""
        self.df.columns = self.col_names     # lowercase col names
        # Check entry
        self.df = self._entry()
        # Add name & class column
        self.df = self._add_name_class(class_name=class_name, replace=replace)
        # Split if transmembrane format
        self.df = self._col_split(col_name=col_split)
        # Add seq or tmd and jmd
        if 'jmd_n' in self.col_names:
            self._adjust_jmd()
        if 'jmd_n' not in self.col_names:
            self._add_tmd_jmd()
        elif 'sequence' not in self.col_names:
            self._add_seq()
        return self.df

    @staticmethod
    def drop_superfluous_col(df):
        """Drop superfluous columns at end of modification and order dataframe"""
        col_names = ['entry', 'name', 'protein names', 'class', 'sequence', 'start', 'stop', 'jmd_n', 'tmd', 'jmd_c']
        col_names_old = [col.lower() for col in list(df)]
        col_drop = [name for name in col_names_old if name not in col_names]
        col_keep = [name for name in col_names_old if name not in col_drop]
        df = df.drop(columns=col_drop)
        df = df[col_keep]
        return df


# 3. TMD Parser
class TmdParser:
    """Parser class to read tmd data frames using Adding class for unifying data.
    UniprotData class is used for read_uniprot_ids function to retrieve Uniprot and
    Mutater class is used to read artificial data or data to mutated with given mutations.
    In: a) jmd_len: length of jmd"""

    def __init__(self, jmd_len=12, ext_len=4):
        self.ext_len = ext_len
        self.jmd_len = jmd_len

    def read_df(self, df=None, class_name=None, drop_col=False):
        """Read data frame to converting it as db input
        In: a) df: data frame with tmd info
            b) class_name: name of class of parsed data"""
        df_unified = Adder(df=df.copy(), jmd_len=self.jmd_len, ext_len=self.ext_len).unified(class_name=class_name,
                                                                                             replace=True)
        if drop_col:
            df_unified = Adder.drop_superfluous_col(df_unified)
        return df_unified

    def read_df_list(self, df_list=None, class_name_list=None):
        """Read list of data frames and list of names for converting them as db input
        In: a) df_list: list of data frames
            b) class_name_list: list of class names for data frames (in same order as df_list)"""
        # Get list of unified data frames
        df_list_unified = [Adder(df=df.copy(), jmd_len=self.jmd_len, ext_len=self.ext_len).
                               unified(class_name=class_name_list[i]) for i, df in enumerate(df_list)]
        # Merge unified data
        df_merged = pd.concat(df_list_unified, axis=0, ignore_index=True, sort=False)
        df_unified = Adder.drop_superfluous_col(df_merged)
        return df_unified

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()

