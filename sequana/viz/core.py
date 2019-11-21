"""Core function for the plotting tools"""
import pandas as pd


__all__ = ["VizInput2D"]


class VizInputSquare(object):
    def __init__(self, x, verbose=False):
        self.verbose = verbose
        self.df = pd.DataFrame(x)


class VizInput2D(object):
    def __init__(self, x, y=None, verbose=False):
        self.verbose = verbose

        self.xy_names = ['x', 'y']
        if isinstance(x, pd.DataFrame) is True:
            self.df = x.copy()
            columns = list(self.df.columns)
            columns[0] = 'x'
            columns[1] = 'y'
            self.xy_names = self.df.columns[0:2]
            self.df.columns = columns
        elif y is None:
            # could be a list of lists, a pandas-compatible dictionary
            self.df = pd.DataFrame(x)
            if self.df.shape[1] != 2:
                if self.df.shape[0] == 2:
                    print("warning transposing data")
                    self.df = self.df.transpose()
        elif x is not None and y is not None:
            self.df = pd.DataFrame({'x':x, 'y':y})
        else:
            raise ValueError('not a dataframe or list of items or dictionary.')




