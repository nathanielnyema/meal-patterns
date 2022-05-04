import pandas as pd
import numpy as np
import datetime
from datetime import timedelta
import csv

def compute_inter(x):
    """
    when we exclude data points the inter event intervals 
    automatically generated are no longer correct so we
    recalculate them here
    """
    diff = (x.start_dts.iloc[1:].reset_index(drop=True) - \
            x.end_dts.iloc[:-1].reset_index(drop=True)).dt.total_seconds()/60
    diff.index = x.iloc[1:].index.get_level_values('Event')
    return diff

def load_labmaster(fpath, excl):
    """
    load the labmaster data
    """
    d=[]
    with open(fpath) as f:
        reader=csv.reader(f,delimiter=';')
        for row in reader:
            d.append(row)
    header_idx = d.index(list(filter(lambda x:'Date' in x, d))[0])
    df = pd.DataFrame(d[header_idx+1:], columns=d[header_idx])
    df['start_dts'] = pd.to_datetime(df.Date + ' ' + df.Time)
    df.set_index('start_dts', inplace = True)
    df = df.drop(columns = ['Date','Time'])
    sp = df.columns.str.split(': ')
    boxes = sp.map(lambda x: x[0]).str.extract('Box(\d)')[0].map(int)
    ev = sp.map(lambda x: x[-1]).tolist()
    df.columns  = pd.MultiIndex.from_tuples(list(zip(boxes, ev)), names = ['Cage', 'feeder'])
    df = df.astype(float)
    feeders = {}
    for f in df.columns.get_level_values('feeder').unique():
        feeder = df.stack('Cage')[f]
        feeder.name = 'In_g'
        feeder = feeder.reset_index()
        feeder['end_dts'] = feeder.start_dts + timedelta(seconds = 10)
        feeders[f] = feeder.set_index(['Cage', 'start_dts'])
    df = pd.concat(feeders, names = ['feeders'])
    return df

def load_data(fpath, sheets, labmaster, excl):
    """
    general function for loading data and prepping
    for downstream analyses.

    Parameters
    ----------
    fpath: str
        the path to the data file
    sheets: list
        a list of the names of the sheets with relevant data
    labmaster: bool
        whether or not this data was collected with the labmaster
    
    Returns
    -------
    df: pd.DataFrame
    """
    if labmaster:
        df = load_labmaster(fpath, excl)
    else:
        df = pd.read_excel(fpath , engine = 'openpyxl', sheet_name = sheets)
        df = pd.concat(df, names = ['feeders'])
        df['start_dts'] = pd.to_datetime(df.StartDate + ' ' +  df.StartTime)
        df['end_dts'] = pd.to_datetime(df.EndDate + ' ' +  df.EndTime)
        df = df.reset_index().set_index(['feeders','Cage', 'start_dts'])
   
    def label_evs(x):
        x['Event'] = np.arange(1,len(x) + 1)
        return x

    df = df.loc[df.In_g>=excl].sort_index()
    df = df.groupby(['feeders', 'Cage']).apply(label_evs)
    df = df.reset_index().set_index(['feeders','Cage', 'Event'])
    diff = df.groupby(['feeders','Cage']).apply(compute_inter)
    df['InterIn_min'] = diff
    df = df[['In_g','InterIn_min','start_dts', 'end_dts']]
    return df

def get_bouts(df, inter_thresh):

    """
    identify bouts and compute relevant statistics on them
    we will use the pandas groupby function to do this so
    we also define a few nested functions here to be called 
    on sub dataframes corresponding to a given cage or a given
    bout for a given cage. our approach will be as follows

      1. initialize a column with NaNs
      2. threshold the inter-event intervals to identify the starts
         of meals and label them from 0-n bouts for each cage
      3. forward fill the NaNs then compute relevant stats for all bouts
    
    Parameters
    ----------
    df : pd.DataFrame
        the dataframe with the event data. this dataframe
        must have the following structure
        
        * the index must be multi-index with levels ['Cage', 'Event']
        * there must be a column start_dts with the start datetimes of the events
        * there must be a column end_dts with the start datetimes of the events
        * there must be a column InterIn_min with the time interval between events in minutes
        * there must be a column In_g with the amount consumed in grams during the event
    
    inter_thresh: float, int
        the threshold for defining a new meal.

    Returns
    -------
    
    bout_stats: pd.DataFrame

    """

    def find_bouts(x):
        """
        use the InterIn_min field of the dataframe to
        identify bouts given our inter-bout threshold 
        """
        bout_starts = x.InterIn_min>inter_thresh
        try:
            bout_starts.loc[pd.IndexSlice[:,1]] = True
        except KeyError:
            pass
        x.loc[bout_starts,'bout_n'] = np.arange(bout_starts.sum())
        return x
    
    def get_bout_stats(x):
        """
        get the duration and size of a given bout
        """
        st = x.start_dts.iloc[0]
        en =  x.end_dts.iloc[-1]
        diff = (en - st).total_seconds()
        df = {
            'start' : st,
            'end' : en,
            'dur_s' : diff,
            'size' : x.In_g.sum()
        }
        return pd.Series(df)
    
    def get_imis(x):
        """
        get the inter-meal intervals. we'll use a convention of
        assigning IMIs to the meal that ends the inter-meal interval
        """
        imis = (x.start.iloc[1:].reset_index(drop=True) \
                - x.end.iloc[:-1].reset_index(drop=True)).dt.total_seconds()
        imis.index = x.iloc[1:].index
        return imis.droplevel(0)
    
    df['bout_n'] = np.ones(df.shape[0])*np.nan 
    df = df.groupby('Cage').apply(find_bouts).fillna(method = 'ffill')
    bout_stats = df.groupby(['Cage','bout_n']).apply(get_bout_stats)
    bout_stats['imis_s'] = bout_stats.groupby('Cage').apply(get_imis)
    bout_stats.loc[bout_stats.imis_s<0,'imis_s'] = 0. # force negative imis to 0. this happens when theres only 1 meal
    return bout_stats.fillna(0), df

def bin_stats(bout_stats, bins):
    """
    given the output of get_bouts, bin the bout statistics
    using pivot table averaging

    Parameters:
    ----------
    bout_stats: pd.DataFrame
        the output from get_bouts
    bins_start: datetime.datetime
        the timestamp of the start of the first bin
    bins_end: datetime.datetime
        the timestamp of the end of the last bin
    binsize_hr: float, int
        the bin size in hours for binning the data
    """
    
    bout_stats['bins'] = pd.cut(bout_stats.start, bins).apply(lambda x: x.left).astype(np.datetime64)
    binned_stats = bout_stats.pivot_table(index = ['Cage', 'bins'], 
                                          values = ['dur_s', 'size', 'imis_s'])
    binned_stats['meal_n'] = bout_stats.groupby(['Cage', 'bins']).apply(len)
    return binned_stats


def run_analysis(df_in, inter_thresh, binsize_hr, start):
    """
    Parameters:
    ----------
    df_in: pd.DataFrame
        the output from load_data
    inter_thresh: float, int
        the threshold for defining a new meal.
    excl: float, int
        exclusion criteria for events. i.e. events where less than this amount
        was consumed will be excluded
    binsize_hr: float, int
        the bin size in hours for binning the data
    start: datetime.datetime
        the timestamp of for the first bin. if not specified we automatically
        use the timestamp of the earliest event as the start time
    
    Returns
    -------
    bout_stats: pd.DataFrame
    binned_stats: pd.DataFrame
    """

    bout_stats = {}
    binned_stats = {}
    for s in df_in.index.get_level_values('feeders').unique():
        df = df_in.loc[s].copy()
        bout_stats[s], _ = get_bouts(df, inter_thresh)
        end = df.end_dts.max()
        bins = pd.date_range(start, end + timedelta(hours = binsize_hr), 
                             freq = f'{binsize_hr}H')
        binned_stats[s] = bin_stats(bout_stats[s], bins)
    
    #combine the dataframes from the feeder events
    bout_stats = pd.concat(bout_stats, names = ['feeder'])
    binned_stats = pd.concat(binned_stats, names = ['feeder'])
    bins = bins[:-1].to_numpy()
    feeders = binned_stats.index.get_level_values("feeder").unique().to_numpy()
    cages = binned_stats.index.get_level_values("Cage").unique().to_numpy()
    new_idx = list(zip(feeders.repeat(bins.size*cages.size), 
                    np.tile(cages.repeat(bins.size), feeders.size), 
                    np.tile(bins, feeders.size*cages.size)))
    binned_stats = binned_stats.reindex(pd.MultiIndex.from_tuples(new_idx, names = ['feeder', 'Cage', 'bins']))
    return bout_stats, binned_stats
    