import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from pathlib import Path
import datetime


def pull_data(fname:str, key = 'box'):
    """
    this function extracts the data from the data files
    produced by medpc in the stripped with variable
    identification format. the result is a dictionary
    where the keys are the boxes and the corresponding
    values are lists of dictionaries containing all data
    associated with that box stored in the given file.
    each of these sub-dictionaries maps the medpc variable
    names (the letters A-B) to the associated values.
    """

    more=True
    all_data={}
    var_names=list(map(chr, range(97, 123)))
    with open(fname) as f: d=f.read().split('\n')
    while more:
        st_ind = d.index('.3')
        st_date = '/'.join(d[:3])
        end_date = '/'.join(d[3:6])
        data = { 'subject': d[6],
                 'experiment': d[7],
                 'group': d[8],
                 'box': int(d[st_ind-9])}
        data['start'] = datetime.datetime.strptime(f"{st_date} {':'.join(d[st_ind-8:st_ind-5])}", 
                                                   "%m/%d/%y %H:%M:%S")
        data['end'] = datetime.datetime.strptime(f"{end_date} {':'.join(d[st_ind-5:st_ind-2])}", 
                                                 "%m/%d/%y %H:%M:%S") 
        d = d[st_ind+1:]
        arr_lens = list(map(int,d[:26]))
        del d[:26]
        data.update({v: np.empty((1, arr_lens[i])) for i, v in enumerate(var_names)})
        for i,length in enumerate(arr_lens):
             data[var_names[i]][:]=d[:length]
             del d[:length]
        if key == 'box':
            if data['box'] in all_data:
                all_data[data['box']].append(data)
            else:
                all_data[data['box']] = [data]
        if key == 'subject':
            if data['subject'] in all_data:
                all_data[data['subject']].append(data)
            else:
                all_data[data['subject']] = [data]
        more=True if len(d)>2 else False
    return all_data

def pull_multiple(fld:str, key = 'box', use_most_recent = True , log = False):
    """
    given a folder with data files 
    """
    data={}
    for i in os.listdir(fld):
        if i == '.DS_Store': continue
        else:
            file = os.path.join(fld,i)
            d = pull_data(file, key)
            for k in d:
                if len(d[k])>1 and log: 
                    print(f"WARNING: file {file} has multiple sessions stored for {key} {k}") 
                if use_most_recent:
                    d_to_append = d[k][-1]
                else:
                    d_to_append = d[k] 
                if k in data:
                    if log:
                        print(f"WARNING: multiple files have data for {key} {k}") 
                    data[k].append(d_to_append)
                else:
                    data[k] = [d_to_append]
    return data

def write_data(fname , d, bx):
    with open(fname,'w') as f:
        f.write(d['start'].strftime("%m \n%d \n%y \n"))
        f.write(d['end'].strftime("%m \n%d \n%y \n"))
        f.write(f'{bx} \n')
        f.write(d['start'].strftime("%H \n%M \n%S \n"))
        f.write(d['end'].strftime("%H \n%M \n%S \n"))
        del d['start']
        del d['end']
        for i in range(2): f.write('NaN \n')
        f.write('.3\n')
        for i in d: f.write(f'{d[i].size}\n')
        for i in d: 
            for j in d[i].flatten(): 
                f.write(f'{j}\n')

if __name__ == '__main__':
    """
    TODO
    """
    pass