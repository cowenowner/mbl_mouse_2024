#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Plot_graph(data):
    import networkx as nx
    import pandas as pd
    data = {'Source': ['1', '2','3'],
    'Target': ['2', '3','1'],
    'Type': ['Directed', 'Directed','Directed'],
    'weight': [3, 1, -3]}
    df = pd.DataFrame(data)
    G = nx.from_pandas_edgelist(df,source = 'Source', target = 'Target', edge_attr = 'weight')
    nx.draw_circular(G)
    nx.draw(G)
