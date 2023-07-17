import spikeinterface.full as si
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import os

# Define the folder paths
npx1_path = '/Users/manishm/Dropbox (Dartmouth College)/NSB2023/23242/01_7_6_23/HC101_23242_g0'
npx2_path = '/Users/manishm/Dropbox (Dartmouth College)/NSB2023/23242/02_7_7_23/HC_2_Neuro2_g0/'

cur_streams = si.get_neo_streams('spikeglx', npx2_path, load_channel_location=True)
raw_rec = si.read_spikeglx(npx2_path, stream_name='imec0.ap')
ttls = si.read_spikeglx(npx2_path, stream_name='nidq')

this_ids = ['nidq#XA1', 'nidq#XA2', 'nidq#XA4']
this_gain = [ttls.get_channel_property(channel_id=x, key='gain_to_uV') for x in this_ids]
this_offset = [ttls.get_channel_property(channel_id=x, key='offset_to_uV') for x in this_ids]

this_signal = ttls.get_traces(channel_ids=[this_ids[2]])*this_gain[2] + this_offset[2]
plt.plot(this_signal)
plt.axhline(np.mean(this_signal) + np.std(this_signal), color='black')
plt.axhline(np.mean(this_signal) - np.std(this_signal), color='black')

plt.show()