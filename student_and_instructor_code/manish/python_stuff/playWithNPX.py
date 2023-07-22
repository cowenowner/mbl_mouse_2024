import spikeinterface.full as si
import numpy as np
import matplotlib.pyplot as plt
import probeinterface as pif
import scipy.io as scio
import os

# Define the folder paths
npx1_path = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '01_7_6_23', 'HC101_23242_g0')
npx2_path = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '02_7_7_23', 'HC_2_Neuro2_g0')

meta_filename = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '02_7_7_23', 'HC_2_Neuro2_g0', 'HC_2_Neuro2_g0_imec0', 'HC_2_Neuro2_g0_t0.imec0.ap.meta')
probe = pif.read_spikeglx(meta_filename)

# cur_streams = si.get_neo_streams('spikeglx', npx2_path)
# raw_rec = si.read_spikeglx(npx2_path, stream_name='imec0.ap')

# this_probe = raw_rec.get_probe()

# ttls = si.read_spikeglx(npx2_path, stream_name='nidq')

# this_ids = ['nidq#XA1', 'nidq#XA2', 'nidq#XA4']
# this_gain = [ttls.get_channel_property(channel_id=x, key='gain_to_uV') for x in this_ids]
# this_offset = [ttls.get_channel_property(channel_id=x, key='offset_to_uV') for x in this_ids]

# this_signal = ttls.get_traces(channel_ids=[this_ids[2]])*this_gain[2] + this_offset[2]
# plt.plot(this_signal)
# plt.axhline(np.mean(this_signal) + np.std(this_signal), color='black')
# plt.axhline(np.mean(this_signal) - np.std(this_signal), color='black')

# plt.show()