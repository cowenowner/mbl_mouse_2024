import spikeinterface.full as si
import probeinterface as pif
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import os

# Define the folder paths
npx1_path = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '01_7_6_23', 'HC101_23242_g0')
npx2_path = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '02_7_7_23', 'HC_2_Neuro2_g0')

meta_filename = os.path.join('E:', 'Dropbox (Dartmouth College)', 'NSB2023', '23242', '02_7_7_23', 'HC_2_Neuro2_g0', 'HC_2_Neuro2_g0_imec0', 'HC_2_Neuro2_g0_t0.imec0.ap.meta')
raw_rec = si.read_spikeglx(npx2_path, stream_name='imec0.ap')
ttls = si.read_spikeglx(npx2_path, stream_name='nidq')

rec1 = si.highpass_filter(raw_rec, freq_min=400.)
bad_channel_ids, channel_labels = si.detect_bad_channels(rec1)
rec2 = rec1.remove_channels(bad_channel_ids)
print('bad_channel_ids', bad_channel_ids)

rec3 = si.phase_shift(rec2)
rec4 = si.common_reference(rec3, operator="median", reference="global")

# here we use static plot using matplotlib backend
fig, axs = plt.subplots(ncols=3, figsize=(20, 10))

si.plot_timeseries(rec1, backend='matplotlib',  clim=(-50, 50), ax=axs[0])
si.plot_timeseries(rec3, backend='matplotlib',  clim=(-50, 50), ax=axs[1])
si.plot_timeseries(rec4, backend='matplotlib',  clim=(-50, 50), ax=axs[2])
for i, label in enumerate(('filter', 'cmr', 'final')):
    axs[i].set_title(label)
    
plt.show()

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