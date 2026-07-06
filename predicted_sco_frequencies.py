import numpy as np
from scipy.io import loadmat
import pandas as pd
from pandas import read_csv, read_table

# Load data, standardise to 100s at dt = 0.01
theta_human = loadmat('./dependencies/theta_human')['Theta_Down'].T[:, :10000]
theta_macaque = loadmat('./dependencies/theta_macaque')['Theta'].T[:, :-1:10]
theta_aal = loadmat('./dependencies/theta_AAL')['Theta_ds'].T[:,:-1]

# Human network parcellation (Schaefer 1000)
labels_n = [name[11:] for name in read_csv('./dependencies/LUT_schaefer.csv')['ROI Name']]
networks_human = {
    name: [i for i, lbl in enumerate(labels_n) if name in lbl]
    for name in ('TempPar', 'Default')
}

for network_name in list(networks_human.keys()):
    networks_human[f'{network_name} LH'] = [i for i, lbl in enumerate(labels_n) if f'LH_{network_name}' in lbl]
    networks_human[f'{network_name} RH'] = [i for i, lbl in enumerate(labels_n) if f'RH_{network_name}' in lbl]

networks_human['all'] = np.arange(1000, dtype=int)
networks_human['all LH'] = np.arange(500, dtype=int)
networks_human['all RH'] = np.arange(500, 1000, dtype=int)

# Macaque network parcellation (RM82)
lut_macaque = read_table('./dependencies/LUT_shen.txt')
tp_regions = ['secondary auditory cortex', 'inferior parietal cortex',
                'superior temporal cortex', 'central temporal cortex']
dmn_regions = ['medial prefrontal cortex', 'medial parietal cortex',
                'retrosplenial cingulate cortex', 'posterior cingulate cortex',
                'dorsolateral prefrontal cortex']
networks_macaque = {
    'TempPar': np.where(lut_macaque['name'].isin(tp_regions))[0],
    'Default': np.where(lut_macaque['name'].isin(dmn_regions))[0],
}

for network_name, region_list in [('TempPar', tp_regions), ('Default', dmn_regions)]:
    networks_macaque[f'{network_name} LH'] = np.where(
        (lut_macaque['name'].isin(region_list)) &
        (lut_macaque['hem'] == 'left')
    )[0]

    networks_macaque[f'{network_name} RH'] = np.where(
        (lut_macaque['name'].isin(region_list)) &
        (lut_macaque['hem'] == 'right')
    )[0]

networks_macaque['all'] = np.arange(len(lut_macaque), dtype=int)
networks_macaque['all LH'] = np.where(lut_macaque['hem'] == 'left')[0]
networks_macaque['all RH'] = np.where(lut_macaque['hem'] == 'right')[0]

# AAL 90 network parcellation

lut_aal = loadmat('./dependencies/LUT_AAL.mat')

labels_aal = [lbl.strip() for lbl in lut_aal['label90'].tolist()]

temporoparietalRegions = {
    'L Angular', 'R Angular',
    'L SupraMarginal', 'R SupraMarginal',
    'L Parietal Inf', 'R Parietal Inf',
    'L Temporal Sup', 'R Temporal Sup'
}

DMNRegions = {
    'L Front Sup Med', 'R Front Sup Med',
    'L Cingulum Post', 'R Cingulum Post',
    'L Precuneus', 'R Precuneus',
    'L Angular', 'R Angular'
}

networks_aal = {
    'TempPar': np.array(
        [i for i, lbl in enumerate(labels_aal) if lbl in temporoparietalRegions],
        dtype=int
    ),

    'Default': np.array(
        [i for i, lbl in enumerate(labels_aal) if lbl in DMNRegions],
        dtype=int
    ),
}

# Hemisphere subnetworks
for network_name, region_set in [('TempPar', temporoparietalRegions),
                                 ('Default', DMNRegions)]:

    networks_aal[f'{network_name} LH'] = np.array(
        [i for i, lbl in enumerate(labels_aal)
         if lbl in region_set and lbl.startswith('L ')],
        dtype=int
    )

    networks_aal[f'{network_name} RH'] = np.array(
        [i for i, lbl in enumerate(labels_aal)
         if lbl in region_set and lbl.startswith('R ')],
        dtype=int
    )

# Whole brain
networks_aal['all'] = np.arange(90, dtype=int)

networks_aal['all LH'] = np.array(
    [i for i, lbl in enumerate(labels_aal) if lbl.startswith('L ')],
    dtype=int
)

networks_aal['all RH'] = np.array(
    [i for i, lbl in enumerate(labels_aal) if lbl.startswith('R ')],
    dtype=int
)

def find_peak_frequency(signal, sampling_rate):
    signal = np.asarray(signal)
    signal = signal - np.mean(signal)

    N = len(signal)

    fft_values = np.abs(np.fft.fft(signal))
    fft_values = fft_values[:N // 2 + 1]

    frequencies = (sampling_rate / N) * np.arange(0, N // 2 + 1)

    peak_index = np.argmax(fft_values)
    peak_frequency = frequencies[peak_index]

    return peak_frequency

parcellations = {'Schaefer': [theta_human, networks_human], 'RM82': [theta_macaque, networks_macaque], 'AAL': [theta_aal, networks_aal]}
data = []

for parcel_name in parcellations:

    theta = parcellations[parcel_name][0]
    networks = parcellations[parcel_name][1]

    for network_name in ['TempPar', 'Default', 'all']:
        network_all = networks[network_name]
        network_lh = networks[f'{network_name} LH']
        network_rh = networks[f'{network_name} RH']

        lh_freq = 0
        for i in network_lh:
            lh_freq += find_peak_frequency(np.sin(theta[i,1000:]),100)
        lh_freq = lh_freq/len(network_lh)

        rh_freq = 0
        for i in network_rh:
            rh_freq += find_peak_frequency(np.sin(theta[i,1000:]),100)
        rh_freq = rh_freq/len(network_rh)

        r_network = np.abs(np.mean(np.exp(1j * theta[network_all,1000:]), axis=0))
        sco_freq = find_peak_frequency(r_network, 100)

        data_temp = {
            'Parcellation': parcel_name,
            'Network': network_name,
            'Sync. Freq. LH': np.round(lh_freq,3),
            'Sync. Freq. RH': np.round(rh_freq,3),
            'Predicted SCO Freq. (LH-RH)': np.round(np.abs(lh_freq - rh_freq),3),
            'Observed SCO Freq.': np.round(sco_freq,3)
            }

        data.append(data_temp)

df = pd.DataFrame(data)
df.to_csv('./figure_data/predicted_sco_frequencies_result.csv')