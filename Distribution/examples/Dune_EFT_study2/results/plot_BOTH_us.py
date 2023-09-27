import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interp1d
import numpy as np

# Enabling LaTeX rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def convert_name_to_latex(name):
    last_part = name.split('_')[-1].lower()
    latex_mapping = {
        'ee': '\mathrm{ee}',
        'emu': '\mathrm{e} \mu',
        'etau': '\mathrm{e} \mathrm{\\tau}',
        'mue': '\mu \mathrm{e}',
        'mumu': '\mu \mu',
        'mutau': '\mu \mathrm{\\tau}',
        'taue': '\\tau \mathrm{e}',
        'taumu': '\\tau \mu',
        'tautau': '\\tau \\tau'
    }
    return r'$\varepsilon_P^{u d}\left\{' + latex_mapping[last_part] + r'\right\}$'

file_path = 'CC_P_us_new_resultsBOTH.txt' # Update the path if needed

interpolation_results = {}

with open(file_path, 'r') as file:
    current_name = None
    x_values = []
    y_values = []
    
    for line in file:
        if not line.strip():
            continue
            
        name, param, chi2 = line.strip().split()
        param, chi2 = float(param), float(chi2)

        if name != current_name and current_name is not None:
            interp_func = interp1d(y_values, x_values, kind='linear', bounds_error=False, fill_value='extrapolate')
            param_value_at_chi2_2_71 = interp_func(2.71)
            interpolation_results[current_name] = param_value_at_chi2_2_71
            x_values, y_values = [], []
            
        current_name = name
        x_values.append(param)
        y_values.append(chi2)

    interp_func = interp1d(y_values, x_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    param_value_at_chi2_2_71 = interp_func(2.71)
    interpolation_results[current_name] = param_value_at_chi2_2_71

names_latex = [convert_name_to_latex(name) for name in interpolation_results.keys()]
param_values = [value for value in interpolation_results.values()]

cmap = LinearSegmentedColormap.from_list('custom_cmap', ['lightblue', 'blue'])

plt.figure(figsize=(10, 6))
for idx, value in enumerate(param_values):
    plt.barh(names_latex[idx], 1e3 - value, left=value, color=cmap(idx / len(param_values)), log=True)

plt.xscale('log')
plt.xlim(1e3, 1e-5)

plt.title('')
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')
plt.grid(axis='x', linestyle='--')
plt.xticks(fontsize=16)
plt.xlabel(r'$\varepsilon$', fontsize=16)
plt.tight_layout()
pdf_path = 'CC_P_us_Both.pdf'
plt.savefig(pdf_path, format='pdf')

