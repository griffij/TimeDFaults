"""Construct chronologies for a given fault
"""

from QuakeRates.dataman.parse_params import parse_param_file, get_event_sets
from glob import glob
#paramfiles = ['../params/Dunstan4event_VanDissen2007_simple.txt',
#              '../params/Dunstan5event_VanDissen2007_simple.txt',
#              '../params/Dunstan6event_VanDissen2007_simple.txt']
#paramfiles = ['../params/Dunstan6eventOxcal_devonshire.csv',
#              '../params/Dunstan5eventOxcal_devonshire.csv',
#              '../params/Dunstan5eventOxcalv2_devonshire.csv',
#              '../params/Dunstan4eventOxcal_devonshire.csv']
#paramfiles = ['../params/Hyde3event_lugcreek.txt',
#              '../params/Hyde4event_lugcreek.txt']
paramfiles = ['../params/Hyde4event_lugcreek_v3.txt']
# Number of sample chronologies to generate
n_samples = 10000

names, event_sets, event_certainties, num_events, \
    tect_regions, fault_styles = get_event_sets(paramfiles, ['all'], ['all'], 1)
print(names)
for i, event_set in enumerate(event_sets):
    event_set.gen_chronologies(n_samples, observation_end=2020, min_separation=1)
    chron_filename = '../data/chronologies_all/' + event_set.name + \
        '_%i_chronologies.csv' % n_samples
    event_set.write_chronology(chron_filename)
    
