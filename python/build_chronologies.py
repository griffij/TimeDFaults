"""Construct chronologies for a given fault
"""

from QuakeRates.dataman.parse_params import parse_param_file, get_event_sets

paramfiles = ['../params/Dunstan_GNS_unpub_simple.txt']
# Number of sample chronologies to generate
n_samples = 10

names, event_sets, event_certainties, num_events, \
    tect_regions, fault_styles = get_event_sets(paramfiles, ['all'], ['all'], 1)
print(names)
for i, event_set in enumerate(event_sets):
    event_set.gen_chronologies(n_samples, observation_end=2020, min_separation=1)
    chron_filename = '../data/chronologies/' + event_set.name + \
        '_%i_chronologies.csv' % n_samples
    event_set.write_chronology(chron_filename)
    
