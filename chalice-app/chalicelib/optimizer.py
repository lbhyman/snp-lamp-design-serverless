from chalicelib.probe import Probe
from chalicelib.LSH import LSH
import os

DEFAULT_PARAMS = {
    'WT': '',
    'SNP': '',
    'minlength': 6,
    'mut_rate': 0.5,
    'beta': [],
    'truncations': [],
    'params': {'temperature': 20.0, 
              'sodium': 0.05, 
              'magnesium': 0.008},
    'concentrations': {'non_mut_target' : 1e-7,
            'mut_target': 1e-7,
            'probeF' : 1e-7,
            'probeQ' : 1e-7,
            'sink' : 1e-7,
            'sinkC' : 1e-7} 
}

class Optimizer():
    
    def __init__(self):
        filename = os.path.join(os.path.dirname(__file__), 'projections.txt')
        self.LSH = LSH(4, 5, 60, filename)
    
    def process_input(self, item):
        WT, SNP = item['WT'].upper(), item['SNP'].upper()
        SNP_base, SNP_index = None, None
        for i in range(len(WT)):
            if item['SNP'][i] != WT[i]:
                SNP_base = item['SNP'][i]
                SNP_index = i
                break
        return {'WT': WT, 'SNP': SNP, 'SNP_index': SNP_index, 'SNP_base': SNP_base}
    
    def optimize(self, item):
        item = self.process_input(item)
        WT, SNP = item['WT'], item['SNP']
        results = self.LSH.get(item)
        if len(results) == 0:
            return None
        result = results[0]
        curr_params = DEFAULT_PARAMS
        curr_params['WT'] = WT
        curr_params['SNP'] = SNP
        truncs = [int(result['trunc_'+str(n)]) for n in range(1,10)]
        curr_params['truncations'] = truncs
        probe = Probe(curr_params)
        return probe.sequences