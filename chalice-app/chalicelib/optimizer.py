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
        WT = item['WT']
        SNP_base, SNP_index = None, None
        for i in range(len(WT)):
            if item['SNP'][i] != WT[i]:
                SNP_base = item['SNP'][i]
                SNP_index = i
                break
        return {'WT': WT, 'SNP_index': SNP_index, 'SNP_base': SNP_base}
    
    def optimize(self, item):
        item = self.process_input(item)
        WT = item['WT']
        SNP = WT[:item['SNP_index']] + item['SNP_base'] + WT[item['SNP_index']+1:]
        results = self.LSH.get(item)
        best_probe = None
        for i in range(len(results)):
            result = results[i]
            curr_params = DEFAULT_PARAMS
            curr_params['WT'] = WT
            curr_params['SNP'] = SNP
            truncs = [int(result['trunc_'+str(n)]) for n in range(1,10)]
            curr_params['truncations'] = truncs
            probe = Probe(curr_params)
            beta = probe.calc_beta()
            if best_probe == None:
                best_probe = probe
            elif probe.beta[0] > best_probe.beta[0]:
                best_probe = probe
            if beta[0] >= 1.5:
                return probe.sequences
            if i >= 9:
                return best_probe.sequences
        if best_probe != None:
            return best_probe.sequences
        return None