import numpy as np
import random as rnd
import functools
import pickle as pk
import thermo_utils as tu

# Designates a fluorogenic probe set, including sequences and useful methods during GA and hill-climbing optimization
@functools.total_ordering
class Probe:
    # Initialize from dict object
    def __init__(self, obj):
        self.sequences = {'non_mut_target': obj['WT'], 'mut_target': obj['SNP']}
        self.minlength = obj['minlength']
        self.concentrations = obj['concentrations']
        self.params = obj['params']
        self.mut_rate = obj['mut_rate']
        self.beta = obj['beta']
        if len(obj['truncations']) == 0:
            self.generate_truncations()
        else:
            self.truncations = obj['truncations']
            self.id_to_sequence(self.truncations)
    
    # Generate random truncations from termini if no truncations are provided    
    def generate_truncations(self):
        self.truncations = tu.generate_truncations(len(self.sequences['mut_target']),self.minlength)
        self.id_to_sequence(self.truncations)
        while (not self.screen()):
            self.truncations = tu.generate_truncations(len(self.sequences['mut_target']),self.minlength)
            self.id_to_sequence(self.truncations)
    
    # Fitness-based sorting method for GA        
    def __lt__ (self, other):
        return self.beta[0] < other.beta[0]

    # Fitness-based sorting method for GA 
    def __eq__ (self, other):
        return self.beta[0] == other.beta[0]
    
    def get_sequences(self):
        return self.sequences
    
    def get_truncations(self):
        return self.truncations
    
    def get_key(self):
        return ','.join(str(elem) for elem in self.truncations)+','+self.sequences['mut_target']+','+self.sequences['non_mut_target']
    
    def get_beta(self):
        return self.beta
    
    def set_beta(self, b):
        self.beta = b
    
    # Convert truncation representation to sequence representation    
    def id_to_sequence(self,seq_id):
        SNP_comp = tu.reverse_complement(self.sequences['mut_target'])
        self.sequences['probeF'] = tu.truncate(tu.truncate(SNP_comp,5,seq_id[0]),3,seq_id[1])
        self.sequences['probeQ'] = tu.truncate(tu.truncate(self.sequences['mut_target'],5,seq_id[2]),3,seq_id[3])
        WT_comp = tu.reverse_complement(self.sequences['non_mut_target'])
        self.sequences['sink'] = tu.truncate(tu.truncate(WT_comp,5,seq_id[5]),3,seq_id[6])
        self.sequences['sinkC'] = tu.truncate(tu.truncate(self.sequences['non_mut_target'],5,seq_id[7]),3,seq_id[8])
    
    # Calculate probe fitness    
    def calc_beta(self):
        try:
            background, SNP, WT = tu.get_activation(self.sequences, self.concentrations, self.params)
            b = SNP * np.log10(SNP/max(WT,background,0.0001))
            self.beta = [b, SNP, WT, background]
            return [b, SNP, WT, background]
        except:
            self.beta = [0,0,0,0]
            return [0,0,0,0]

    # Randomly change one truncation when probe is mutated during GA
    def mutate(self):
        if rnd.random() < self.mut_rate:
            index = rnd.choice([0,1,2,3,5,6,7,8])
            blunt_end = self.truncations[4]
            change = rnd.choice([-1,1])
            if blunt_end == 5 and (index == 0 or index == 3):
                self.truncations[0] = self.truncations[0] + change
                self.truncations[3] = self.truncations[3] + change
            elif blunt_end == 3 and (index == 1 or index == 2):
                self.truncations[1] = self.truncations[1] + change
                self.truncations[2] = self.truncations[2] + change
            else:
                self.truncations[index] = self.truncations[index] + change
            for i in range(len(self.truncations)):
                if self.truncations[i] < 0:
                    self.truncations[i] = 0
            self.id_to_sequence(self.truncations)
    
    # Cross the probe with another probe during GA    
    def cross(self,other_probe):
        probe_options = [self.truncations[:5],other_probe.truncations[:5]]
        sink_options = [self.truncations[5:],other_probe.truncations[5:]]
        choice = rnd.choice([0,1])
        new_truncations = probe_options[choice] + sink_options[1-choice]
        new_probe_params = self.__dict__
        new_probe_params['truncations'] = new_truncations
        new_probe_params['beta'] = [0,0,0,0]
        child = Probe(new_probe_params)
        child.mutate()
        return child
    
    # Determine possible next steps for hill climbing
    def next_iteration(self):
        output = []
        truncations = self.truncations
        for i in range(len(truncations)):
            new_truncations = truncations.copy()
            new_truncations[i] = new_truncations[i] + 1
            if not new_truncations in output:
                output.append(new_truncations)
            new_truncations = truncations.copy()
            if truncations[i] >= 1:
                new_truncations[i] = new_truncations[i] - 1
                if not new_truncations in output:
                    output.append(new_truncations)
        output_probes = []
        for trunc in output:
            new_probe_params = self.__dict__
            new_probe_params['truncations'] = trunc
            new_probe_params['beta'] = [0,0,0,0]
            output_probes.append(Probe(new_probe_params).__dict__)
        return output_probes
    
    # TODO: Add pkl file
    # Determine whether the probe falls within the high-fitness design space by PCA
    def screen(self):
        probe = []
        SNP, WT, probe_seq, probe_duplex, sink_seq, sink_duplex = self.sequences['mut_target'], self.sequences['non_mut_target'], \
            self.sequences['probeF'], self.sequences['probeQ'], self.sequences['sink'], self.sequences['sinkC']
        tokenized = self.truncations
        # Add Tms
        probe = probe + [tu.TM(probe_duplex), tu.TM(probe_seq), tu.TM(sink_duplex), tu.TM(sink_seq)]
        # Add Tm ratios
        probe = probe + [tu.TM(probe_duplex)-tu.TM(sink_duplex), tu.TM(probe_seq)-tu.TM(sink_seq), \
            tu.TM(probe_seq)-tu.TM(probe_duplex), tu.TM(sink_seq)-tu.TM(sink_duplex)]
        # Add SNP proximity to ends
        SNP_pos = 0
        for i in range(len(SNP)):
            if SNP[i] != WT[i]:
                SNP_pos = i
        probe_3_dist = len(SNP) - 1 - tokenized[0] - SNP_pos
        probe_5_dist = SNP_pos - tokenized[1]
        probe_SNP_dist = min(probe_3_dist,probe_5_dist)
        probeC_3_dist = len(SNP) - 1 - tokenized[3] - SNP_pos
        probeC_5_dist = SNP_pos - tokenized[2]
        probeC_SNP_dist = min(probeC_3_dist,probeC_5_dist)
        sink_3_dist = len(SNP) - 1 - tokenized[0] - SNP_pos
        sink_5_dist = SNP_pos - tokenized[1]
        sink_SNP_dist = min(sink_3_dist,sink_5_dist)
        sinkC_3_dist = len(SNP) - 1 - tokenized[3] - SNP_pos
        sinkC_5_dist = SNP_pos - tokenized[2]
        sinkC_SNP_dist = min(sinkC_3_dist,sinkC_5_dist)
        probe = probe + [probe_SNP_dist, probeC_SNP_dist, sink_SNP_dist, sinkC_SNP_dist]
        # Test if the probe falls within the high-fitness PCA region
        pca = pk.load(open("../data/pca.pkl",'rb'))
        probe = np.array(probe)
        params_trans = pca.transform(probe.reshape(1,-1))
        if (params_trans[0][0] < 15 and params_trans[0][1] > -20 and params_trans[0][0] > -60 and \
            params_trans[0][1] < 60):
            return True
        return False