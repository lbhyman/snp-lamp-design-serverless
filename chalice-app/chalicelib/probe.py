import chalicelib.thermo_utils as tu

# Designates a fluorogenic probe set, including sequences and useful methods during GA and hill-climbing optimization

class Probe:
    # Initialize from dict object
    def __init__(self, obj):
        self.sequences = {'WT': obj['WT'], 'SNP': obj['SNP']}
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
        self.truncations = tu.generate_truncations(len(self.sequences['SNP']),self.minlength)
        self.id_to_sequence(self.truncations)
    
    def get_key(self):
        return ','.join(str(elem) for elem in self.truncations)+','+self.sequences['SNP']+','+self.sequences['WT']
    
    # Convert truncation representation to sequence representation    
    def id_to_sequence(self,seq_id):
        SNP_comp = tu.reverse_complement(self.sequences['SNP'])
        self.sequences['probeF'] = tu.truncate(tu.truncate(SNP_comp,5,seq_id[0]),3,seq_id[1])
        self.sequences['probeQ'] = tu.truncate(tu.truncate(self.sequences['SNP'],5,seq_id[2]),3,seq_id[3])
        WT_comp = tu.reverse_complement(self.sequences['WT'])
        self.sequences['sink'] = tu.truncate(tu.truncate(WT_comp,5,seq_id[5]),3,seq_id[6])
        self.sequences['sinkC'] = tu.truncate(tu.truncate(self.sequences['WT'],5,seq_id[7]),3,seq_id[8])