import numpy as np

class ShieldingTensorAnalyzer:
    '''
    inputs: 
    data_dict:output data from read_orca
    reference_eigvalues: reference compound isotropic shielding value
    '''
    def __init__(self, data_dict, reference_eigenvalues):
        self.data_dict = data_dict
        self.ref_eigenvalues = reference_eigenvalues
        self.modified_eigvalues = self.calculate_modified_eigvalues()
    
    def calculate_modified_eigvalues(self):
        modified_eigvalues = {}
        for nucleus, tensor in self.data_dict.items():
            eigenvalues = np.linalg.eigvals(tensor).real
            modified_eig = self.ref_eigenvalues-eigenvalues 
            modified_eigvalues[nucleus] = modified_eig
        return modified_eigvalues
    
    def eigvalues(self, nucleus):
        tensor = self.data_dict.get(nucleus)
        if tensor is not None:
            return self.modified_eigvalues.get(nucleus)
        else:
            return None
    
    def isotropic_shift(self, nucleus):
        eigenvalues = self.modified_eigvalues.get(nucleus)
        if eigenvalues is not None:
            sum_eigs = sum(eigenvalues)
            return sum_eigs / 3
        else:
            return None
    
    def sorted_eigvalues(self, nucleus):
        eigenvalues = self.modified_eigvalues.get(nucleus)
        if eigenvalues is not None:
            isotropic = self.isotropic_shift(nucleus)
            abs_diff = np.abs(eigenvalues - isotropic)
            sorted_indices = np.argsort(abs_diff)[::-1]  # Sort in descending order
            sorted_eigenvalues = eigenvalues[sorted_indices]
            return sorted_eigenvalues
        else:
            return None
    
    def delta_value(self, nucleus):
        sorted_eigenvalues = self.sorted_eigvalues(nucleus)
        isotropic = self.isotropic_shift(nucleus)
        if sorted_eigenvalues is not None and isotropic is not None:
            deltazz = sorted_eigenvalues[0] - isotropic
            return deltazz
        else:
            return None
    
    def ita_value(self, nucleus):
        sorted_eigenvalues = self.sorted_eigvalues(nucleus)
        delta = self.delta_value(nucleus)
        if sorted_eigenvalues is not None and delta is not None and len(sorted_eigenvalues) >= 3:
            deltaxx = sorted_eigenvalues[1]
            deltayy = sorted_eigenvalues[2]
            ita = (deltayy - deltaxx) / delta
            return ita
        else:
            return None


