input_params = {
    'V': 1e-6,
    'n': 7,
    'p': 3/float(4),
    'y': 1/float(3), 
    'a_in': 0.99, 
    'a_out': 0.1, 
    'theta': {
        'mean': 0.075, # [0.05, 0.2) % CHANGE?
        'variance': 0.01 # [0, 0.02) % CHANGE?    
    },
    'theta_in': {
        'mean': 0.0003, # Sm_base, # [0, 0.0006) % CHANGE?
        'variance': 0
    },
    'theta_out': {
        'mean': 0.0003, # [0, 0.0006) % CHANGE?
        'variance': 0
    },
    'phi': {
        'mean': 3, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': 0.001 # [0, 0.5)
    },
}