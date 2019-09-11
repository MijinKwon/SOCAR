from extract_ACGs import *
from extract_RMGs import *
from run_RWR import *
from evaluate_performance import *

### Get resistance module genes (RMGs) ###
extract_ACGs()
extract_RMGs()

### Simulate RWR for test drugs ###
process_rwr_input_file()
# os.system("Rscript random_walk_with_restart.R")

### Calculate potential sensitizer scores ###
calculate_PSS('DEGACG')
calculate_PSS('DEG')

### Evaluate model performance ###
calculate_AUROC('DEGACG')
calculate_AUROC('DEG')

### Predict sensitization mechanisms of candidate sensitizers ###
predict_sensitization_mechanisms()