import numpy as np
import random
output = "./instances/random.tsp"
num_targets = 30

r_min = 65.7
density = 0.1

square_side = r_min * np.sqrt(num_targets / density)

with open(output, 'w') as f:
    f.write("NODE_COORD_SECTION\n")
    for i in range(1, num_targets + 1):
        f.write(f"{i} {random.uniform(0, square_side)} {random.uniform(0,square_side)}\n")
    f.write("EOF")
    