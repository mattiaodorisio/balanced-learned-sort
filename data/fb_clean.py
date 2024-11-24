# This removes the huge values at the end of facebook dataset
import numpy as np
import struct

DATASET_NAME = 'fb_200M_uint64'

data = np.fromfile(DATASET_NAME, dtype=np.uint64)
print("FB data read")

data = data[1:]
data = data[(data < np.quantile(data, 0.99999))]

with open(DATASET_NAME + "_cleaned", "wb") as f:
    f.write(struct.pack("Q", len(data)))
    data.tofile(f)

print(f"The FB dataset has been cleaned and now contains {len(data)} items")
