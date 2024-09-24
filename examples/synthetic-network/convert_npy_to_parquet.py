import numpy

# read in all .npy files in the subdirectory data
a_i = numpy.load('data/A_i.npy', allow_pickle=True)
a_p = numpy.load('data/A_p.npy', allow_pickle=True)
a_c = numpy.load('data/A_c.npy', allow_pickle=True)
l_i = numpy.load('data/L_i.npy', allow_pickle=True)
q_c = numpy.load('data/Q_c.npy', allow_pickle=True)
position = numpy.load('data/rel_positions.npy', allow_pickle=True)

# save as parquet files
import pandas as pd

pd.DataFrame(a_i).to_parquet('data/A_i.parquet')
pd.DataFrame(a_p).to_parquet('data/A_p.parquet')
pd.DataFrame(a_c).to_parquet('data/A_c.parquet')
pd.DataFrame(l_i).to_parquet('data/L_i.parquet')
pd.DataFrame(q_c).T.to_parquet('data/Q_c.parquet')
pd.DataFrame(position).to_parquet('data/rel_positions.parquet')
