#!/usr/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

np="${1:-90}"
export xyz_dir="/home/victor/Data/MD_Data/AgCo/Dataset1/XYZ"
export data_file="/home/victor/Data/MD_Data/AgCo/Dataset1/data.dat"
export substrate_dir="/home/victor/Data/carbon_substrates"
export images_data="images.dat"
export descriptors_data="data.dat"
export max_files=-1 # Put -1 to use all files in the dataset
export atom_typ1="Ag"
export atom_typ2="Co"
export n_px=256
export nz=20
export ht=200.0



start_time=$(date +%s)

mpirun --use-hwthread-cpus -np "$np" ./main

end_time=$(date +%s)
echo "Execution time: $((end_time - start_time)) seconds"
