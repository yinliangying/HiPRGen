if [ -d ./scratch ]; then
    rm -r ./scratch
fi

mkdir scratch

python run_bucketing.py ./data/ronald_LIBE.json ./scratch/ronald_mol_entries.pickle ./scratch/buckets.sqlite ./scratch/report.tex

mpiexec --use-hwthread-cpus -n 8 python run_network_generation.py ./scratch/ronald_mol_entries.pickle ./scratch/buckets.sqlite ./scratch/rn.sqlite ./scratch/report.tex

python run_initial_state.py ./scratch/rn.sqlite ./scratch/ronald_mol_entries.pickle

RNMC --database=./scratch/rn.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --logging=true

python test.py