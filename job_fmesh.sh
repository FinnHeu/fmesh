#!/bin/bash
#SBATCH --account=oze.oze
#SBATCH --job-name=modify_DARS
#SBATCH --partition=mpp
#SBATCH --time=06:00:00
#SBATCH --qos=12h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Load necessary modules
module load conda

# Activate the conda environment
conda activate fmesh

rm ./jigsaw/*
rm ./log/*.log
rm ./temp/*.pkl
rm ./mesh_files/*

rm *.pkl
rm *.vtk
rm *.jpg
rm *.out


echo "Mesh generation started."

# Run the mesh generation script
python fmesh.py > "./log/mesh_generation.log" 2>&1 &

# Wait for all background jobs to complete
wait

echo "Mesh generation completed."