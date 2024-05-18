# setup conda environment

# must pass name of environment
if [ -z "$*" ]; then echo "must pass name of environment"; fi

conda install -y jupyter
#conda install -y sos sos-pbs -c conda-forge
#conda install -y sos-notebook jupyterlab-sos sos-papermill -c conda-forge
#conda install -y sos-r sos-python sos-bash -c conda-forge
conda install -y matplotlib
conda install -y numpy
conda install -y scipy
conda install -y pandas
conda install -y bioconda::pysam
conda install -y seaborn
conda install -y python-levenshtein
pip3 install notebook
pip3 install ipykernel
python -m ipykernel install --user --name=$@
jupyter kernelspec list