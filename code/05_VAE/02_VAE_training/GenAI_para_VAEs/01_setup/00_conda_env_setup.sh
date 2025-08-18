conda create -n GenAI_VAE python=3.9 --yes
conda activate GenAI_VAE
conda install pytorch pytorch-cuda=11.7 -c pytorch -c nvidia --yes
pip install "ray[tune]"
pip install "ray[data]"
pip install "ray[air]"
conda install -c conda-forge packaging  --yes
pip install transformers
conda install -c pytorch torchtext --yes
pip install hpbandster ConfigSpace
conda install -c anaconda openpyxl --yes
conda install -c conda-forge matplotlib --yes
conda install -c anaconda seaborn  --yes
conda install -c plotly plotly=5.11.0  --yes
pip install chardet
pip install charset-normalizer==2.1.0