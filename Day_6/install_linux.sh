# Download conda installation script
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh

# Run installation script
bash Anaconda3-2019.03-Linux-x86_64.sh

# Initialize conda
source ~/.bashrc

# Installa packages
conda install pybind11 cmake eigen
