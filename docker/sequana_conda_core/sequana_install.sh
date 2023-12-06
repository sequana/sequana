cd /home/sequana

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh -b -f -p /home/sequana/miniconda3

source /home/sequana/miniconda3/bin/activate

conda config --add channels r
conda config --add channels bioconda
conda install numpy matplotlib pandas scipy ipython graphviz -y


cp /home/sequana/miniconda3/lib/python3.5/site-packages/matplotlib/mpl-data/matplotlibrc .
sed -i -e 's/Qt4Agg/Agg/'g matplotlibrc
sed -i -e 's/Qt5Agg/Agg/'g matplotlibrc
