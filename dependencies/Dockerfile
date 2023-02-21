FROM nvidia/cuda:11.2.1-runtime-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN \
  # Update package list
  apt update -y && \
  # Install...
  apt install -y \
  git \
  vim \
  wget \
  python3 \
  python3-dev \
  python3-pip \
  tmux \
  # Needed for fragmenstein
  libxrender1 \
  # Remove package lists
  && rm -rf /var/lib/apt/lists/*

# # Install conda/mamba
# RUN wget -O Mambaforge.sh  "https://github.com/conda-forge/miniforge/releases/download/22.11.1-1/Mambaforge-22.11.1-1-$(uname)-$(uname -m).sh"
# RUN bash Mambaforge.sh -b -p $HOME/conda

# COPY cache_3_env.yml /tmp/cache_3_env.yml
# RUN $HOME/conda/bin/mamba env create -f /tmp/cache_3_env.yml

# RUN echo "export PATH=\$HOME/conda/bin:\$PATH" >> $HOME/.bashrc
# RUN echo "conda init bash" >> $HOME/.bashrc
# RUN echo "conda activate cache_3" >> $HOME/.bashrc

COPY installs /installs
RUN cd /installs/PyRosetta4.Release.python38.ubuntu.release-337/setup && \
    python3 setup.py install


RUN pip3 install \
      fragmenstein==0.9.12.6 \
      rdkit==2022.09.4 \
      biopython==1.80 \
      numpy==1.24.2 \
      pathos==0.3.0


ENTRYPOINT ["/bin/bash"]