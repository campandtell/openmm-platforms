# Use NVIDIA CUDA base image
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    git \
    wget \
    eog \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------
# Install Miniconda (so we can get a GPU-enabled OpenMM via conda)
# ------------------------------------------------------------------
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Use bash for subsequent RUN commands
SHELL ["/bin/bash", "-c"]

# Create a new conda environment with Python 3.10
# Install GPU-enabled OpenMM plus any other Python packages
RUN conda create -y -n openmm-gpu python=3.10 \
    && conda run -n openmm-gpu conda install -y -c conda-forge \
       openmm cudatoolkit=11.8 \
       numpy matplotlib mdtraj

# (Optional) Automatically activate openmm-gpu env on container start
RUN echo "source activate openmm-gpu" >> ~/.bashrc

# Copy project files
COPY . /workspace/

# By default, launch a bash shell
CMD ["/bin/bash"]
