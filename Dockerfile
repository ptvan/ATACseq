FROM ubuntu:24.04
WORKDIR /configs

SHELL ["bash", "-l" ,"-c"]
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN apt update && apt install -y wget curl bzip2 && \
    wget --no-check-certificate  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh 
RUN bash miniforge.sh -b -p /root/miniforge
ENV PATH="$PATH:/root/miniforge/bin"
COPY environment.yaml /configs
COPY *.nf /pipeline/
COPY nextflow.config /configs
RUN mamba env create -f /configs/environment.yaml && \
    mamba init
RUN source /root/.bashrc 
RUN echo "mamba activate atacseq" >> ~/.bashrc
ENTRYPOINT ["bash", "-l", "-c"]