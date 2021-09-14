FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install git curl wget zsh -y
RUN chsh -s $(which zsh)
SHELL ["/bin/zsh", "-c"]
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"
WORKDIR /root
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/usr/local/miniconda/bin:$PATH
RUN conda init zsh
RUN conda update conda -y
RUN conda update --all -y
RUN apt-get install build-essential clustalw probcons libboost-all-dev pkg-config bzip2 -y
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install numpy matplotlib pandas seaborn scipy biopython locarna mafft contrafold rnastructure -y
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN source /root/.cargo/env
ENV PATH=/root/.cargo/bin:$PATH
# If the following command fails, please remove RUSTFLAGS='...' to disable advanced installation options
RUN RUSTFLAGS='--emit asm -C target-feature=+avx -C target-feature=+ssse3 -C target-feature=+mmx' cargo install consprob consalifold
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz
RUN tar -xf ViennaRNA-2.4.18.tar.gz
WORKDIR /root/ViennaRNA-2.4.18
RUN ./configure && make -j$(nproc) && make install
WORKDIR /root
RUN rm -rf ViennaRNA-2.4.18
RUN rm ViennaRNA-2.4.18.tar.gz
RUN wget https://github.com/satoken/centroid-rna-package/archive/refs/tags/v0.0.16.tar.gz
RUN tar -xf v0.0.16.tar.gz
WORKDIR /root/centroid-rna-package-0.0.16
RUN ./configure && make -j$(nproc) && make install
WORKDIR /root
RUN rm -rf centroid-rna-package-0.0.16
RUN rm v0.0.16.tar.gz
RUN wget https://rth.dk/resources/petfold/download/PETfold2.1.tar.gz
RUN tar -xf PETfold2.1.tar.gz
WORKDIR /root/PETfold/src
RUN make clean && make -j$(nproc)
RUN cp ../bin/* /usr/local/bin
RUN echo -e "export PETFOLDBIN=/usr/local/bin" >> /root/.zshrc
WORKDIR /root
RUN rm -rf PETfold
RUN rm PETfold2.1.tar.gz
RUN git clone https://github.com/heartsh/raf-fixed
WORKDIR /root/raf-fixed
RUN make -j$(nproc)
RUN cp raf /usr/local/bin
WORKDIR /root
RUN rm -rf raf-fixed
RUN git clone https://github.com/heartsh/contralign-fixed
WORKDIR /root/contralign-fixed
RUN make -j$(nproc)
RUN cp contralign /usr/local/bin
WORKDIR /root
RUN rm -rf contralign-fixed
RUN echo -e "\nexport CONTRAFOLD_DIR=/usr/local/miniconda/bin/contrafold" >> /root/.zshrc
RUN echo -e "\nexport CONTRALIGN_DIR=/usr/local/bin" >> /root/.zshrc
RUN apt-get clean -y && apt-get autoremove -y
RUN conda clean --all -y
