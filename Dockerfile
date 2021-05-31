# The pupose of this image is mostly to show how to setup the environment
# for FEniCS_ii. You can use the docker image either as 
# (build you own local one)
#
#  `docker build --no-cache -t xii-minimal .`
#
# and (run the container environment)
#
#  `docker run -it -v $(pwd):/home/fenics/shared xii-minimal`
#

# Base image
FROM quay.io/fenicsproject/dev

USER fenics

# fenics_ii dependency 
RUN git clone https://github.com/nschloe/quadpy.git && \
    cd quadpy && \
    git checkout v0.12.10 && \
    python3 setup.py install --user && \
    cd ..

# fenics_ii dependency 
RUN git clone https://mirok-w-simula@bitbucket.org/mirok-w-simula/cbc.block.git && \
    cd cbc.block && \
    python3 setup.py install --user && \
    cd ..


# fenics_ii 
RUN git clone https://github.com/MiroK/fenics_ii.git && \
    cd fenics_ii && \
    git fetch --all  && \
    git checkout hsfrac-minimal       && \
    cd ..
    
ENV PYTHONPATH="/home/fenics/fenics_ii/":"$PYTHONPATH"

USER root

# Release
# Might need to do docker login 
# docker build --no-cache -t xii-minimal .
# docker tag xii-minimal:latest mirok/xii-minimal
# docker push mirok/xii-minimal