# https://github.com/anaconda/docker-images
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /src

# Install unzip
RUN apt-get update && apt-get install -y unzip

# Create and activate the environment
COPY lhasa.yml /src/
RUN conda env create --file lhasa.yml

RUN echo "source activate lhasa" > ~/.bashrc
ENV PATH /opt/conda/envs/lhasa/bin:$PATH

# Download and unzip static.zip
RUN wget https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip -O /tmp/static.zip && \
    unzip /tmp/static.zip && \
    rm /tmp/static.zip

# Download and unzip exposure.zip
RUN wget https://gpm.nasa.gov/sites/default/files/data/landslides/exposure.zip -O /tmp/exposure.zip && \
    unzip /tmp/exposure.zip && \
    rm /tmp/exposure.zip

# Download and unzip ref_data.zip
RUN wget https://gpm.nasa.gov/sites/default/files/data/landslides/ref_data.zip -O /tmp/ref_data.zip && \
    unzip /tmp/ref_data.zip -d pfdf && \
    rm /tmp/ref_data.zip

# Necessary directories/files
RUN mkdir imerg && mkdir smap
RUN touch ~/.urs_cookies && touch ~/.dodsrc
RUN echo "HTTP.NETRC=~/.netrc" >> ~/.dodsrc && \
    echo "HTTP.COOKIEJAR=~/.urs_cookies" >> ~/.dodsrc

# Copy LHASA model and script
COPY lhasa.py model.json /src/

# Set the entrypoint; write all files to the output directory (volume)
ENTRYPOINT ["python", "lhasa.py", "--output_path=output"]
