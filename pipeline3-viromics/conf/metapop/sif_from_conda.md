to create a sif image from a conda environment

# Option 1

1- activate your conda env
```bash
conda activate your_env
conda env export > environment.yml
```

2- create a docker file
```bash
Bootstrap: docker

From: continuumio/miniconda3

%files
    environment.yml

%post
    /opt/conda/bin/conda env create -f environment.yml

%runscript
    exec /opt/conda/envs/$(head -n 1 environment.yml | cut -f 2 -d ' ')/bin/"$@"
```

3- Build the image
```
singularity build conda.sif Singularity
singularity build --force conda.sif docker_file_
```


# Option 2
1- Install `conda-pack`
```
conda-pack -n <MY_ENV> -o packed_environment.tar.gz
```

2- Create Docker file
```bash
Bootstrap: docker

From: continuumio/miniconda3

%files
    packed_environment.tar.gz /packed_environment.tar.gz

%post
    tar xvzf /packed_environment.tar.gz -C /opt/conda
    conda-unpack
    rm /packed_environment.tar.gz
```

3- Build the image
```bash
singularity build --fakeroot <OUTPUT_CONTAINER.sif> Singularity
```
