# test_three_processes
Nextflow workflow development

Usage:
`nextflow run main.nf -params-file inputs_local.json`

## Changes
- The apptainer directory contains the Dockerfile used to build Apptainer inside Docker.  You can replace the .sif with any singularity image.  Please contact rcc@odu.edu if you want a copy of bamtools.sif. Also, since I'm not sure if squashfs is installed on Omics/ECS, I included the unsquash command for compatibility.
- To run Apptainer inside Docker, I think you must run docker with the --privileged option.  This is added to nextflow.config. Otherwise, you will get an error about building apptainer with suid or allowing creation of user namespaces.  Again, this was added for compatibility since I don't know how Omics runs or if it can be customized.
- I'm not sure why but this container does not work with symlinks, which is the default behavior with NextFlow while staging input/output files in the working directory.  To workaround this, set stageInMode to copy in main.nf.
- I published the bamtools container to dockerhub. It can be found at https://hub.docker.com/r/tstilwel/apptainer-in-docker/tags
