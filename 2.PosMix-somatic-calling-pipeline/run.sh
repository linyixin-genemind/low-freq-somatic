#!/bin/bash
Sopt="-B /home:/home -B /mnt:/mnt -B /mnt/nfsp1/gitlab/biobase:/biobase"
sif="/mnt/nfsp1/gitlab/biobase/image/biobase.sif"
env="/mnt/nfsp1/gitlab/biobase/apps/somatic_snv_posmix/env"

singularity exec $Sopt $sif mamba run -p $env $@