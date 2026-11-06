#!/bin/bash
# dvrun.sh TOOL ARGS...: run a DeepVariant 1.10.0 tool (make_examples,
# call_variants, postprocess_variants) from the image unpacked in
# ~/vc/dv/rootfs, with the image's own loader, libraries and Python (a2 allows
# no containers). Paths given as arguments are host paths. Each tool's zip is
# unpacked once into ~/vc/dv/x/TOOL (dv_unzip.py: symlinks kept, as unzip
# does); its main module runs with the Python path the Bazel launcher would
# set (the runfiles root, its repositories and com_google_protobuf/python),
# since the launcher always unpacks itself again.
R=~/vc/dv/rootfs
tool=$1; shift
X=~/vc/dv/x/$tool
if [ ! -e $X/.done ]; then
  rm -rf $X; mkdir -p $X
  python3 -I ~/vc/scripts/dv_unzip.py $R/opt/deepvariant/bin/$tool.zip $X > /dev/null || exit 1
  touch $X/.done
fi
M=$X/runfiles
PP=$M:$M/com_google_protobuf/python
for d in $M/*/; do PP=$PP:${d%/}; done
PP=$PP:$R/usr/local/lib/python3.10/dist-packages:$R/usr/lib/python3/dist-packages
L=$R/usr/lib/x86_64-linux-gnu:$R/lib/x86_64-linux-gnu:$R/usr/local/lib:$R/usr/local/lib/python3.10/dist-packages/tensorflow
exec env -i HOME=$HOME PATH=$R/opt/conda/envs/bio/bin:/usr/bin:/bin TF_USE_LEGACY_KERAS=1 \
  PYTHONHOME=$R/usr PYTHONNOUSERSITE=1 PYTHONPATH=$PP RUNFILES_DIR=$M \
  OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} TMPDIR=${TMPDIR:-$HOME/vc/tmp} \
  $R/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2 --library-path $L \
  $R/usr/bin/python3.10 $M/com_google_deepvariant/deepvariant/$tool.py "$@"
