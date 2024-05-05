#!/bin/sh

# assumes mamba installation

if command -v mamba &>/dev/null; then
    mamba env update -p $CONDA_PREFIX -f abcenv.yml
else
    conda env update -p $CONDA_PREFIX -f abcenv.yml
fi

$

wget -P ./bin/ https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
chmod +x ./bin/juicer_tools_1.22.01.jar

mkdir -p $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/envs/abc_lite/etc/conda/deactivate.d

echo "#!/bin/sh" >> $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d/env_vars.sh
chmod +x $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d/env_vars.sh
echo "export JUICERTOOLS=$(pwd)/bin/juicer_tools_1.22.01.jar" >> $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d/env_vars.sh
