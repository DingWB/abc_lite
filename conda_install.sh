#!/bin/sh

# assumes mamba installation

if [[ $(which mamba) ]]; then
    mamba env update -p $CONDA_PREFIX -f abcenv.yml
else
    conda env update -p $CONDA_PREFIX -f abcenv.yml
fi
wget -P ./bin/ https://github.com/aidenlab/JuicerTools/releases/download/v3.0.0/juicer_tools.3.0.0.jar
chmod +x ./bin/juicer_tools.3.0.0.jar

mkdir -p $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/envs/abc_lite/etc/conda/deactivate.d

echo "#!/bin/sh" >> $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d/env_vars.sh
echo "export JUICERTOOLS=./bin/juicer_tools.3.0.0.jar" >> $CONDA_PREFIX/envs/abc_lite/etc/conda/activate.d/env_vars.sh
