#!/bin/bash

# 运行教程中所有示例的脚本，方便将成图过程放在readthedocs服务器上


if [[ $1 == '1' || $1 == '' ]]; then
    cd dynamic/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi


if [[ $1 == '2' || $1 == '' ]]; then
    cd static/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '3' || $1 == '' ]]; then
    cd dynamic/run_upar
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '4' || $1 == '' ]]; then
    cd static/run_upar
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '5' || $1 == '' ]]; then
    cd integ_converg/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi