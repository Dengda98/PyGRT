#!/bin/bash
set -e

# 运行教程中所有示例的脚本，方便将成图过程放在readthedocs服务器上


if [[ $1 == '1' || $1 == '' ]]; then
    cd Tutorial/dynamic/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi


if [[ $1 == '2' || $1 == '' ]]; then
    cd Tutorial/static/run
    chmod +x *.sh
    ./run.sh
    # python run.py
    cd -
fi

if [[ $1 == '3' || $1 == '' ]]; then
    cd Tutorial/dynamic/run_upar
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '4' || $1 == '' ]]; then
    cd Tutorial/static/run_upar
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '5' || $1 == '' ]]; then
    cd Advanced/integ_converg/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '6' || $1 == '' ]]; then
    cd Advanced/kernel/run
    chmod +x *.sh
    ./run.sh
    python run.py
    cd -
fi

if [[ $1 == '7' || $1 == '' ]]; then
    cd Advanced/filon/run
    chmod +x *.sh
    ./run.sh
    python run.py
    python plot.py
    cd -
fi

if [[ $1 == '8' || $1 == '' ]]; then
    cd Lamb_problem/run
    chmod +x *.sh
    ./run.sh
    python lamb1_plot_time.py
    python lamb1_plot_freq_time.py
    cd -
fi

if [[ $1 == '9' || $1 == '' ]]; then
    cd Advanced/waveform_drift/run
    python run.py
    cd -
fi