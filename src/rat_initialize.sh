#! /bin/bash

INIT=0
RUN=1

DEV=1


if [ $INIT -eq 1 ]
    if [ $DEV -eq 1 ]
    then
        # use locally installed rat packagels
        python models/RAT/src/rat/cli/rat_cli.py init -d ./ -g
    else
        rat init -d ./ -g
    fi
fi

if [ $RUN -eq 1 ]
    if [ $DEV -eq 1 ]
    then
        # use locally installed rat packagels
        python models/RAT/src/rat/cli/rat_cli.py run -d ./ -g
    else
        rat run -d ./ -g
    fi
fi

