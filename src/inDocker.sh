#!/bin/bash
if [ -f /.dockerinit ]; then
    echo "True";
else
    echo "False";
fi
