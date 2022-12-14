#! /bin/bash
# [[file:trajectory.note::*install.sh][install.sh:1]]
version=v0.2.0
cargo im --offline

#cargo im
install -D -t bin/$version ~/.cargo/bin/topo
install -D -t bin/$version ~/.cargo/bin/lindemann
upx bin/$version/*
# install.sh:1 ends here
