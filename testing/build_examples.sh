#!/bin/bash
cd ../example/build
cmake ..
make
cp natcsi ../../testing/
cd ../../example2/build
cmake ..
make
cp legendre ../../testing/
cd ../../testing
