#!/bin/bash

for i in input/ComsolAges/revised_Tp_qx_qy_qz_*; do (echo ./comsolAges $i ${i/input/output} &); done

