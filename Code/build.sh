#!/bin/bash
mkdir -p $PREFIX/bin
mkdir -p $PREFIX/share
mkdir -p $PREFIX/share/TransposonAnnotator_reasonaTE
cp $RECIPE_DIR/*.py $PREFIX/share/TransposonAnnotator_reasonaTE
cp -r $RECIPE_DIR/ncbi_cdd_candidates $PREFIX/share/TransposonAnnotator_reasonaTE
cp $RECIPE_DIR/reasonaTE $PREFIX/bin
chmod +x $PREFIX/bin/reasonaTE
