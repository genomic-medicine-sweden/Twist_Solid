#!/usr/bin/env bash
set -e


TAG_OR_BRANCH="${TAG_OR_BRANCH:-develop}"

echo "Build container ...."
docker build --build-arg TAG_OR_BRANCH=${TAG_OR_BRANCH} -f Twist_Solid.dockerfile . --tag gms/twist_solid:$TAG_OR_BRANCH --no-cache

echo "Export container to twist_solid_temp.tar.gz"
docker save gms/twist_solid:${TAG_OR_BRANCH} > twist_solid_temp.tar.gz

echo "Build singularirty..."
singularity build --fakeroot twist_solid_${TAG_OR_BRANCH}.sif Twist_Solid.def

echo "Remove twist_solid_temp.tar.gz"
rm twist_solid_temp.tar.gz
