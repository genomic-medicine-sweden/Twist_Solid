
BRANCH_OR_TAG="add-validation-ref-yaml" bash build/build_conda.sh

TAG_OR_BRANCH="add-validation-ref-yaml" bash build_container.sh

#docker build --build-arg TwIST_SOLID_VERSION=develop . -f docker/Twist_DNA_Solid.dockerfile --tag twist_solid:develop

#singularity build twist_solid_develop.sif docker-daemon://twist_solid:develop
