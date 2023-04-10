INC_DIR=$1

GIT_VERSION=$(git describe --always --tags)

echo "#define MASHMAP_GIT_VERSION" \"$GIT_VERSION\" > "$INC_DIR"/mashmap_git_version.hpp
