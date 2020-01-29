build_dir="build"

mkdir -p $build_dir
cd $build_dir

python3 -m http.server 8000 --bind localhost
