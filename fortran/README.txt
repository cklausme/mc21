Tested on MacOS 13.4:

1) Install Xcode command-line tools <https://developer.apple.com/download/all/>
2) Install gfortran using Homebrew <https://docs.brew.sh/Installation>
3) To compile:
gfortran -o mc21lsodes -O3 -ftree-vectorize -ffixed-line-length-0 '-std=legacy' -Wl,-framework -Wl,Accelerate opkda1.f opkda2.f opkdmain.f mc21lsodes.f

Run using RunSimLSODEs function in mc21.nb Mathematica notebook.

