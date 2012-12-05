set -e
version=$( awk '/#define BAM_VERSION/ { print $3 }' bam.h | tr -d '"')
version="${version}-`git describe --always`"
echo "Building and installing samtools-${version}"
make
make install prefix=/home/public/usr/stow/samtools-$version
cd /home/public/usr/stow
stow -D samtools-*
stow -v samtools-$version
cd -
